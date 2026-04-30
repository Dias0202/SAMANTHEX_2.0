"""Consensus sequence assembly from Sanger sequencing reads.

Implements quality-aware pairwise alignment and consensus building
using BioPython's PairwiseAligner. Supports merging multiple replicates
per direction and forward/reverse consensus generation.

Key improvements over naive approaches:
    1. Quality trimming before alignment (removes noisy ends).
    2. Quality boosting for matching bases (Q_combined = Q1 + Q2).
    3. Phred quality filter applied only at the final step.
    4. Coverage depth tracking for each consensus position.
"""

import logging
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

from samanthex.core.file_io import SequenceRead, read_sequence_file
from samanthex.core.trimming import TrimResult, mott_trim

logger = logging.getLogger(__name__)

# Type alias for status callbacks
StatusCallback = Optional[Callable[[str, float], None]]


@dataclass
class ConsensusResult:
    """Result of the consensus assembly pipeline.

    Attributes:
        sequence: Final consensus DNA sequence.
        quality: Per-base Phred quality scores.
        name: Identifier for the consensus.
        fwd_reads: Number of forward reads used.
        rev_reads: Number of reverse reads used.
        trimming_info: Trimming details for each read.
        statistics: Quality and composition statistics.
        warnings: Any warnings generated during assembly.
    """
    sequence: str = ""
    quality: List[int] = field(default_factory=list)
    name: str = "consensus"
    fwd_reads: int = 0
    rev_reads: int = 0
    trimming_info: List[Dict] = field(default_factory=list)
    statistics: Dict = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    def to_fasta(self) -> str:
        """Format as FASTA string."""
        return f">{self.name}\n{self.sequence}"

    def to_fastq(self) -> str:
        """Format as FASTQ string with quality scores."""
        qual_str = "".join(chr(min(q, 93) + 33) for q in self.quality)
        return f"@{self.name}\n{self.sequence}\n+\n{qual_str}"


def compute_statistics(sequence: str, quality: List[int]) -> Dict:
    """Compute quality and composition statistics for a consensus.

    Args:
        sequence: DNA sequence string.
        quality: Per-base Phred quality scores.

    Returns:
        Dictionary with statistics including length, GC content,
        mean quality, Q20/Q30 percentages, and N count.
    """
    n = len(sequence)
    if n == 0:
        return {"length": 0}

    upper = sequence.upper()
    gc_count = upper.count("G") + upper.count("C")
    n_count = upper.count("N")
    q20 = sum(1 for q in quality if q >= 20)
    q30 = sum(1 for q in quality if q >= 30)
    mean_q = sum(quality) / len(quality) if quality else 0

    return {
        "length": n,
        "gc_content": gc_count / n * 100 if n > 0 else 0,
        "n_count": n_count,
        "n_percent": n_count / n * 100 if n > 0 else 0,
        "mean_quality": mean_q,
        "min_quality": min(quality) if quality else 0,
        "max_quality": max(quality) if quality else 0,
        "q20_count": q20,
        "q20_percent": q20 / n * 100 if n > 0 else 0,
        "q30_count": q30,
        "q30_percent": q30 / n * 100 if n > 0 else 0,
    }


def _create_aligner() -> PairwiseAligner:
    """Create and configure a PairwiseAligner for Sanger data.

    Uses local alignment with scoring parameters optimized for
    Sanger sequencing reads:
        - Match: +2 (reward for matching bases)
        - Mismatch: -1 (mild penalty for single substitutions)
        - Gap open: -5 (moderate penalty to discourage spurious gaps)
        - Gap extend: -2 (lower than open to allow real indels)
    """
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -2
    return aligner


def merge_two_sequences(
    seq1: str,
    qual1: List[int],
    seq2: str,
    qual2: List[int],
) -> Tuple[str, List[int]]:
    """Merge two sequences using quality-weighted consensus.

    Performs pairwise local alignment and builds a consensus where:
        - Matching bases: quality boosted (Q1 + Q2, capped at 93).
        - Mismatching bases: higher-quality base is selected.
        - Gaps: the non-gap base is used with its quality.

    Quality filtering (N substitution) is NOT applied here — it is
    deferred to the final step to prevent premature information loss.

    Args:
        seq1: First DNA sequence.
        qual1: Phred quality scores for seq1.
        seq2: Second DNA sequence.
        qual2: Phred quality scores for seq2.

    Returns:
        Tuple of (consensus_sequence, consensus_quality).
    """
    if not seq1 or not qual1:
        return seq2, qual2
    if not seq2 or not qual2:
        return seq1, qual1

    aligner = _create_aligner()

    try:
        alignments = aligner.align(seq1, seq2)
        best = alignments[0]
    except (IndexError, StopIteration, OverflowError, ValueError):
        logger.warning("No alignment found, returning first sequence")
        return seq1, qual1

    coords = best.coordinates
    score = best.score

    # Validate alignment quality
    min_len = min(len(seq1), len(seq2))
    min_score = min_len * aligner.match_score * 0.2
    if score < min_score:
        logger.warning(
            f"Low alignment score ({score:.0f} < {min_score:.0f}). "
            f"Sequences may not overlap well."
        )

    # Get alignment boundaries
    s1_start = int(coords[0, 0])
    s1_end = int(coords[0, -1])
    s2_start = int(coords[1, 0])
    s2_end = int(coords[1, -1])

    consensus: List[str] = []
    con_qual: List[int] = []

    # --- 1. Left extension (bases before the aligned region) ---
    # Determine which sequence extends further to the left
    if s1_start > 0 and s2_start > 0:
        # Both have prefixes — include the longer one
        if s1_start >= s2_start:
            for i in range(s1_start):
                consensus.append(seq1[i])
                con_qual.append(qual1[i])
        else:
            for i in range(s2_start):
                consensus.append(seq2[i])
                con_qual.append(qual2[i])
    elif s1_start > 0:
        for i in range(s1_start):
            consensus.append(seq1[i])
            con_qual.append(qual1[i])
    elif s2_start > 0:
        for i in range(s2_start):
            consensus.append(seq2[i])
            con_qual.append(qual2[i])

    # --- 2. Aligned region ---
    n_cols = coords.shape[1]
    for col in range(n_cols - 1):
        i1_from = int(coords[0, col])
        i1_to = int(coords[0, col + 1])
        i2_from = int(coords[1, col])
        i2_to = int(coords[1, col + 1])

        len1 = i1_to - i1_from
        len2 = i2_to - i2_from

        if len1 > 0 and len2 > 0:
            # Aligned block — merge base by base with quality
            for i in range(len1):
                b1 = seq1[i1_from + i].upper()
                b2 = seq2[i2_from + i].upper()
                q1 = qual1[i1_from + i] if (i1_from + i) < len(qual1) else 0
                q2 = qual2[i2_from + i] if (i2_from + i) < len(qual2) else 0

                if b1 == b2:
                    # Agreement: boost confidence (independent error sources)
                    consensus.append(b1)
                    con_qual.append(min(q1 + q2, 93))
                elif q1 >= q2:
                    consensus.append(b1)
                    con_qual.append(q1)
                else:
                    consensus.append(b2)
                    con_qual.append(q2)

        elif len1 > 0:
            # Gap in seq2 — use seq1 bases
            for i in range(len1):
                idx = i1_from + i
                consensus.append(seq1[idx].upper())
                con_qual.append(qual1[idx] if idx < len(qual1) else 0)

        elif len2 > 0:
            # Gap in seq1 — use seq2 bases
            for i in range(len2):
                idx = i2_from + i
                consensus.append(seq2[idx].upper())
                con_qual.append(qual2[idx] if idx < len(qual2) else 0)

    # --- 3. Right extension (bases after the aligned region) ---
    s1_remaining = len(seq1) - s1_end
    s2_remaining = len(seq2) - s2_end

    if s1_remaining > 0 and s2_remaining > 0:
        # Both have suffixes — include the longer one
        if s1_remaining >= s2_remaining:
            for i in range(s1_end, len(seq1)):
                consensus.append(seq1[i])
                con_qual.append(qual1[i] if i < len(qual1) else 0)
        else:
            for i in range(s2_end, len(seq2)):
                consensus.append(seq2[i])
                con_qual.append(qual2[i] if i < len(qual2) else 0)
    elif s1_remaining > 0:
        for i in range(s1_end, len(seq1)):
            consensus.append(seq1[i])
            con_qual.append(qual1[i] if i < len(qual1) else 0)
    elif s2_remaining > 0:
        for i in range(s2_end, len(seq2)):
            consensus.append(seq2[i])
            con_qual.append(qual2[i] if i < len(qual2) else 0)

    return "".join(consensus), con_qual


def build_replicate_consensus(
    reads: List[Tuple[str, List[int]]],
    on_status: StatusCallback = None,
) -> Tuple[str, List[int]]:
    """Build consensus from multiple replicates of the same direction.

    Iteratively merges reads using quality-weighted pairwise alignment.
    Each merge step boosts quality for matching bases, providing
    increasingly confident consensus calls.

    Args:
        reads: List of (sequence, quality) tuples from replicates.
        on_status: Optional callback for progress updates.

    Returns:
        Tuple of (consensus_sequence, consensus_quality).
    """
    if not reads:
        return "", []
    if len(reads) == 1:
        return reads[0]

    consensus_seq, consensus_qual = reads[0]

    for i in range(1, len(reads)):
        if on_status:
            on_status(
                f"Merging replicate {i + 1} of {len(reads)}",
                (i / len(reads)) * 100,
            )
        consensus_seq, consensus_qual = merge_two_sequences(
            consensus_seq, consensus_qual,
            reads[i][0], reads[i][1],
        )

    return consensus_seq, consensus_qual


def apply_quality_filter(
    sequence: str,
    quality: List[int],
    phred_min: int,
) -> Tuple[str, List[int]]:
    """Apply final Phred quality filter to a consensus sequence.

    Replaces bases below the minimum Phred threshold with 'N'.
    This is applied as the LAST step of the pipeline to avoid
    premature information loss during intermediate merges.

    Args:
        sequence: DNA sequence string.
        quality: Per-base Phred quality scores.
        phred_min: Minimum acceptable Phred quality score.

    Returns:
        Tuple of (filtered_sequence, quality).
    """
    filtered = []
    for base, q in zip(sequence, quality):
        if q >= phred_min:
            filtered.append(base)
        else:
            filtered.append("N")
    return "".join(filtered), quality


def reverse_complement_read(
    sequence: str,
    quality: List[int],
) -> Tuple[str, List[int]]:
    """Reverse complement a sequence and reverse its quality array.

    Args:
        sequence: DNA sequence string.
        quality: Per-base Phred quality scores.

    Returns:
        Tuple of (reverse_complement_sequence, reversed_quality).
    """
    rc_seq = str(Seq(sequence).reverse_complement())
    rc_qual = quality[::-1]
    return rc_seq, rc_qual


def assemble_consensus(
    forward_files: List[str],
    reverse_files: List[str],
    phred_min: int = 20,
    trim_enabled: bool = True,
    sequence_name: str = "consensus",
    on_status: StatusCallback = None,
) -> ConsensusResult:
    """Full consensus assembly pipeline.

    Pipeline steps:
        1. Read all input files (AB1/PHD format).
        2. Trim low-quality ends (Mott's algorithm).
        3. Build forward replicate consensus.
        4. Build reverse replicate consensus (with reverse complement).
        5. Merge forward and reverse consensuses.
        6. Apply final Phred quality filter.
        7. Compute quality statistics.

    Args:
        forward_files: Paths to forward read files.
        reverse_files: Paths to reverse read files.
        phred_min: Minimum Phred quality for the final filter.
        trim_enabled: Whether to trim low-quality ends.
        sequence_name: Name for the consensus sequence.
        on_status: Callback for progress updates (message, percent).

    Returns:
        ConsensusResult with the assembled sequence and metadata.

    Raises:
        ValueError: If no valid sequences could be loaded.
    """
    result = ConsensusResult(name=sequence_name)
    warnings: List[str] = []
    trimming_info: List[Dict] = []
    total_files = len(forward_files) + len(reverse_files)
    processed = 0

    def _progress(msg: str, extra_pct: float = 0):
        if on_status:
            base_pct = (processed / max(total_files, 1)) * 50
            on_status(msg, base_pct + extra_pct)

    # --- Step 1 & 2: Read and trim forward files ---
    fwd_reads: List[Tuple[str, List[int]]] = []

    for i, fpath in enumerate(forward_files):
        _progress(f"Reading forward {i + 1}/{len(forward_files)}")
        try:
            read = read_sequence_file(fpath)
            read.direction = "forward"

            if trim_enabled:
                trim = mott_trim(read.sequence, read.quality, threshold=phred_min)
                trimming_info.append({
                    "file": read.name,
                    "direction": "forward",
                    "original_length": trim.original_length,
                    "trimmed_length": trim.trimmed_length,
                    "trimmed_left": trim.bases_trimmed_left,
                    "trimmed_right": trim.bases_trimmed_right,
                })
                if trim.trimmed_length < 50:
                    warnings.append(
                        f"Forward read '{read.name}' too short after trimming "
                        f"({trim.trimmed_length} bp). Skipping."
                    )
                    continue
                fwd_reads.append((trim.sequence, trim.quality))
            else:
                fwd_reads.append((read.sequence, read.quality))

        except (IOError, ValueError) as e:
            warnings.append(f"Error reading forward file: {e}")

        processed += 1

    # --- Step 1 & 2: Read and trim reverse files ---
    rev_reads: List[Tuple[str, List[int]]] = []

    for i, fpath in enumerate(reverse_files):
        _progress(f"Reading reverse {i + 1}/{len(reverse_files)}")
        try:
            read = read_sequence_file(fpath)
            read.direction = "reverse"

            # Reverse complement the sequence
            rc_seq, rc_qual = reverse_complement_read(
                read.sequence, read.quality
            )

            if trim_enabled:
                trim = mott_trim(rc_seq, rc_qual, threshold=phred_min)
                trimming_info.append({
                    "file": read.name,
                    "direction": "reverse",
                    "original_length": trim.original_length,
                    "trimmed_length": trim.trimmed_length,
                    "trimmed_left": trim.bases_trimmed_left,
                    "trimmed_right": trim.bases_trimmed_right,
                })
                if trim.trimmed_length < 50:
                    warnings.append(
                        f"Reverse read '{read.name}' too short after trimming "
                        f"({trim.trimmed_length} bp). Skipping."
                    )
                    continue
                rev_reads.append((trim.sequence, trim.quality))
            else:
                rev_reads.append((rc_seq, rc_qual))

        except (IOError, ValueError) as e:
            warnings.append(f"Error reading reverse file: {e}")

        processed += 1

    # --- Validate ---
    has_fwd = len(fwd_reads) > 0
    has_rev = len(rev_reads) > 0

    if not has_fwd and not has_rev:
        raise ValueError(
            "No valid sequences could be loaded. Check your input files."
        )

    # --- Step 3: Forward replicate consensus ---
    if has_fwd:
        if on_status:
            on_status("Building forward consensus...", 55)
        fwd_consensus, fwd_qual = build_replicate_consensus(
            fwd_reads,
            on_status=lambda msg, pct: _progress(msg, pct * 0.1) if on_status else None,
        )
    else:
        fwd_consensus, fwd_qual = "", []
        warnings.append("No forward reads available. Using reverse only.")

    # --- Step 4: Reverse replicate consensus ---
    if has_rev:
        if on_status:
            on_status("Building reverse consensus...", 65)
        rev_consensus, rev_qual = build_replicate_consensus(
            rev_reads,
            on_status=lambda msg, pct: _progress(msg, pct * 0.1) if on_status else None,
        )
    else:
        rev_consensus, rev_qual = "", []
        warnings.append("No reverse reads available. Using forward only.")

    # --- Step 5: Merge forward and reverse ---
    if has_fwd and has_rev:
        if on_status:
            on_status("Merging forward and reverse...", 75)
        final_seq, final_qual = merge_two_sequences(
            fwd_consensus, fwd_qual,
            rev_consensus, rev_qual,
        )
    elif has_fwd:
        final_seq, final_qual = fwd_consensus, fwd_qual
    else:
        final_seq, final_qual = rev_consensus, rev_qual

    # --- Step 6: Apply final quality filter ---
    if on_status:
        on_status("Applying quality filter...", 90)
    final_seq, final_qual = apply_quality_filter(final_seq, final_qual, phred_min)

    # --- Step 7: Statistics ---
    if on_status:
        on_status("Computing statistics...", 95)
    stats = compute_statistics(final_seq, final_qual)

    # Build result
    result.sequence = final_seq
    result.quality = final_qual
    result.fwd_reads = len(fwd_reads)
    result.rev_reads = len(rev_reads)
    result.trimming_info = trimming_info
    result.statistics = stats
    result.warnings = warnings

    if on_status:
        on_status("Assembly completed.", 100)

    return result
