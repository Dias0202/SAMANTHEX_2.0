"""NCBI BLAST integration for consensus sequence identification.

Provides functions to query NCBI's BLAST service and format results
for display. Supports configurable database, program, and result limits.
"""

import socket
from dataclasses import dataclass, field
from typing import Callable, List, Optional

from Bio.Blast import NCBIWWW, NCBIXML

# Set a generous timeout for BLAST queries
socket.setdefaulttimeout(600)

# Available BLAST databases for nucleotide queries
BLAST_DATABASES = {
    "nt": "Nucleotide collection (nt)",
    "refseq_rna": "RefSeq RNA",
    "refseq_genomic": "RefSeq Genomic",
    "est": "Expressed Sequence Tags",
}

# Available BLAST programs for nucleotide queries
BLAST_PROGRAMS = {
    "blastn": "Standard nucleotide BLAST",
    "megablast": "Megablast (highly similar sequences)",
    "dc-megablast": "Discontiguous megablast",
}


@dataclass
class BlastHit:
    """A single BLAST alignment hit.

    Attributes:
        title: Full title of the matching sequence.
        accession: GenBank accession number.
        score: Alignment bit score.
        evalue: Expect value (statistical significance).
        identities: Number of identical positions.
        align_length: Total alignment length.
        identity_pct: Percentage identity.
        query_start: Start position in query.
        query_end: End position in query.
        subject_start: Start position in subject.
        subject_end: End position in subject.
        query_seq: Query alignment string.
        match_seq: Match indicator string.
        subject_seq: Subject alignment string.
        gaps: Number of gap positions.
    """
    title: str = ""
    accession: str = ""
    score: float = 0.0
    evalue: float = 0.0
    identities: int = 0
    align_length: int = 0
    identity_pct: float = 0.0
    query_start: int = 0
    query_end: int = 0
    subject_start: int = 0
    subject_end: int = 0
    query_seq: str = ""
    match_seq: str = ""
    subject_seq: str = ""
    gaps: int = 0


def run_blast(
    sequence: str,
    program: str = "blastn",
    database: str = "nt",
    max_hits: int = 10,
    evalue_threshold: float = 0.05,
    on_status: Optional[Callable[[str], None]] = None,
) -> List[BlastHit]:
    """Run a BLAST query against NCBI's servers.

    Submits the sequence to NCBI BLAST and parses the XML results.

    Args:
        sequence: DNA sequence to search.
        program: BLAST program ('blastn', 'megablast', 'dc-megablast').
        database: Target database ('nt', 'refseq_rna', etc.).
        max_hits: Maximum number of hits to return.
        evalue_threshold: E-value cutoff for filtering results.
        on_status: Optional callback for status updates.

    Returns:
        List of BlastHit objects sorted by score (descending).

    Raises:
        ConnectionError: If NCBI cannot be reached.
        RuntimeError: If BLAST query fails.
    """
    if on_status:
        on_status("Submitting BLAST query to NCBI...")

    try:
        result_handle = NCBIWWW.qblast(
            program,
            database,
            sequence,
            descriptions=max_hits,
            alignments=max_hits,
            hitlist_size=max_hits,
            expect=evalue_threshold,
        )
    except Exception as e:
        if "urlopen" in str(e).lower() or "timeout" in str(e).lower():
            raise ConnectionError(
                "Could not connect to NCBI. Check your internet connection."
            ) from e
        raise RuntimeError(f"BLAST query failed: {e}") from e

    if on_status:
        on_status("Parsing BLAST results...")

    try:
        blast_record = NCBIXML.read(result_handle)
    except Exception as e:
        raise RuntimeError(f"Error parsing BLAST results: {e}") from e

    hits: List[BlastHit] = []

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            # Extract accession from title
            accession = ""
            title_parts = alignment.title.split("|")
            if len(title_parts) >= 4:
                accession = title_parts[3].strip()

            hit = BlastHit(
                title=alignment.title,
                accession=accession,
                score=hsp.bits,
                evalue=hsp.expect,
                identities=hsp.identities,
                align_length=hsp.align_length,
                identity_pct=(
                    hsp.identities / hsp.align_length * 100
                    if hsp.align_length > 0 else 0
                ),
                query_start=hsp.query_start,
                query_end=hsp.query_end,
                subject_start=hsp.sbjct_start,
                subject_end=hsp.sbjct_end,
                query_seq=hsp.query,
                match_seq=hsp.match,
                subject_seq=hsp.sbjct,
                gaps=hsp.gaps,
            )
            hits.append(hit)
            # Only take the best HSP per alignment
            break

    return hits


def format_blast_results(hits: List[BlastHit], max_align_display: int = 150) -> str:
    """Format BLAST hits into a human-readable text report.

    Args:
        hits: List of BlastHit objects.
        max_align_display: Maximum alignment characters to display.

    Returns:
        Formatted text string suitable for display in the GUI.
    """
    if not hits:
        return "No similar sequences found in the database.\n"

    lines = []
    lines.append(f"{'=' * 70}")
    lines.append(f"  BLAST Results — {len(hits)} hit(s) found")
    lines.append(f"{'=' * 70}\n")

    for i, hit in enumerate(hits, 1):
        lines.append(f"--- Hit {i} ---")
        lines.append(f"  Title:    {hit.title}")
        if hit.accession:
            lines.append(f"  Accession: {hit.accession}")
        lines.append(f"  Score:    {hit.score:.1f} bits")
        lines.append(f"  E-value:  {hit.evalue:.2e}")
        lines.append(
            f"  Identity: {hit.identities}/{hit.align_length} "
            f"({hit.identity_pct:.1f}%)"
        )
        if hit.gaps:
            lines.append(f"  Gaps:     {hit.gaps}")
        lines.append(
            f"  Query:    {hit.query_start} - {hit.query_end}"
        )
        lines.append(
            f"  Subject:  {hit.subject_start} - {hit.subject_end}"
        )
        lines.append("")

        # Show alignment (truncated if very long)
        q = hit.query_seq[:max_align_display]
        m = hit.match_seq[:max_align_display]
        s = hit.subject_seq[:max_align_display]
        suffix = "..." if len(hit.query_seq) > max_align_display else ""

        lines.append(f"  Query:   {q}{suffix}")
        lines.append(f"  Match:   {m}{suffix}")
        lines.append(f"  Subject: {s}{suffix}")
        lines.append("")

    return "\n".join(lines)
