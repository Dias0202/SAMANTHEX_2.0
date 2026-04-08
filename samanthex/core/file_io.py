"""File I/O for Sanger sequencing data formats (AB1, PHD).

Supports reading ABI chromatogram files (.ab1) and Phred output files (.phd.1),
extracting both nucleotide sequences and per-base Phred quality scores.
"""

import os
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq


@dataclass
class SequenceRead:
    """Represents a single sequencing read with quality data.

    Attributes:
        name: Identifier for the read (from file or filename).
        sequence: DNA sequence string (uppercase ACGTN).
        quality: Per-base Phred quality scores.
        file_path: Absolute path to the source file.
        direction: 'forward' or 'reverse' (set by the caller).
        file_format: Detected file format ('ab1' or 'phd').
    """
    name: str
    sequence: str
    quality: List[int]
    file_path: str
    direction: str = ""
    file_format: str = ""

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def mean_quality(self) -> float:
        if not self.quality:
            return 0.0
        return sum(self.quality) / len(self.quality)

    def __repr__(self) -> str:
        return (f"SequenceRead(name={self.name!r}, length={self.length}, "
                f"mean_q={self.mean_quality:.1f}, direction={self.direction!r})")


def read_ab1(file_path: str) -> SequenceRead:
    """Read an ABI (.ab1) chromatogram file.

    Uses BioPython SeqIO to parse the binary ABI format, extracting
    the base-called sequence and Phred quality scores.

    Args:
        file_path: Path to the .ab1 file.

    Returns:
        SequenceRead with sequence and quality data.

    Raises:
        IOError: If the file cannot be read or parsed.
    """
    try:
        record = SeqIO.read(file_path, "abi")
        seq_str = str(record.seq).upper()
        quality = list(record.letter_annotations.get("phred_quality", []))

        if not seq_str:
            raise ValueError("Empty sequence in ABI file")

        if len(quality) != len(seq_str):
            raise ValueError(
                f"Sequence length ({len(seq_str)}) does not match "
                f"quality length ({len(quality)})"
            )

        return SequenceRead(
            name=record.name or os.path.splitext(os.path.basename(file_path))[0],
            sequence=seq_str,
            quality=quality,
            file_path=file_path,
            file_format="ab1",
        )
    except Exception as e:
        raise IOError(
            f"Error reading ABI file '{os.path.basename(file_path)}': {e}"
        ) from e


def read_phd(file_path: str) -> SequenceRead:
    """Read a PHD (.phd.1) Phred output file.

    Parses the text-based PHD format, extracting bases and quality
    values from the BEGIN_DNA/END_DNA section.

    Args:
        file_path: Path to the .phd.1 file.

    Returns:
        SequenceRead with sequence and quality data.

    Raises:
        IOError: If the file cannot be read or parsed.
    """
    seq_parts: List[str] = []
    qual_parts: List[int] = []
    in_dna = False

    try:
        with open(file_path, "r") as f:
            for line in f:
                stripped = line.strip()
                if stripped == "BEGIN_DNA":
                    in_dna = True
                    continue
                elif stripped == "END_DNA":
                    break
                if in_dna and stripped:
                    parts = stripped.split()
                    if len(parts) >= 2:
                        seq_parts.append(parts[0].upper())
                        qual_parts.append(int(parts[1]))

        if not seq_parts:
            raise ValueError("No sequence data found in PHD file")

        name = os.path.basename(file_path)
        # Remove extensions like .phd.1
        for ext in (".phd.1", ".phd"):
            if name.lower().endswith(ext):
                name = name[: -len(ext)]
                break

        return SequenceRead(
            name=name,
            sequence="".join(seq_parts),
            quality=qual_parts,
            file_path=file_path,
            file_format="phd",
        )
    except Exception as e:
        raise IOError(
            f"Error reading PHD file '{os.path.basename(file_path)}': {e}"
        ) from e


def read_sequence_file(file_path: str) -> SequenceRead:
    """Auto-detect file format and read a sequencing file.

    Determines the format from the file extension and delegates to
    the appropriate reader function.

    Args:
        file_path: Path to the sequencing file (.ab1 or .phd.1).

    Returns:
        SequenceRead with sequence and quality data.

    Raises:
        IOError: If the file cannot be read.
        ValueError: If the file format is not recognized.
    """
    lower = file_path.lower()

    if lower.endswith(".ab1"):
        return read_ab1(file_path)
    elif lower.endswith(".phd.1") or lower.endswith(".phd"):
        return read_phd(file_path)
    else:
        # Try AB1 first (binary format), then PHD (text format)
        try:
            return read_ab1(file_path)
        except IOError:
            try:
                return read_phd(file_path)
            except IOError:
                raise ValueError(
                    f"Unrecognized file format: '{os.path.basename(file_path)}'. "
                    f"Supported formats: .ab1, .phd.1"
                )


def get_supported_filetypes() -> List[Tuple[str, str]]:
    """Return file type filters for file dialogs."""
    return [
        ("Sequencing Files", "*.ab1 *.phd.1"),
        ("ABI Files", "*.ab1"),
        ("PHD Files", "*.phd.1"),
        ("All Files", "*.*"),
    ]
