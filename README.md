# SAMANTHEX v3.0 — Sanger Consensus Assembler

**SAMANTHEX** is a graphical bioinformatics tool for assembling high-quality consensus sequences from Sanger sequencing data. It provides an intuitive interface for loading, trimming, aligning, and generating consensus sequences from forward and reverse reads, with integrated NCBI BLAST search.

Developed at Universidade Federal de Minas Gerais (UFMG).

---

## Features

| Feature | Description |
|---|---|
| **Multi-replicate support** | Load multiple forward and/or reverse read files (replicates) per sample |
| **Quality trimming** | Mott's modified algorithm (same as phred) for automatic removal of low-quality ends |
| **Quality-aware consensus** | Bases are selected using Phred quality scores; matching bases get boosted confidence |
| **Configurable Phred threshold** | User chooses the minimum quality score (10–40) for the final consensus |
| **NCBI BLAST integration** | Search consensus against NCBI databases directly from the GUI |
| **Quality visualization** | Color-coded per-base quality bar chart with Q20/Q30 reference lines |
| **Quality statistics** | Length, GC content, mean Phred, %Q20, %Q30, N count |
| **Multiple export formats** | FASTA, quality report (TXT), BLAST results (TXT) |
| **Standalone executable** | Single-file `.exe` for Windows — no Python installation required |
| **Modern dark GUI** | Professional ttk-based interface optimized for lab environments |

---

## Supported File Formats

- **ABI (.ab1)** — Binary chromatogram files from Applied Biosystems sequencers
- **PHD (.phd.1)** — Text-based Phred output files

---

## Assembly Pipeline

SAMANTHEX follows established best practices for Sanger consensus generation:

```
1. READ         Load AB1/PHD files for forward and reverse directions
                ↓
2. TRIM         Remove low-quality ends using Mott's algorithm
                (Kadane's maximum subarray to find best quality region)
                ↓
3. FWD MERGE    Merge forward replicates with quality-weighted consensus
                ↓
4. REV MERGE    Reverse-complement reverse reads, then merge replicates
                ↓
5. F/R MERGE    Align and merge forward + reverse consensuses
                (Local pairwise alignment, quality-based conflict resolution)
                ↓
6. FILTER       Apply Phred quality threshold (low-quality bases → N)
                ↓
7. REPORT       Compute statistics and display results
```

### Quality-Aware Base Calling

During sequence merging, SAMANTHEX uses Phred quality scores to resolve conflicts:

- **Matching bases** → Quality is boosted: `Q = min(Q1 + Q2, 93)`. Agreement from independent reads increases confidence.
- **Mismatching bases** → The base with the higher quality score is selected.
- **Gaps** → The non-gap base is used with its original quality.
- **Final filter** → Only applied at the last step. Bases below the Phred threshold are replaced with `N`.

### Mott's Trimming Algorithm

Low-quality ends are trimmed using Mott's modified algorithm, which finds the contiguous region maximizing `Σ(Qi - threshold)`. This is the same approach used by phred/phrap and is equivalent to Kadane's maximum subarray algorithm.

**Reference:** Ewing B, Green P. *Base-calling of automated sequencer traces using phred. II. Error probabilities.* Genome Research. 1998;8(3):186-194.

---

## Installation

### Option A: Standalone Executable (Recommended for End Users)

Download `Samanthex.exe` from the [Releases page](https://github.com/Dias0202/SAMANTHEX_2.0/releases/latest) — **no installation required**. Just double-click to run. The executable is not stored in the Git repository; it is built with `build_exe.py` and attached to each release.

### Option B: From Source

**Requirements:** Python 3.9+ with tkinter

```bash
# 1. Clone the repository
git clone https://github.com/Dias0202/SAMANTHEX_2.0.git
cd SAMANTHEX_2.0

# 2. Install dependencies
pip install -r requirements.txt

# 3. Run the application
python run.py
```

### Building the Executable

To build the standalone `.exe` yourself:

```bash
pip install pyinstaller
python build_exe.py
```

The executable will be created in `dist/Samanthex.exe`.

---

## Usage Guide

### 1. Select Input Files

Click **"Select Forward Files"** and **"Select Reverse Files"** to load your sequencing data. You can select multiple files per direction (replicates of the same sample).

> **Tip:** You can use only forward OR only reverse files if bidirectional sequencing is not available.

### 2. Configure Settings

| Setting | Description | Default |
|---|---|---|
| **Sequence Name** | Identifier for the output FASTA header | Auto-generated timestamp |
| **Min. Phred** | Minimum quality score (bases below this become N) | 20 |
| **Trim low-quality ends** | Enable/disable Mott's algorithm trimming | Enabled |

### 3. Assemble Consensus

Click **"Assemble Consensus"**. The progress bar shows pipeline advancement. Results appear in three tabs:

- **Consensus Sequence** — The final FASTA output. Click "Copy Sequence" to copy to clipboard.
- **Quality Report** — Detailed statistics including length, GC%, quality distribution, trimming details, and per-read information.
- **BLAST Results** — (After running BLAST) Shows top database matches with scores, E-values, and alignments.

### 4. Run BLAST (Optional)

After assembly, click **"Run BLAST"** to search the consensus against NCBI. Configure:

- **Database**: `nt` (nucleotide collection), `refseq_rna`, `refseq_genomic`, `est`
- **Program**: `blastn`, `megablast` (for highly similar sequences), `dc-megablast`
- **Max hits**: Number of results to display (5–50)

> **Note:** BLAST queries require an internet connection and may take several minutes.

### 5. Export Results

Click **"Export Results"** to save all outputs to a chosen directory:

- `{name}.fasta` — Consensus sequence
- `{name}_quality_report.txt` — Detailed quality statistics
- `{name}_blast_results.txt` — BLAST results (if available)

---

## Project Structure

```
SAMANTHEX_2.0/
├── samanthex/                  # Main package
│   ├── __init__.py             # Version and metadata
│   ├── main.py                 # Application entry point
│   ├── core/                   # Bioinformatics modules
│   │   ├── file_io.py          # AB1/PHD file reading
│   │   ├── trimming.py         # Mott's quality trimming algorithm
│   │   ├── consensus.py        # Alignment and consensus assembly
│   │   └── blast.py            # NCBI BLAST integration
│   └── gui/                    # User interface modules
│       ├── theme.py            # Dark theme, colors, fonts, ttk styles
│       └── app.py              # Main GUI application class
├── run.py                      # Launch script
├── build_exe.py                # Standalone executable builder
├── docs/
│   └── generate_pdf.py         # PDF documentation generator
├── requirements.txt            # Python dependencies
├── README.md                   # This file
└── LICENSE                     # MIT License
```

### Module Descriptions

| Module | Responsibility |
|---|---|
| `core/file_io.py` | Reads AB1 and PHD sequencing files. `SequenceRead` dataclass encapsulates sequence + quality + metadata. |
| `core/trimming.py` | Mott's algorithm implementation with sliding-window fallback. Returns `TrimResult` with trim positions and statistics. |
| `core/consensus.py` | Full pipeline: read → trim → replicate merge → F/R merge → filter. Uses BioPython `PairwiseAligner` for alignment. Quality boosting for concordant bases. |
| `core/blast.py` | NCBI BLAST query submission and XML result parsing. Configurable database, program, and hit count. |
| `gui/theme.py` | Centralized dark color palette (Catppuccin-inspired), font definitions, and ttk style configuration. |
| `gui/app.py` | `SamanthexApp(tk.Tk)` — full GUI with file selection, settings, tabbed output (Sequence/Quality/BLAST), Canvas quality plot, progress bar, export, and clipboard support. |

---

## Dependencies

| Package | Version | Purpose |
|---|---|---|
| **biopython** | ≥ 1.80 | Sequence I/O, pairwise alignment (`PairwiseAligner`), BLAST interface |
| **numpy** | ≥ 1.21 | Numerical operations (used internally by BioPython) |
| **tkinter** | (built-in) | GUI framework |

---

## Technical Notes

### Why PairwiseAligner instead of pairwise2?

BioPython deprecated `pairwise2` in version 1.80. SAMANTHEX v3.0 uses the modern `Bio.Align.PairwiseAligner` class, which is faster, more memory-efficient, and returns alignment coordinates that enable precise quality tracking.

### Why defer N-substitution to the final step?

In the legacy version, bases below the Phred threshold were replaced with `N` during intermediate merges. This caused problems because:

1. `N` characters degrade alignment quality in subsequent merge steps.
2. A low-quality base in one replicate could mask a high-quality base from another replicate.

In v3.0, quality scores are tracked throughout the pipeline, and the Phred filter is only applied at the very end.

### Alignment Parameters

| Parameter | Value | Rationale |
|---|---|---|
| Mode | Local | Finds the best overlapping region without forcing end-to-end alignment |
| Match | +2 | Standard reward for correct base calls |
| Mismatch | -1 | Mild penalty; sequencing errors are expected |
| Gap open | -5 | Moderate penalty to prevent spurious gaps |
| Gap extend | -2 | Lower than open penalty to accommodate real insertions/deletions |

---

## License

This project is licensed under the **MIT License**. See [LICENSE](LICENSE) for details.

Copyright (c) 2025 Gabriel Dias Moreira

## Financial Support

This study was funded by the National Institute of Allergy and Infectious Diseases (NIAID) through the Centers for Research in Emerging Infectious Diseases (CREID) "The Coordinating Research on Emerging Arboviral Threats Encompassing the Neotropics (CREATE-NEO)" grant U01 AI151807 to Nikos Vasilakis and Kathryn A. Hanley. Betânia Paiva Drumond's research is also supported by Conselho Nacional de Desenvolvimento Científico e Tecnológico do Ministério da Ciência e Tecnologia e Inovação (CNPq). GDM is a student from the Graduation Program in Microbiology of Universidade Federal de Minas Gerais (UFMG), which is supported by the Fundação de Amparo à Pesquisa do Estado de Minas Gerais (FAPEMIG).
