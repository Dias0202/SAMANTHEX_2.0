#!/usr/bin/env python3
"""Generate PDF documentation for SAMANTHEX.

Creates a professional PDF user manual from structured content.
Uses fpdf2 for PDF generation (no LaTeX required).

Usage:
    pip install fpdf2
    python docs/generate_pdf.py

Output:
    docs/SAMANTHEX_User_Manual.pdf
"""

import os
import sys

try:
    from fpdf import FPDF
except ImportError:
    print("fpdf2 is required to generate PDF documentation.")
    print("Install it with: pip install fpdf2")
    sys.exit(1)


class ManualPDF(FPDF):
    """Custom PDF class with header/footer for the manual."""

    def header(self):
        self.set_font("Helvetica", "I", 9)
        self.set_text_color(120, 120, 120)
        self.cell(0, 8, "SAMANTHEX v3.0 - User Manual", align="L")
        self.ln(10)

    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 8)
        self.set_text_color(120, 120, 120)
        self.cell(0, 10, f"Page {self.page_no()}/{{nb}}", align="C")

    def chapter_title(self, title):
        self.set_font("Helvetica", "B", 16)
        self.set_text_color(30, 30, 46)
        self.cell(0, 12, title, new_x="LMARGIN", new_y="NEXT")
        self.set_draw_color(137, 180, 250)
        self.line(10, self.get_y(), 200, self.get_y())
        self.ln(6)

    def section_title(self, title):
        self.set_font("Helvetica", "B", 13)
        self.set_text_color(50, 50, 80)
        self.cell(0, 10, title, new_x="LMARGIN", new_y="NEXT")
        self.ln(2)

    def body_text(self, text):
        self.set_font("Helvetica", "", 10)
        self.set_text_color(40, 40, 40)
        self.multi_cell(0, 6, text)
        self.ln(3)

    def bullet(self, text):
        self.set_font("Helvetica", "", 10)
        self.set_text_color(40, 40, 40)
        self.cell(8, 6, "-")
        self.multi_cell(0, 6, text)
        self.ln(1)

    def code_block(self, text):
        self.set_font("Courier", "", 9)
        self.set_fill_color(240, 240, 245)
        self.set_text_color(50, 50, 50)
        y_start = self.get_y()
        self.multi_cell(0, 5, text, fill=True)
        self.ln(3)

    def table_row(self, cells, bold=False, fill=False):
        style = "B" if bold else ""
        self.set_font("Helvetica", style, 9)
        if fill:
            self.set_fill_color(230, 235, 245)
        col_width = (self.w - 20) / len(cells)
        for cell in cells:
            self.cell(col_width, 7, str(cell), border=1, fill=fill)
        self.ln()


def generate_manual():
    """Generate the complete PDF user manual."""
    pdf = ManualPDF()
    pdf.alias_nb_pages()
    pdf.set_auto_page_break(auto=True, margin=20)

    # --- Cover Page ---
    pdf.add_page()
    pdf.ln(50)
    pdf.set_font("Helvetica", "B", 32)
    pdf.set_text_color(30, 30, 46)
    pdf.cell(0, 15, "SAMANTHEX", align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 14)
    pdf.set_text_color(100, 100, 120)
    pdf.cell(0, 10, "Sanger Consensus Assembler", align="C",
             new_x="LMARGIN", new_y="NEXT")
    pdf.ln(5)
    pdf.set_font("Helvetica", "B", 12)
    pdf.set_text_color(137, 180, 250)
    pdf.cell(0, 10, "Version 3.0 - User Manual", align="C",
             new_x="LMARGIN", new_y="NEXT")
    pdf.ln(20)
    pdf.set_font("Helvetica", "", 10)
    pdf.set_text_color(80, 80, 80)
    pdf.cell(0, 7, "Gabriel Dias Moreira", align="C",
             new_x="LMARGIN", new_y="NEXT")
    pdf.cell(0, 7, "Universidade Federal de Minas Gerais (UFMG)", align="C",
             new_x="LMARGIN", new_y="NEXT")
    pdf.cell(0, 7, "Graduation Program in Microbiology", align="C",
             new_x="LMARGIN", new_y="NEXT")

    # --- Chapter 1: Introduction ---
    pdf.add_page()
    pdf.chapter_title("1. Introduction")
    pdf.body_text(
        "SAMANTHEX is a bioinformatics tool designed to assemble high-quality "
        "consensus sequences from Sanger sequencing data. It provides an "
        "intuitive graphical interface for processing chromatogram files, "
        "with integrated quality trimming, multi-replicate consensus assembly, "
        "and NCBI BLAST search."
    )
    pdf.body_text(
        "The tool accepts ABI (.ab1) and PHD (.phd.1) file formats, supports "
        "multiple forward and reverse replicates per sample, and uses Phred "
        "quality scores for intelligent base calling during consensus generation."
    )

    pdf.section_title("Key Features")
    pdf.bullet("Load multiple forward and reverse read files (replicates)")
    pdf.bullet("Automatic quality trimming using Mott's algorithm (same as phred)")
    pdf.bullet("Quality-aware consensus with Phred score boosting for concordant bases")
    pdf.bullet("Configurable minimum Phred quality threshold (10-40)")
    pdf.bullet("NCBI BLAST integration with selectable database and program")
    pdf.bullet("Per-base quality visualization with Q20/Q30 reference lines")
    pdf.bullet("Detailed quality statistics (length, GC%, mean Q, %Q20, %Q30, N count)")
    pdf.bullet("Export: FASTA, quality report, BLAST results")
    pdf.bullet("Standalone executable - no Python installation required")

    # --- Chapter 2: Installation ---
    pdf.add_page()
    pdf.chapter_title("2. Installation")

    pdf.section_title("Option A: Standalone Executable (Recommended)")
    pdf.body_text(
        "Download the Samanthex.exe file from the releases page. "
        "No installation is required - simply double-click to run. "
        "This is the recommended option for users who do not have "
        "Python installed."
    )

    pdf.section_title("Option B: From Source Code")
    pdf.body_text("Requirements: Python 3.9 or later with tkinter.")
    pdf.code_block(
        "# Clone the repository\n"
        "git clone https://github.com/Dias0202/SAMANTHEX_2.0.git\n"
        "cd SAMANTHEX_2.0\n\n"
        "# Install dependencies\n"
        "pip install -r requirements.txt\n\n"
        "# Run the application\n"
        "python run.py"
    )

    pdf.section_title("Building the Executable Yourself")
    pdf.code_block(
        "pip install pyinstaller\n"
        "python build_exe.py"
    )
    pdf.body_text("The executable will be created in the dist/ directory.")

    # --- Chapter 3: Usage Guide ---
    pdf.add_page()
    pdf.chapter_title("3. Usage Guide")

    pdf.section_title("Step 1: Select Input Files")
    pdf.body_text(
        'Click "Select Forward Files" to load your forward sequencing reads, '
        'and "Select Reverse Files" for reverse reads. You can select multiple '
        "files per direction - these are treated as replicates of the same sample."
    )
    pdf.body_text(
        "Supported formats: ABI chromatogram files (.ab1) and Phred output "
        "files (.phd.1). You can mix formats within the same direction."
    )
    pdf.body_text(
        "Note: You can use only forward OR only reverse files if bidirectional "
        "sequencing is not available."
    )

    pdf.section_title("Step 2: Configure Settings")
    pdf.table_row(["Setting", "Description", "Default"], bold=True, fill=True)
    pdf.table_row(["Sequence Name", "Identifier for the output", "Auto timestamp"])
    pdf.table_row(["Min. Phred", "Quality threshold (bases below = N)", "20"])
    pdf.table_row(["Trim ends", "Enable Mott's quality trimming", "Enabled"])
    pdf.ln(3)

    pdf.section_title("Step 3: Assemble Consensus")
    pdf.body_text(
        'Click "Assemble Consensus" to start the pipeline. The progress bar '
        "shows the current stage. Results are displayed in three tabs:"
    )
    pdf.bullet(
        "Consensus Sequence - The final FASTA output with a 'Copy' button."
    )
    pdf.bullet(
        "Quality Report - Detailed statistics including length, GC content, "
        "quality distribution, trimming details, and warnings."
    )
    pdf.bullet(
        "BLAST Results - Populated after running a BLAST search."
    )

    pdf.section_title("Step 4: Run BLAST (Optional)")
    pdf.body_text(
        'After assembly, click "Run BLAST" to search the consensus against '
        "NCBI databases. You can configure the database (nt, refseq_rna, etc.), "
        "program (blastn, megablast), and number of hits. An internet "
        "connection is required. BLAST queries typically take 1-5 minutes."
    )

    pdf.section_title("Step 5: Export Results")
    pdf.body_text(
        'Click "Export Results" to save all outputs to a chosen directory: '
        "the consensus FASTA, quality report, and BLAST results (if available)."
    )

    # --- Chapter 4: Pipeline Details ---
    pdf.add_page()
    pdf.chapter_title("4. Assembly Pipeline")

    pdf.section_title("Overview")
    pdf.code_block(
        "1. READ    - Load AB1/PHD files\n"
        "2. TRIM    - Remove low-quality ends (Mott's algorithm)\n"
        "3. MERGE F - Merge forward replicates (quality-weighted)\n"
        "4. MERGE R - Reverse-complement + merge reverse replicates\n"
        "5. MERGE   - Align and merge forward + reverse consensuses\n"
        "6. FILTER  - Apply Phred threshold (low quality -> N)\n"
        "7. REPORT  - Compute statistics and display results"
    )

    pdf.section_title("Quality Trimming: Mott's Algorithm")
    pdf.body_text(
        "Sanger sequencing reads typically have low-quality bases at both ends. "
        "Mott's algorithm finds the contiguous region that maximizes the "
        "cumulative sum of (Qi - threshold) across all positions. This is "
        "mathematically equivalent to Kadane's maximum subarray algorithm."
    )
    pdf.body_text(
        "This is the same approach used by phred, one of the most widely "
        "cited programs in bioinformatics."
    )
    pdf.body_text(
        "Reference: Ewing B, Green P. Base-calling of automated sequencer "
        "traces using phred. II. Error probabilities. Genome Research. "
        "1998;8(3):186-194."
    )

    pdf.section_title("Quality-Aware Consensus")
    pdf.body_text("During sequence merging, quality scores determine base selection:")
    pdf.bullet(
        "Matching bases: Quality is boosted (Q = Q1 + Q2, capped at 93). "
        "Agreement from independent reads increases confidence."
    )
    pdf.bullet(
        "Mismatching bases: The base with the higher Phred score is selected."
    )
    pdf.bullet("Gaps: The non-gap base is used with its original quality.")
    pdf.bullet(
        "N-substitution is deferred to the final step to prevent premature "
        "information loss during intermediate merges."
    )

    pdf.section_title("Alignment Parameters")
    pdf.table_row(["Parameter", "Value", "Rationale"], bold=True, fill=True)
    pdf.table_row(["Mode", "Local", "Finds best overlap region"])
    pdf.table_row(["Match", "+2", "Standard reward for correct calls"])
    pdf.table_row(["Mismatch", "-1", "Mild penalty (errors expected)"])
    pdf.table_row(["Gap open", "-5", "Prevents spurious gaps"])
    pdf.table_row(["Gap extend", "-2", "Accommodates real indels"])

    # --- Chapter 5: Interpreting Results ---
    pdf.add_page()
    pdf.chapter_title("5. Interpreting Results")

    pdf.section_title("Quality Statistics")
    pdf.table_row(["Metric", "Good Value", "Description"], bold=True, fill=True)
    pdf.table_row(["Length", ">400 bp", "Total consensus length"])
    pdf.table_row(["Mean Quality", ">30", "Average Phred score"])
    pdf.table_row(["%Q20", ">95%", "Bases with Q>=20 (99% accuracy)"])
    pdf.table_row(["%Q30", ">80%", "Bases with Q>=30 (99.9% accuracy)"])
    pdf.table_row(["N bases", "<2%", "Bases below quality threshold"])
    pdf.table_row(["GC Content", "40-60%*", "Depends on organism"])
    pdf.ln(3)
    pdf.body_text("* Expected GC content varies by organism and gene.")

    pdf.section_title("Quality Profile Visualization")
    pdf.body_text(
        "The Quality Report tab displays a per-base quality bar chart with "
        "color coding:"
    )
    pdf.bullet("Green bars: Quality >= 30 (excellent, 99.9% base call accuracy)")
    pdf.bullet("Yellow bars: Quality 20-29 (good, 99% accuracy)")
    pdf.bullet("Red bars: Quality < 20 (poor, consider reviewing)")
    pdf.body_text(
        "Dashed lines indicate Q20 and Q30 thresholds. A high-quality "
        "consensus should have predominantly green bars."
    )

    pdf.section_title("BLAST Results")
    pdf.body_text("For each BLAST hit, the following information is shown:")
    pdf.bullet("Title: Full description of the matching database sequence")
    pdf.bullet("Score: Alignment bit score (higher = better match)")
    pdf.bullet("E-value: Statistical significance (lower = more significant)")
    pdf.bullet("Identity: Percentage of matching bases in the alignment")
    pdf.bullet("Alignment: Visual display of query vs. subject sequences")

    # --- Chapter 6: Troubleshooting ---
    pdf.add_page()
    pdf.chapter_title("6. Troubleshooting")

    pdf.section_title("Common Issues")

    pdf.set_font("Helvetica", "B", 10)
    pdf.cell(0, 8, "High N content in consensus", new_x="LMARGIN", new_y="NEXT")
    pdf.body_text(
        "If the consensus has many N bases, try: (1) lowering the Phred "
        "threshold, (2) checking input file quality in the original "
        "sequencing software, (3) ensuring forward and reverse reads are "
        "from the same sample."
    )

    pdf.set_font("Helvetica", "B", 10)
    pdf.cell(0, 8, "Short consensus after trimming", new_x="LMARGIN", new_y="NEXT")
    pdf.body_text(
        "If reads are very short after trimming, the raw data may be of "
        "poor quality. Try disabling trimming to see the full (untrimmed) "
        "consensus and evaluate the quality profile."
    )

    pdf.set_font("Helvetica", "B", 10)
    pdf.cell(0, 8, "BLAST timeout or connection error", new_x="LMARGIN", new_y="NEXT")
    pdf.body_text(
        "BLAST queries require an internet connection. If you experience "
        "timeouts, try again later - NCBI servers may be under heavy load. "
        "Consider using megablast for faster searches."
    )

    pdf.set_font("Helvetica", "B", 10)
    pdf.cell(0, 8, "Application does not start (executable)", new_x="LMARGIN", new_y="NEXT")
    pdf.body_text(
        "Some antivirus software may flag PyInstaller executables. Add an "
        "exception for Samanthex.exe in your antivirus settings. The "
        "application is safe and open-source."
    )

    # --- Save ---
    output_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(output_dir, "SAMANTHEX_User_Manual.pdf")
    pdf.output(output_path)
    print(f"PDF manual generated: {output_path}")
    return output_path


if __name__ == "__main__":
    generate_manual()
