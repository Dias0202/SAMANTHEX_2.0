"""Main GUI application for SAMANTHEX.

Provides a modern tkinter/ttk-based interface with:
    - File selection for forward and reverse reads.
    - Configurable assembly settings (Phred threshold, trimming).
    - Tabbed output (Consensus, Quality Report, BLAST Results).
    - Progress bar and status updates.
    - Quality profile visualization on Canvas.
    - Export and clipboard support.
"""

import os
import sys
import threading
import tkinter as tk
from datetime import datetime
from tkinter import filedialog, messagebox, ttk
from typing import List, Optional

from samanthex import __version__
from samanthex.core.blast import (
    BLAST_DATABASES,
    BLAST_PROGRAMS,
    format_blast_results,
    run_blast,
)
from samanthex.core.consensus import ConsensusResult, assemble_consensus
from samanthex.core.file_io import get_supported_filetypes
from samanthex.gui.theme import Colors, Fonts, configure_styles


class SamanthexApp(tk.Tk):
    """Main application window for SAMANTHEX."""

    def __init__(self):
        super().__init__()

        # --- State ---
        self.forward_files: List[str] = []
        self.reverse_files: List[str] = []
        self.last_result: Optional[ConsensusResult] = None
        self._processing = False

        # --- Window setup ---
        self.title(f"SAMANTHEX v{__version__} — Sanger Consensus Assembler")
        self.geometry("1100x750")
        self.minsize(900, 650)
        self.configure(bg=Colors.BG_PRIMARY)

        # App icon (ignore if not found)
        try:
            icon_path = self._resource_path("icon.ico")
            if os.path.exists(icon_path):
                self.iconbitmap(icon_path)
        except Exception:
            pass

        # --- Styles ---
        self.style = configure_styles(self)

        # --- Variables ---
        self.phred_var = tk.StringVar(value="20")
        self.trim_var = tk.BooleanVar(value=True)
        self.seq_name_var = tk.StringVar()
        self.blast_db_var = tk.StringVar(value="nt")
        self.blast_prog_var = tk.StringVar(value="blastn")
        self.blast_hits_var = tk.StringVar(value="10")

        # --- Build UI ---
        self._create_widgets()

    @staticmethod
    def _resource_path(relative_path: str) -> str:
        """Get absolute path to resource (works with PyInstaller)."""
        try:
            base = sys._MEIPASS
        except AttributeError:
            base = os.path.abspath(".")
        return os.path.join(base, relative_path)

    # ------------------------------------------------------------------ #
    #  UI Construction
    # ------------------------------------------------------------------ #

    def _create_widgets(self):
        """Build all GUI widgets."""
        # Main container with padding
        self.main = ttk.Frame(self, padding=12)
        self.main.pack(fill="both", expand=True)
        self.main.columnconfigure(0, weight=1)
        self.main.rowconfigure(4, weight=1)

        self._create_header()
        self._create_file_section()
        self._create_settings_section()
        self._create_action_bar()
        self._create_output_section()
        self._create_status_bar()

    def _create_header(self):
        """Application title bar."""
        header = ttk.Frame(self.main)
        header.grid(row=0, column=0, sticky="ew", pady=(0, 8))

        ttk.Label(
            header, text="SAMANTHEX",
            style="Title.TLabel",
        ).pack(side="left")

        ttk.Label(
            header,
            text=f"v{__version__} — Sanger Consensus Assembler",
            style="Small.TLabel",
        ).pack(side="left", padx=(8, 0), pady=(6, 0))

    def _create_file_section(self):
        """File selection panel."""
        frame = ttk.LabelFrame(self.main, text="  Input Files  ", padding=10)
        frame.grid(row=1, column=0, sticky="ew", pady=(0, 6))
        frame.columnconfigure(1, weight=1)

        # Forward files
        ttk.Button(
            frame, text="Select Forward Files",
            command=self._on_select_forward, style="TButton", width=22,
        ).grid(row=0, column=0, padx=(0, 8), pady=3, sticky="w")

        self.fwd_label = ttk.Label(
            frame, text="No files selected", style="Small.TLabel"
        )
        self.fwd_label.grid(row=0, column=1, sticky="w")

        # Reverse files
        ttk.Button(
            frame, text="Select Reverse Files",
            command=self._on_select_reverse, style="TButton", width=22,
        ).grid(row=1, column=0, padx=(0, 8), pady=3, sticky="w")

        self.rev_label = ttk.Label(
            frame, text="No files selected", style="Small.TLabel"
        )
        self.rev_label.grid(row=1, column=1, sticky="w")

    def _create_settings_section(self):
        """Assembly settings panel."""
        frame = ttk.LabelFrame(self.main, text="  Assembly Settings  ", padding=10)
        frame.grid(row=2, column=0, sticky="ew", pady=(0, 6))

        # Row 0: Sequence name + Phred
        ttk.Label(frame, text="Sequence Name:").grid(
            row=0, column=0, padx=(0, 6), pady=3, sticky="w"
        )
        name_entry = ttk.Entry(
            frame, textvariable=self.seq_name_var, width=30
        )
        name_entry.grid(row=0, column=1, padx=(0, 20), pady=3, sticky="w")

        ttk.Label(frame, text="Min. Phred:").grid(
            row=0, column=2, padx=(0, 6), pady=3, sticky="w"
        )
        phred_combo = ttk.Combobox(
            frame, textvariable=self.phred_var, width=6,
            values=["10", "15", "20", "25", "30", "35", "40"],
            state="readonly",
        )
        phred_combo.grid(row=0, column=3, padx=(0, 20), pady=3, sticky="w")

        # Trim checkbox
        ttk.Checkbutton(
            frame, text="Trim low-quality ends (Mott's algorithm)",
            variable=self.trim_var,
        ).grid(row=0, column=4, padx=(0, 6), pady=3, sticky="w")

    def _create_action_bar(self):
        """Action buttons and progress bar."""
        frame = ttk.Frame(self.main)
        frame.grid(row=3, column=0, sticky="ew", pady=(0, 6))
        frame.columnconfigure(2, weight=1)

        # Assemble button
        self.assemble_btn = ttk.Button(
            frame, text="Assemble Consensus",
            command=self._on_assemble, style="Success.TButton",
        )
        self.assemble_btn.grid(row=0, column=0, padx=(0, 6), sticky="w")

        # BLAST button
        self.blast_btn = ttk.Button(
            frame, text="Run BLAST",
            command=self._on_blast, style="Primary.TButton", state="disabled",
        )
        self.blast_btn.grid(row=0, column=1, padx=(0, 6), sticky="w")

        # Export button
        self.export_btn = ttk.Button(
            frame, text="Export Results",
            command=self._on_export, state="disabled",
        )
        self.export_btn.grid(row=0, column=3, padx=(6, 0), sticky="e")

        # Progress bar
        self.progress = ttk.Progressbar(
            frame, mode="determinate", length=300,
            style="Green.Horizontal.TProgressbar",
        )
        self.progress.grid(row=0, column=2, padx=12, sticky="ew")

    def _create_output_section(self):
        """Tabbed output area."""
        self.notebook = ttk.Notebook(self.main)
        self.notebook.grid(row=4, column=0, sticky="nsew", pady=(0, 6))

        # --- Tab 1: Consensus Sequence ---
        tab_seq = ttk.Frame(self.notebook, padding=6)
        self.notebook.add(tab_seq, text="  Consensus Sequence  ")
        tab_seq.rowconfigure(0, weight=1)
        tab_seq.columnconfigure(0, weight=1)

        self.seq_text = tk.Text(
            tab_seq, font=Fonts.CODE, wrap="word",
            bg=Colors.BG_TERTIARY, fg=Colors.SUCCESS,
            insertbackground=Colors.FG_PRIMARY,
            selectbackground=Colors.ACCENT,
            selectforeground=Colors.BG_PRIMARY,
            borderwidth=0, padx=8, pady=8,
        )
        self.seq_text.grid(row=0, column=0, sticky="nsew")

        seq_scroll = ttk.Scrollbar(tab_seq, orient="vertical",
                                   command=self.seq_text.yview)
        seq_scroll.grid(row=0, column=1, sticky="ns")
        self.seq_text.configure(yscrollcommand=seq_scroll.set)

        # Copy button below text
        seq_btns = ttk.Frame(tab_seq)
        seq_btns.grid(row=1, column=0, columnspan=2, sticky="w", pady=(4, 0))
        ttk.Button(
            seq_btns, text="Copy Sequence",
            command=self._copy_sequence,
        ).pack(side="left", padx=(0, 6))

        # --- Tab 2: Quality Report ---
        tab_quality = ttk.Frame(self.notebook, padding=6)
        self.notebook.add(tab_quality, text="  Quality Report  ")
        tab_quality.rowconfigure(1, weight=1)
        tab_quality.columnconfigure(0, weight=1)

        # Quality canvas for visualization
        self.quality_canvas = tk.Canvas(
            tab_quality, height=120,
            bg=Colors.BG_TERTIARY, highlightthickness=0,
        )
        self.quality_canvas.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        self.quality_canvas.bind("<Configure>", self._on_canvas_resize)

        # Quality text report
        self.quality_text = tk.Text(
            tab_quality, font=Fonts.CODE, wrap="word",
            bg=Colors.BG_TERTIARY, fg=Colors.FG_PRIMARY,
            insertbackground=Colors.FG_PRIMARY, borderwidth=0,
            padx=8, pady=8,
        )
        self.quality_text.grid(row=1, column=0, sticky="nsew")

        qual_scroll = ttk.Scrollbar(tab_quality, orient="vertical",
                                    command=self.quality_text.yview)
        qual_scroll.grid(row=1, column=1, sticky="ns")
        self.quality_text.configure(yscrollcommand=qual_scroll.set)

        # --- Tab 3: BLAST Results ---
        tab_blast = ttk.Frame(self.notebook, padding=6)
        self.notebook.add(tab_blast, text="  BLAST Results  ")
        tab_blast.rowconfigure(1, weight=1)
        tab_blast.columnconfigure(0, weight=1)

        # BLAST settings bar
        blast_settings = ttk.Frame(tab_blast)
        blast_settings.grid(row=0, column=0, columnspan=2, sticky="ew", pady=(0, 6))

        ttk.Label(blast_settings, text="Database:").pack(side="left", padx=(0, 4))
        blast_db = ttk.Combobox(
            blast_settings, textvariable=self.blast_db_var, width=18,
            values=list(BLAST_DATABASES.keys()), state="readonly",
        )
        blast_db.pack(side="left", padx=(0, 12))

        ttk.Label(blast_settings, text="Program:").pack(side="left", padx=(0, 4))
        blast_prog = ttk.Combobox(
            blast_settings, textvariable=self.blast_prog_var, width=14,
            values=list(BLAST_PROGRAMS.keys()), state="readonly",
        )
        blast_prog.pack(side="left", padx=(0, 12))

        ttk.Label(blast_settings, text="Max hits:").pack(side="left", padx=(0, 4))
        blast_hits = ttk.Combobox(
            blast_settings, textvariable=self.blast_hits_var, width=5,
            values=["5", "10", "20", "50"], state="readonly",
        )
        blast_hits.pack(side="left")

        # BLAST output text
        self.blast_text = tk.Text(
            tab_blast, font=Fonts.CODE, wrap="word",
            bg=Colors.BG_TERTIARY, fg=Colors.FG_PRIMARY,
            insertbackground=Colors.FG_PRIMARY, borderwidth=0,
            padx=8, pady=8,
        )
        self.blast_text.grid(row=1, column=0, sticky="nsew")

        blast_scroll = ttk.Scrollbar(tab_blast, orient="vertical",
                                     command=self.blast_text.yview)
        blast_scroll.grid(row=1, column=1, sticky="ns")
        self.blast_text.configure(yscrollcommand=blast_scroll.set)

    def _create_status_bar(self):
        """Status bar at the bottom."""
        self.status_frame = ttk.Frame(self.main)
        self.status_frame.grid(row=5, column=0, sticky="ew")

        self.status_label = ttk.Label(
            self.status_frame, text="Ready",
            style="Status.TLabel",
        )
        self.status_label.pack(fill="x")

    # ------------------------------------------------------------------ #
    #  Event Handlers
    # ------------------------------------------------------------------ #

    def _on_select_forward(self):
        """Handle forward file selection."""
        files = filedialog.askopenfilenames(
            title="Select Forward Sequencing Files",
            filetypes=get_supported_filetypes(),
        )
        if files:
            self.forward_files = list(files)
            names = [os.path.basename(f) for f in files]
            self.fwd_label.configure(
                text=f"{len(files)} file(s): {', '.join(names)}",
                style="TLabel",
            )
            self._update_status(f"{len(files)} forward file(s) selected")

    def _on_select_reverse(self):
        """Handle reverse file selection."""
        files = filedialog.askopenfilenames(
            title="Select Reverse Sequencing Files",
            filetypes=get_supported_filetypes(),
        )
        if files:
            self.reverse_files = list(files)
            names = [os.path.basename(f) for f in files]
            self.rev_label.configure(
                text=f"{len(files)} file(s): {', '.join(names)}",
                style="TLabel",
            )
            self._update_status(f"{len(files)} reverse file(s) selected")

    def _on_assemble(self):
        """Start the assembly pipeline."""
        if not self.forward_files and not self.reverse_files:
            messagebox.showwarning(
                "No Files",
                "Please select at least one forward or reverse sequencing file.",
            )
            return

        if self._processing:
            return

        self._set_processing(True)
        thread = threading.Thread(target=self._assemble_thread, daemon=True)
        thread.start()

    def _on_blast(self):
        """Start a BLAST search."""
        if not self.last_result or not self.last_result.sequence:
            messagebox.showwarning(
                "No Sequence",
                "Please assemble a consensus sequence before running BLAST.",
            )
            return

        # Check for N-only sequence
        clean = self.last_result.sequence.replace("N", "")
        if len(clean) < 20:
            messagebox.showwarning(
                "Sequence Too Short",
                "The consensus sequence has too few called bases for BLAST.",
            )
            return

        if self._processing:
            return

        self._set_processing(True)
        self.notebook.select(2)  # Switch to BLAST tab

        self.blast_text.delete("1.0", tk.END)
        self.blast_text.insert(
            tk.END,
            "Submitting BLAST query to NCBI...\n"
            "This may take several minutes depending on server load.\n\n"
            "Please wait...\n"
        )

        thread = threading.Thread(target=self._blast_thread, daemon=True)
        thread.start()

    def _on_export(self):
        """Export results to files."""
        if not self.last_result:
            return

        save_dir = filedialog.askdirectory(title="Select Export Directory")
        if not save_dir:
            return

        name = self.last_result.name
        errors = []

        # Save FASTA
        try:
            fasta_path = os.path.join(save_dir, f"{name}.fasta")
            with open(fasta_path, "w") as f:
                f.write(self.last_result.to_fasta() + "\n")
        except Exception as e:
            errors.append(f"FASTA: {e}")

        # Save quality report
        try:
            report_path = os.path.join(save_dir, f"{name}_quality_report.txt")
            with open(report_path, "w") as f:
                f.write(self._format_quality_report())
        except Exception as e:
            errors.append(f"Quality report: {e}")

        # Save BLAST results if available
        blast_content = self.blast_text.get("1.0", tk.END).strip()
        if blast_content and "BLAST Results" in blast_content:
            try:
                blast_path = os.path.join(save_dir, f"{name}_blast_results.txt")
                with open(blast_path, "w") as f:
                    f.write(blast_content)
            except Exception as e:
                errors.append(f"BLAST results: {e}")

        if errors:
            messagebox.showwarning(
                "Export Partially Failed",
                "Some files could not be saved:\n" + "\n".join(errors),
            )
        else:
            messagebox.showinfo(
                "Export Complete",
                f"Files saved to:\n{save_dir}",
            )

    def _copy_sequence(self):
        """Copy consensus sequence to clipboard."""
        if not self.last_result:
            return
        self.clipboard_clear()
        self.clipboard_append(self.last_result.sequence)
        self._update_status("Sequence copied to clipboard")

    # ------------------------------------------------------------------ #
    #  Processing Threads
    # ------------------------------------------------------------------ #

    def _assemble_thread(self):
        """Run assembly in background thread."""
        try:
            phred = int(self.phred_var.get())
            name = self.seq_name_var.get().strip()
            if not name:
                name = f"consensus_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

            def on_status(msg, pct):
                self.after(0, self._update_status, msg)
                self.after(0, self._update_progress, pct)

            result = assemble_consensus(
                forward_files=self.forward_files,
                reverse_files=self.reverse_files,
                phred_min=phred,
                trim_enabled=self.trim_var.get(),
                sequence_name=name,
                on_status=on_status,
            )

            self.last_result = result
            self.after(0, self._show_results, result)

        except Exception as e:
            msg = str(e)
            self.after(
                0, lambda: messagebox.showerror("Assembly Error", msg)
            )
        finally:
            self.after(0, self._set_processing, False)

    def _blast_thread(self):
        """Run BLAST in background thread."""
        try:
            def on_status(msg):
                self.after(0, self._update_status, msg)

            hits = run_blast(
                sequence=self.last_result.sequence,
                program=self.blast_prog_var.get(),
                database=self.blast_db_var.get(),
                max_hits=int(self.blast_hits_var.get()),
                on_status=on_status,
            )

            formatted = format_blast_results(hits)
            self.after(0, self._show_blast_results, formatted)

        except ConnectionError as e:
            msg = str(e)
            self.after(
                0, lambda: messagebox.showerror("Connection Error", msg)
            )
        except Exception as e:
            msg = str(e)
            self.after(
                0, lambda: messagebox.showerror("BLAST Error", msg)
            )
        finally:
            self.after(0, self._set_processing, False)

    # ------------------------------------------------------------------ #
    #  UI Updates
    # ------------------------------------------------------------------ #

    def _set_processing(self, state: bool):
        """Enable/disable UI during processing."""
        self._processing = state
        new_state = "disabled" if state else "normal"
        self.assemble_btn.configure(state=new_state)

        if state:
            self.blast_btn.configure(state="disabled")
            self.export_btn.configure(state="disabled")
        else:
            has_result = self.last_result is not None
            self.blast_btn.configure(
                state="normal" if has_result else "disabled"
            )
            self.export_btn.configure(
                state="normal" if has_result else "disabled"
            )

        if not state:
            self.progress["value"] = 0

    def _update_status(self, message: str):
        """Update the status bar text."""
        self.status_label.configure(text=message)
        self.update_idletasks()

    def _update_progress(self, percent: float):
        """Update the progress bar."""
        self.progress["value"] = min(percent, 100)
        self.update_idletasks()

    def _show_results(self, result: ConsensusResult):
        """Display assembly results in all tabs."""
        # --- Consensus tab ---
        self.seq_text.delete("1.0", tk.END)
        self.seq_text.insert(tk.END, result.to_fasta() + "\n")

        # Show warnings if any
        if result.warnings:
            self.seq_text.insert(tk.END, "\n--- Warnings ---\n")
            for w in result.warnings:
                self.seq_text.insert(tk.END, f"  * {w}\n")

        # --- Quality tab ---
        self._draw_quality_profile(result.quality)
        self.quality_text.delete("1.0", tk.END)
        self.quality_text.insert(tk.END, self._format_quality_report())

        # --- Switch to consensus tab ---
        self.notebook.select(0)

        # --- Auto-save FASTA ---
        try:
            fasta_name = f"{result.name}.fasta"
            with open(fasta_name, "w") as f:
                f.write(result.to_fasta() + "\n")
            self._update_status(
                f"Assembly complete! Saved as '{fasta_name}' "
                f"({result.statistics.get('length', 0)} bp)"
            )
        except Exception as e:
            self._update_status(f"Assembly complete! (Could not auto-save: {e})")

        # Low-alignment merges are not blocked, but the user must see them.
        notify = messagebox.showwarning if result.warnings else messagebox.showinfo
        warning_note = (
            f"\n\n{len(result.warnings)} warning(s) - review them in the "
            f"Consensus tab before using this sequence."
            if result.warnings else ""
        )
        notify(
            "Assembly Complete" + (" - with warnings" if result.warnings else ""),
            f"Consensus: {result.statistics.get('length', 0)} bp\n"
            f"Mean Quality: {result.statistics.get('mean_quality', 0):.1f}\n"
            f"Bases >= Q20: {result.statistics.get('q20_percent', 0):.1f}%\n"
            f"N bases: {result.statistics.get('n_count', 0)} "
            f"({result.statistics.get('n_percent', 0):.1f}%)\n"
            f"GC content: {result.statistics.get('gc_content', 0):.1f}%\n"
            f"Forward reads: {result.fwd_reads} | Reverse reads: {result.rev_reads}"
            f"{warning_note}",
        )

    def _show_blast_results(self, formatted: str):
        """Display BLAST results."""
        self.blast_text.delete("1.0", tk.END)
        self.blast_text.insert(tk.END, formatted)
        self._update_status("BLAST search completed")

    def _format_quality_report(self) -> str:
        """Build a text quality report from the last result."""
        if not self.last_result:
            return ""

        r = self.last_result
        s = r.statistics
        lines = []

        lines.append("=" * 50)
        lines.append("  QUALITY REPORT")
        lines.append("=" * 50)
        lines.append("")
        lines.append(f"  Sequence Name:   {r.name}")
        lines.append(f"  Length:           {s.get('length', 0)} bp")
        lines.append(f"  GC Content:       {s.get('gc_content', 0):.1f}%")
        lines.append("")
        lines.append("  Quality Scores:")
        lines.append(f"    Mean:           {s.get('mean_quality', 0):.1f}")
        lines.append(f"    Min:            {s.get('min_quality', 0)}")
        lines.append(f"    Max:            {s.get('max_quality', 0)}")
        lines.append(
            f"    Bases >= Q20:   {s.get('q20_count', 0)} "
            f"({s.get('q20_percent', 0):.1f}%)"
        )
        lines.append(
            f"    Bases >= Q30:   {s.get('q30_count', 0)} "
            f"({s.get('q30_percent', 0):.1f}%)"
        )
        lines.append(
            f"    N bases:        {s.get('n_count', 0)} "
            f"({s.get('n_percent', 0):.1f}%)"
        )
        lines.append("")
        lines.append(f"  Reads Used:")
        lines.append(f"    Forward:        {r.fwd_reads}")
        lines.append(f"    Reverse:        {r.rev_reads}")

        if r.trimming_info:
            lines.append("")
            lines.append("  Trimming Details:")
            for t in r.trimming_info:
                lines.append(
                    f"    {t['file']} ({t['direction']}): "
                    f"{t['original_length']} -> {t['trimmed_length']} bp "
                    f"(trimmed {t['trimmed_left']}L + {t['trimmed_right']}R)"
                )

        if r.warnings:
            lines.append("")
            lines.append("  Warnings:")
            for w in r.warnings:
                lines.append(f"    * {w}")

        lines.append("")
        return "\n".join(lines)

    # ------------------------------------------------------------------ #
    #  Quality Visualization
    # ------------------------------------------------------------------ #

    def _draw_quality_profile(self, quality: List[int]):
        """Draw a quality bar chart on the Canvas widget."""
        self.quality_canvas.delete("all")
        self.quality_canvas.update_idletasks()

        w = self.quality_canvas.winfo_width()
        h = self.quality_canvas.winfo_height()

        if not quality or w < 10 or h < 10:
            return

        n = len(quality)
        max_q = 60  # Scale: max Phred for y-axis
        margin_bottom = 18
        margin_top = 5
        plot_h = h - margin_bottom - margin_top

        bar_width = max(1, w / n)

        for i, q in enumerate(quality):
            x = i * bar_width
            bar_height = min(q / max_q, 1.0) * plot_h
            y = h - margin_bottom - bar_height

            if q >= 30:
                color = Colors.Q_HIGH
            elif q >= 20:
                color = Colors.Q_MEDIUM
            else:
                color = Colors.Q_LOW

            self.quality_canvas.create_rectangle(
                x, y, x + bar_width, h - margin_bottom,
                fill=color, outline="",
            )

        # Draw Q20 threshold line
        y20 = h - margin_bottom - (20 / max_q) * plot_h
        self.quality_canvas.create_line(
            0, y20, w, y20, fill=Colors.WARNING, dash=(4, 4), width=1,
        )
        self.quality_canvas.create_text(
            w - 4, y20 - 4, text="Q20", fill=Colors.WARNING,
            anchor="se", font=Fonts.CODE_SMALL,
        )

        # Draw Q30 threshold line
        y30 = h - margin_bottom - (30 / max_q) * plot_h
        self.quality_canvas.create_line(
            0, y30, w, y30, fill=Colors.SUCCESS, dash=(4, 4), width=1,
        )
        self.quality_canvas.create_text(
            w - 4, y30 - 4, text="Q30", fill=Colors.SUCCESS,
            anchor="se", font=Fonts.CODE_SMALL,
        )

        # X-axis label
        self.quality_canvas.create_text(
            w // 2, h - 3, text=f"Position (1 - {n})",
            fill=Colors.FG_SECONDARY, anchor="s", font=Fonts.CODE_SMALL,
        )

    def _on_canvas_resize(self, event):
        """Redraw quality profile when canvas is resized."""
        if self.last_result and self.last_result.quality:
            self._draw_quality_profile(self.last_result.quality)
