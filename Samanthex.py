import tkinter as tk
from tkinter import filedialog, OptionMenu, StringVar, messagebox, ttk
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import os
import socket
import sys
import threading
from datetime import datetime

# --- Initial Configurations ---
color1 = "#302f2f"
color2 = "#787775"
color3 = "#3c3c3c"
forward_files = []
reverse_files = []
last_consensus_sequence = ""
socket.setdefaulttimeout(300)  

# --- Utility Functions ---
def resource_path(relative_path):
    """Get absolute path to resource, works for dev and for PyInstaller"""
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

def reverse_complement(sequence):
    """Returns the reverse complement sequence"""
    return str(Seq(sequence).reverse_complement())

# --- Graphical Interface Functions ---
def select_forward_files():
    """Allows the user to select forward files"""
    global forward_files
    files = filedialog.askopenfilenames(
        title="Select Forward Sequences",
        filetypes=(("ABI and PHD Files", "*.ab1 *.phd.1"), ("All files", "*.*"))
    )
    if files:
        forward_files = files
        forward_files_label.config(text=f"{len(forward_files)} forward file(s) selected")
        update_status("Forward files selected")

def select_reverse_files():
    """Allows the user to select reverse files"""
    global reverse_files
    files = filedialog.askopenfilenames(
        title="Select Reverse Sequences",
        filetypes=(("ABI and PHD Files", "*.ab1 *.phd.1"), ("All files", "*.*"))
    )
    if files:
        reverse_files = files
        reverse_files_label.config(text=f"{len(reverse_files)} reverse file(s) selected")
        update_status("Reverse files selected")

# --- File Reading Functions ---
def read_ab1_file(file):
    """Reads ABI (.ab1) files"""
    try:
        record = SeqIO.read(file, "abi")
        return record.seq, record.letter_annotations["phred_quality"]
    except Exception as e:
        print(f"Error reading ABI file {file}: {e}")
        return Seq(""), []

def read_phd_file(file):
    """Reads PHD (.phd.1) files"""
    seq_parts, qual_parts = [], []
    in_dna_section = False
    try:
        with open(file, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("BEGIN_DNA"):
                    in_dna_section = True
                    continue
                elif line.startswith("END_DNA"):
                    in_dna_section = False
                    continue
                if in_dna_section and line:
                    parts = line.split()
                    if len(parts) >= 2:
                        seq_parts.append(parts[0])
                        qual_parts.append(int(parts[1]))
        if not seq_parts:
            return Seq(""), []
        return Seq("".join(seq_parts)), qual_parts
    except Exception as e:
        print(f"Error processing PHD file {file}: {e}")
        return Seq(""), []

# --- CONSENSUS LOGIC ---
def merge_two_alignments(seq1, qual1, seq2, qual2, phred_min):
    """Merges two aligned sequences considering quality"""
    # Performs local alignment
    alignments = pairwise2.align.localms(str(seq1), str(seq2), 2, -1, -5, -2)
    if not alignments:
        return seq1, qual1  # Returns the first sequence if no alignment
    
    # Gets the best alignment
    align1, align2, score, begin, end = alignments[0]
    
    consensus_seq, consensus_qual = [], []
    idx1, idx2 = 0, 0

    for base1, base2 in zip(align1, align2):
        # Gets quality values for current positions
        qual_val1 = qual1[idx1] if base1 != '-' and idx1 < len(qual1) else 0
        qual_val2 = qual2[idx2] if base2 != '-' and idx2 < len(qual2) else 0

        # Decide which base to keep based on quality
        if base1 == base2:
            chosen_base = base1
            chosen_qual = max(qual_val1, qual_val2)
        elif base1 == '-':
            chosen_base = base2
            chosen_qual = qual_val2
        elif base2 == '-':
            chosen_base = base1
            chosen_qual = qual_val1
        else:  # Mismatch - choose base with higher quality
            if qual_val1 >= qual_val2:
                chosen_base = base1
                chosen_qual = qual_val1
            else:
                chosen_base = base2
                chosen_qual = qual_val2

        # Check if quality meets minimum requirement
        if chosen_qual < phred_min:
            consensus_seq.append('N')
            consensus_qual.append(0)
        else:
            consensus_seq.append(chosen_base)
            consensus_qual.append(chosen_qual)

        # Advance indices only if not gap
        if base1 != '-': 
            idx1 += 1
        if base2 != '-': 
            idx2 += 1

    return Seq("".join(consensus_seq)), consensus_qual

def create_replicate_consensus(sequence_list, quality_list, phred_min):
    """Creates a consensus from multiple replicates of the same direction"""
    if not sequence_list:
        return None, None

    consensus_seq = sequence_list[0]
    consensus_qual = quality_list[0]

    for i in range(1, len(sequence_list)):
        update_status(f"Merging replicate {i+1} of {len(sequence_list)}")
        consensus_seq, consensus_qual = merge_two_alignments(
            consensus_seq, consensus_qual,
            sequence_list[i], quality_list[i],
            phred_min
        )
    return consensus_seq, consensus_qual

# --- Main Processing Function ---
def start_processing():
    """Starts sequence processing in a separate thread"""
    if not forward_files or not reverse_files:
        messagebox.showerror("Error", "Select Forward and Reverse files.")
        return
    
    # Disable buttons during processing
    start_btn.config(state=tk.DISABLED)
    blast_btn.config(state=tk.DISABLED)
    
    # Start processing in a separate thread
    thread = threading.Thread(target=process_sequences)
    thread.daemon = True
    thread.start()

def process_sequences():
    """Processes sequences (executed in separate thread)"""
    global last_consensus_sequence
    
    try:
        phred_value = int(phred_minimum.get())
        base_name = sequence_name.get() or f"consensus_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        
        # Update status
        update_status("Reading forward files...")
        
        # Read forward sequences
        fwd_sequences = []
        fwd_qualities = []
        for i, f in enumerate(forward_files):
            update_status(f"Reading forward {i+1}/{len(forward_files)}: {os.path.basename(f)}")
            seq, qual = read_ab1_file(f) if f.endswith('.ab1') else read_phd_file(f)
            if seq:
                fwd_sequences.append(seq)
                fwd_qualities.append(qual)
        
        # Read reverse sequences (and apply reverse complement)
        update_status("Reading reverse files...")
        rev_sequences = []
        rev_qualities = []
        for i, f in enumerate(reverse_files):
            update_status(f"Reading reverse {i+1}/{len(reverse_files)}: {os.path.basename(f)}")
            orig_seq, orig_qual = read_ab1_file(f) if f.endswith('.ab1') else read_phd_file(f)
            if orig_seq:
                rev_sequences.append(reverse_complement(orig_seq))
                rev_qualities.append(orig_qual[::-1])  # Reverse quality order
        
        # Check if sequences were loaded correctly
        if not fwd_sequences or not rev_sequences:
            messagebox.showerror("Error", "Could not load sequences. Check the files.")
            return
        
        # Create consensus for each direction
        update_status("Creating forward consensus...")
        consensus_fwd, qual_fwd = create_replicate_consensus(fwd_sequences, fwd_qualities, phred_value)
        
        update_status("Creating reverse consensus...")
        consensus_rev, qual_rev = create_replicate_consensus(rev_sequences, rev_qualities, phred_value)
        
        # Merge forward and reverse consensuses
        update_status("Merging forward and reverse consensuses...")
        final_consensus, _ = merge_two_alignments(consensus_fwd, qual_fwd, consensus_rev, qual_rev, phred_value)
        
        # Update global variable
        last_consensus_sequence = str(final_consensus)
        final_result = f">{base_name}_consensus\n{last_consensus_sequence}"
        
        # Update interface in main thread
        root.after(0, update_result, final_result, base_name)
        
    except Exception as e:
        error_msg = f"Error during processing: {str(e)}"
        print(error_msg)
        root.after(0, lambda: messagebox.showerror("Error", error_msg))
    finally:
        # Re-enable buttons
        root.after(0, lambda: start_btn.config(state=tk.NORMAL))
        root.after(0, lambda: blast_btn.config(state=tk.NORMAL if last_consensus_sequence else tk.DISABLED))

def update_result(final_result, base_name):
    """Updates the interface with the processing result"""
    result_text.delete(1.0, tk.END)
    result_text.insert(tk.END, final_result)
    blast_btn.config(state=tk.NORMAL)
    
    # Save file
    try:
        with open(f"{base_name}_consensus.fasta", "w") as f:
            f.write(final_result)
        messagebox.showinfo("Success", f"Processing completed!\nSequence saved as '{base_name}_consensus.fasta'")
        update_status("Processing completed successfully")
    except Exception as e:
        messagebox.showerror("Error", f"Error saving file: {str(e)}")

def update_status(message):
    """Updates the status bar"""
    status_bar.config(text=f"Status: {message}")
    root.update_idletasks()

# --- BLAST ---
def run_blast():
    """Runs BLAST search in a separate thread"""
    if not last_consensus_sequence:
        messagebox.showerror("Error", "No consensus sequence available for BLAST")
        return
    
    # Disable button during execution
    blast_btn.config(state=tk.DISABLED)
    result_text.delete(1.0, tk.END)
    result_text.insert(tk.END, "Starting BLAST search...\nThis may take several minutes.\n")
    
    # Start BLAST in a separate thread
    thread = threading.Thread(target=run_blast_thread)
    thread.daemon = True
    thread.start()

def run_blast_thread():
    """Runs BLAST (in separate thread)"""
    try:
        root.after(0, lambda: update_status("Running BLAST..."))
        
        # Run BLAST
        result_handle = NCBIWWW.qblast(
            "blastn", 
            "nt", 
            last_consensus_sequence,
            descriptions=10,
            alignments=10,
            hitlist_size=10
        )
        
        # Read and process results
        root.after(0, lambda: update_status("Processing BLAST results..."))
        blast_records = NCBIXML.read(result_handle)
        
        # Format results
        results = []
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                results.append(f"> {alignment.title}")
                results.append(f"  Score: {hsp.score}, E-value: {hsp.expect}")
                results.append(f"  Identity: {hsp.identities}/{hsp.align_length} ({hsp.identities/hsp.align_length*100:.1f}%)")
                results.append(f"  Query: {hsp.query[0:75]}...")
                results.append(f"  Match: {hsp.match[0:75]}...")
                results.append(f"  Sbjt:  {hsp.sbjct[0:75]}...")
                results.append("")
        
        if not results:
            results = ["No similar sequence found."]
        
        # Update interface in main thread
        root.after(0, lambda: show_blast_results("\n".join(results)))
        
    except Exception as e:
        error_msg = f"Error running BLAST: {str(e)}"
        print(error_msg)
        root.after(0, lambda: messagebox.showerror("Error", error_msg))
    finally:
        root.after(0, lambda: blast_btn.config(state=tk.NORMAL))
        root.after(0, lambda: update_status("Ready"))

def show_blast_results(results):
    """Displays BLAST results in the interface"""
    result_text.delete(1.0, tk.END)
    result_text.insert(tk.END, results)
    update_status("BLAST completed")

# --- GUI Assembly ---
root = tk.Tk()
root.config(bg=color1)
root.title("Sanger Consensus Assembler - Samanthex v2.0")
root.geometry("1280x720")
root.minsize(800, 600)

# Main frame
main_frame = tk.Frame(root, bg=color1)
main_frame.pack(fill="both", expand=True, padx=10, pady=10)
main_frame.grid_columnconfigure(1, weight=1)

# Widgets
font_label = ("Calibri", 12)
font_btn = ("Calibri", 11)
font_text = ("Courier New", 10)

# File selection frame
file_frame = tk.LabelFrame(main_frame, text="File Selection", font=font_label, bg=color1, fg="white")
file_frame.grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

select_fwd_btn = tk.Button(file_frame, text="Select Forward Sequences", command=select_forward_files,
                           font=font_btn, bg=color2, fg="white", relief="flat")
select_fwd_btn.grid(row=0, column=0, padx=5, pady=5, sticky="w")
forward_files_label = tk.Label(file_frame, text="No forward file selected", font=font_label, bg=color1,
                              fg="white")
forward_files_label.grid(row=0, column=1, padx=5, pady=5, sticky="w")

select_rev_btn = tk.Button(file_frame, text="Select Reverse Sequences", command=select_reverse_files,
                           font=font_btn, bg=color2, fg="white", relief="flat")
select_rev_btn.grid(row=1, column=0, padx=5, pady=5, sticky="w")
reverse_files_label = tk.Label(file_frame, text="No reverse file selected", font=font_label, bg=color1,
                              fg="white")
reverse_files_label.grid(row=1, column=1, padx=5, pady=5, sticky="w")

# Settings frame
config_frame = tk.LabelFrame(main_frame, text="Settings", font=font_label, bg=color1, fg="white")
config_frame.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

tk.Label(config_frame, text="Sequence name:", font=font_label, bg=color1, fg="white").grid(row=0, column=0, padx=5, pady=5, sticky="w")
sequence_name = tk.Entry(config_frame, font=font_label, bg=color2, fg="white", insertbackground="white", width=40)
sequence_name.grid(row=0, column=1, padx=5, pady=5, sticky="w")

tk.Label(config_frame, text="Minimum Phred:", font=font_label, bg=color1, fg="white").grid(row=1, column=0, padx=5, pady=5, sticky="w")
phred_minimum = StringVar(root)
phred_minimum.set("20")
phred_options = ["10", "15", "20", "25", "30", "35", "40"]
phred_menu = OptionMenu(config_frame, phred_minimum, *phred_options)
phred_menu.config(bg=color2, fg="white", font=font_btn, relief="flat", highlightthickness=0)
phred_menu["menu"].config(bg=color2, fg="white")
phred_menu.grid(row=1, column=1, padx=5, pady=5, sticky="w")

# Button frame
btn_frame = tk.Frame(main_frame, bg=color1)
btn_frame.grid(row=2, column=0, columnspan=2, pady=10, sticky="w")
start_btn = tk.Button(btn_frame, text="Generate Final Consensus", command=start_processing, font=font_btn,
                      bg="green", fg="white", relief="flat")
start_btn.pack(side="left", padx=5)
blast_btn = tk.Button(btn_frame, text="Run BLAST", command=run_blast, state=tk.DISABLED, font=font_btn, 
                      bg="blue", fg="white", relief="flat")
blast_btn.pack(side="left", padx=5)

# Results area
result_frame = tk.LabelFrame(main_frame, text="Result", font=font_label, bg=color1, fg="white")
result_frame.grid(row=3, column=0, columnspan=2, padx=5, pady=5, sticky="nsew")
main_frame.grid_rowconfigure(3, weight=1)
result_frame.grid_rowconfigure(0, weight=1)
result_frame.grid_columnconfigure(0, weight=1)

result_text = tk.Text(result_frame, height=15, font=font_text, bg=color3, fg="lightgreen", 
                      insertbackground="white", wrap="word")
result_text.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

# Scrollbar for text
scrollbar = tk.Scrollbar(result_frame, orient="vertical", command=result_text.yview)
scrollbar.grid(row=0, column=1, sticky="ns")
result_text.config(yscrollcommand=scrollbar.set)

# Status bar
status_bar = tk.Label(main_frame, text="Status: Ready", relief="sunken", anchor="w", 
                      font=("Calibri", 10), bg=color1, fg="white")
status_bar.grid(row=4, column=0, columnspan=2, sticky="ew", padx=5, pady=2)

root.mainloop()
