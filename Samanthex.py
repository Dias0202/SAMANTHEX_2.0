import tkinter as tk
from tkinter import filedialog, OptionMenu, StringVar, messagebox
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import os
import socket

cor1 = "#302f2f"  # grey
cor2 = "#787775"
arquivos_forward = []
arquivos_reverse = []
ultima_sequencia_consenso = ""

socket.setdefaulttimeout(120)  # 120 seconds

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)


# Function to select forward sequence files
def selecionar_arquivos_forward():
    global arquivos_forward
    arquivos_forward = filedialog.askopenfilenames(
        title="Select the Forward sequences",
        filetypes=(("ABI and PHD Files", "*.ab1 *.phd.1"), ("All Files", "*.*"))
    )
    arquivos_label_fwd.config(text=f"{len(arquivos_forward)} forward file(s) selected")
    return arquivos_forward

# Function to select reverse sequence files
def selecionar_arquivos_reverse():
    global arquivos_reverse
    arquivos_reverse = filedialog.askopenfilenames(
        title="Select the Reverse sequences",
        filetypes=(("ABI and PHD Files", "*.ab1 *.phd.1"), ("All Files", "*.*"))
    )
    arquivos_label_rev.config(text=f"{len(arquivos_reverse)} reverse file(s) selected")
    return arquivos_reverse

# Function to read .ab1 file and extract sequence and quality
def ler_arquivo_ab1(arquivo):
    record = SeqIO.read(arquivo, "abi")
    return record.seq, record.letter_annotations["phred_quality"]

# Function to read .phd.1 file and extract sequence and quality
def ler_arquivo_phd(arquivo):
    seq = []
    qual = []
    in_sequence = False
    in_quality = False
    try:
        with open(arquivo, "r") as file:
            for line in file:
                line = line.strip()

                if line.startswith("BEGIN_SEQUENCE"):
                    in_sequence = True
                    continue
                elif line.startswith("END_SEQUENCE"):
                    in_sequence = False
                    continue

                if line.startswith("BEGIN_QUALITY"):
                    in_quality = True
                    continue
                elif line.startswith("END_QUALITY"):
                    in_quality = False
                    continue

                if in_sequence and line:
                    seq.append(line)

                if in_quality and line:
                    try:
                        qual.append(int(line))  # Convert quality to integer
                    except ValueError:
                        print(f"Invalid quality in file {arquivo}: {line}")
                        qual.append(0)  # Insert 0 for invalid values

        if not seq or not qual:
            print(f"Error reading sequence or quality in file: {arquivo}")
            return Seq(""), []

        return Seq("".join(seq)), qual

    except Exception as e:
        print(f"Error processing file {arquivo}: {e}")
        return Seq(""), []

# Function to reverse complement the sequence (for reverse reads)
def reverter_complementar(sequencia):
    return str(Seq(sequencia).reverse_complement())

# Find the overlap between two sequences
def encontrar_sobreposicao(seq1, seq2):
    """ Find the largest overlap between the end of seq1 and the start of seq2 """
    max_overlap = min(len(seq1), len(seq2))
    for i in range(max_overlap, 0, -1):
        if seq1[-i:] == seq2[:i]:
            return i
    return 0

# Assemble consensus sequence with overlap between forward and reverse
def montar_consenso_com_sobreposicao(seq_fwd, qual_fwd, seq_rev, qual_rev, phred_min):
    sobreposicao = encontrar_sobreposicao(seq_fwd, seq_rev)
    consenso = list(seq_fwd[:-sobreposicao]) if sobreposicao > 0 else list(seq_fwd)

    # Process overlapping region
    for i in range(sobreposicao):
        fwd_base = seq_fwd[-sobreposicao + i]
        rev_base = seq_rev[i]
        fwd_qual = qual_fwd[-sobreposicao + i]
        rev_qual = qual_rev[i]

        if fwd_qual >= phred_min and rev_qual >= phred_min:
            consenso.append(fwd_base if fwd_qual >= rev_qual else rev_base)
        elif fwd_qual >= phred_min:
            consenso.append(fwd_base)
        elif rev_qual >= phred_min:
            consenso.append(rev_base)
        else:
            consenso.append('N')  # Low-quality base

    # remaining part of the reverse sequence
    consenso += list(seq_rev[sobreposicao:])

    return "".join(consenso)

# Align and assemble consensus sequence
def alinhar_sequencias(sequencias_fwd, qualidades_fwd, sequencias_rev, qualidades_rev, phred_min):
    return montar_consenso_com_sobreposicao(sequencias_fwd, qualidades_fwd, sequencias_rev, qualidades_rev, phred_min)

# Display "Starting BLAST query" in the GUI
def exibir_mensagem_blast():
    resultado_text.delete(1.0, tk.END)
    resultado_text.insert(tk.END, "Starting a BLAST query...\n")
    root.update_idletasks()

# Run BLAST on the most recent consensus sequence
def executar_blast():
    global ultima_sequencia_consenso
    if not ultima_sequencia_consenso:
        messagebox.showerror("Error", "No consensus sequence available for BLAST")
        return

    exibir_mensagem_blast()
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", ultima_sequencia_consenso)
        blast_records = NCBIXML.read(result_handle)

        hits = []
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                hit_info = (
                    f"Sequence: {alignment.title}\n"
                    f"Score: {hsp.score}\n"
                    f"Identities: {hsp.identities}/{hsp.align_length}\n"
                    f"E-value: {hsp.expect}\n\n"
                )
                hits.append(hit_info)

        if hits:
            resultado_text.delete(1.0, tk.END)
            resultado_text.insert(tk.END, "\n".join(hits))
        else:
            resultado_text.delete(1.0, tk.END)
            resultado_text.insert(tk.END, "No similar sequences found")

    except Exception as e:
        print(f"Error performing BLAST query: {e}")
        resultado_text.delete(1.0, tk.END)
        resultado_text.insert(tk.END, f"Error performing BLAST query: {e}")

def iniciar_processamento():
    global ultima_sequencia_consenso
    if not arquivos_forward or not arquivos_reverse:
        messagebox.showerror("Error", "Select at least one forward and one reverse sequence")
        return

    nome_final = nome_sequencia.get()
    phred_value = int(phred_minimo.get())  # Minimum Phred score

    sequencias_fwd = []
    qualidades_fwd = []
    sequencias_rev = []
    qualidades_rev = []

    # Load forward files
    for arquivo in arquivos_forward:
        ext = os.path.splitext(arquivo)[1]
        if ext == '.ab1':
            seq, qual = ler_arquivo_ab1(arquivo)
        elif ext == '.phd.1':
            seq, qual = ler_arquivo_phd(arquivo)
        else:
            print(f"Unsupported file: {arquivo}")
            continue
        sequencias_fwd.append(seq)
        qualidades_fwd.append(qual)

    # Load reverse files
    for arquivo in arquivos_reverse:
        ext = os.path.splitext(arquivo)[1]
        if ext == '.ab1':
            seq, qual = ler_arquivo_ab1(arquivo)
        elif ext == '.phd.1':
            seq, qual = ler_arquivo_phd(arquivo)
        else:
            print(f"Unsupported file: {arquivo}")
            continue
        sequencias_rev.append(reverter_complementar(seq))
        qualidades_rev.append(qual[::-1])

    if not sequencias_fwd or not sequencias_rev:
        print("Error: At least one forward and reverse sequence required")
        return

    # Generate consensus sequence
    consenso = alinhar_sequencias(sequencias_fwd[0], qualidades_fwd[0],
                                  sequencias_rev[0], qualidades_rev[0], phred_value)

    ultima_sequencia_consenso = consenso

    blast_btn.config(state=tk.NORMAL)

    print(f"Consensus Sequence ({nome_final}): {consenso}")

    # Save the consensus sequence to a file (optional)
    with open(f"{nome_final}_consensus.fasta", "w") as f:
        f.write(f">{nome_final}_consensus\n")
        f.write(consenso)

    print("Processing completed")

# GUI setup
cor1 = "#302f2f"
cor2 = "#787775"
root = tk.Tk()
root.config(bg="#302f2f")
root.title("Consensus Sequence Assembler - Samanthex")
root.geometry("1280x720")

# Automatic resizing
root.grid_columnconfigure(0, weight=1)
root.grid_rowconfigure(0, weight=1)

# Left panel
frame_esquerda = tk.Frame(root, width=520, height=1080, bg=cor1)
frame_esquerda.pack(side="left", fill="both", expand=True)
frame_esquerda.pack_propagate(False)

# Right panel (logo and controls)
frame_direita = tk.Frame(root, width=1400, height=1080, bg=cor1)
frame_direita.pack(side="right", fill="both", expand=True)
frame_direita.pack_propagate(False)


# Configure grid in right panel
frame_direita.grid_columnconfigure(1, weight=1)

# Button to select forward sequences
selecionar_fwd_btn = tk.Button(frame_direita, text="Select Forward Sequences", command=selecionar_arquivos_forward,
                               font=("Calibri 15"), bg=cor1, fg="white")
selecionar_fwd_btn.grid(row=0, column=0, padx=10, pady=10, sticky="w")

# Label for selected forward sequences
arquivos_label_fwd = tk.Label(frame_direita, text="No forward file selected", font=("Calibri 15"), bg=cor1, fg="white")
arquivos_label_fwd.grid(row=0, column=1, padx=10, pady=10, sticky="w")

# Button to select reverse sequences
selecionar_rev_btn = tk.Button(frame_direita, text="Select Reverse Sequences", command=selecionar_arquivos_reverse,
                               font=("Calibri 15"), bg=cor1, fg="white")
selecionar_rev_btn.grid(row=1, column=0, padx=10, pady=10, sticky="w")

# Label for selected reverse sequences
arquivos_label_rev = tk.Label(frame_direita, text="No reverse file selected", font=("Calibri 15"), bg=cor1, fg="white")
arquivos_label_rev.grid(row=1, column=1, padx=10, pady=10, sticky="w")

# Input for final sequence name
tk.Label(frame_direita, text="Name of the final sequence:", font=("Calibri 15"), bg=cor1, fg="white").grid(row=2, column=0, padx=10, pady=10, sticky="w")
nome_sequencia = tk.Entry(frame_direita, width=50, font=("Calibri 15"), bg=cor1, fg="white")
nome_sequencia.grid(row=2, column=1, padx=10, pady=10, sticky="w")

# Dropdown to select minimum Phred value
tk.Label(frame_direita, text="Minimum value of Phred:", font=("Calibri 15"), bg=cor1, fg="white").grid(row=3, column=0, padx=10, pady=10, sticky="w")
phred_minimo = StringVar(root)
phred_minimo.set("15")
opcoes_phred = ["10", "15", "20", "25", "30", "35", "40"]
phred_menu = OptionMenu(frame_direita, phred_minimo, *opcoes_phred)
phred_menu.config(bg=cor1, fg="white", font=("Calibri 15"))
phred_menu.grid(row=3, column=1, padx=10, pady=10, sticky="w")

# Button to start consensus assembly
iniciar_btn = tk.Button(frame_direita, text="Start", command=iniciar_processamento, font=("Calibri 15"), bg=cor1, fg="white")
iniciar_btn.grid(row=4, column=0, padx=10, pady=10, sticky="w")

# Button to run BLAST (disabled initially)
blast_btn = tk.Button(frame_direita, text="BLAST", command=executar_blast, state=tk.DISABLED, font=("Calibri 15"), bg=cor1, fg="white")
blast_btn.grid(row=4, column=1, padx=10, pady=10, sticky="w")

# Text area for displaying BLAST results
resultado_text = tk.Text(frame_direita, height=20, width=100, font=("Calibri 10"), bg=cor1, fg="white")
resultado_text.grid(row=5, column=0, columnspan=2, padx=10, pady=10)

root.mainloop()
