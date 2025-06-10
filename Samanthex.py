import tkinter as tk
from tkinter import filedialog, OptionMenu, StringVar, messagebox
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import os
import socket
import sys

# --- Configurações Iniciais ---
cor1 = "#302f2f"
cor2 = "#787775"
arquivos_forward = []
arquivos_reverse = []
ultima_sequencia_consenso = ""
socket.setdefaulttimeout(120)


# --- Utility Functions ---
def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)


def reverter_complementar(sequencia):
    return str(Seq(sequencia).reverse_complement())


# --- Graphical Interface Functions ---
def selecionar_arquivos_forward():
    global arquivos_forward
    arquivos_forward = filedialog.askopenfilenames(
        title= "Select Forward sequences",
        filetypes=(("Files ABI e PHD", "*.ab1 *.phd.1"), ("All files", "*.*"))
    )
    arquivos_label_fwd.config(text=f"{len(arquivos_forward)} selected forward file(s)")


def selecionar_arquivos_reverse():
    global arquivos_reverse
    arquivos_reverse = filedialog.askopenfilenames(
        title="Select Reverse sequences",
        filetypes=(("Arquivos ABI e PHD", "*.ab1 *.phd.1"), ("All files", "*.*"))
    )
    arquivos_label_rev.config(text=f"{len(arquivos_reverse)} selected reverse file(s)")


# --- Funções de Leitura de Arquivos ---
def ler_arquivo_ab1(arquivo):
    record = SeqIO.read(arquivo, "abi")
    return record.seq, record.letter_annotations["phred_quality"]


def ler_arquivo_phd(arquivo):
    seq_parts, qual_parts = [], []
    in_dna_section = False
    try:
        with open(arquivo, "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith("BEGIN_DNA"):
                    in_dna_section = True; continue
                elif line.startswith("END_DNA"):
                    in_dna_section = False; continue
                if in_dna_section and line:
                    parts = line.split()
                    if len(parts) >= 2:
                        seq_parts.append(parts[0]);
                        qual_parts.append(int(parts[1]))
        if not seq_parts: return Seq(""), []
        return Seq("".join(seq_parts)), qual_parts
    except Exception as e:
        print(f"Error processing PHD file {arquivo}: {e}")
        return Seq(""), []


# --- CONSENSUS LOGIC ---

def fundir_dois_alinhamentos(seq1, qual1, seq2, qual2, phred_min):

    alignments = pairwise2.align.localms(str(seq1), str(seq2), 2, -1, -5, -2)
    if not alignments: return seq1 + seq2, qual1 + qual2

    alinh1, alinh2, score, begin, end = alignments[0]

    consenso_seq, consenso_qual = [], []
    idx1, idx2 = 0, 0

    for base1, base2 in zip(alinh1, alinh2):
        qual_val1 = qual1[idx1] if base1 != '-' and idx1 < len(qual1) else 0
        qual_val2 = qual2[idx2] if base2 != '-' and idx2 < len(qual2) else 0

        if base1 == base2:
            chosen_base = base1
            chosen_qual = max(qual_val1, qual_val2)
        elif base1 == '-':
            chosen_base = base2
            chosen_qual = qual_val2
        elif base2 == '-':
            chosen_base = base1
            chosen_qual = qual_val1
        else:  # Mismatch
            if qual_val1 >= qual_val2:
                chosen_base = base1
                chosen_qual = qual_val1
            else:
                chosen_base = base2
                chosen_qual = qual_val2

        if chosen_qual < phred_min:
            consenso_seq.append('N')
        else:
            consenso_seq.append(chosen_base)
        consenso_qual.append(chosen_qual)

        if base1 != '-': idx1 += 1
        if base2 != '-': idx2 += 1

    return Seq("".join(consenso_seq)), consenso_qual


def criar_consenso_de_replicatas(lista_sequencias, lista_qualidades, phred_min):

    if not lista_sequencias: return None, None

    consenso_seq = lista_sequencias[0]
    consenso_qual = lista_qualidades[0]

    for i in range(1, len(lista_sequencias)):
        print(f"    ...merging with replicate {i + 1}")
        consenso_seq, consenso_qual = fundir_dois_alinhamentos(
            consenso_seq, consenso_qual,
            lista_sequencias[i], lista_qualidades[i],
            phred_min
        )
    return consenso_seq, consenso_qual


# --- Main Processing Function ---
def iniciar_processamento():

    global ultima_sequencia_consenso
    if not arquivos_forward or not arquivos_reverse:
        messagebox.showerror("Erro", "Select the Forward and Reverse replicate files.")
        return

    phred_value = int(phred_minimo.get())
    nome_base = nome_sequencia.get() or "final_sequence"

    sequencias_fwd = []
    qualidades_fwd = []
    for f in arquivos_forward:
        seq, qual = (ler_arquivo_ab1 if f.endswith('.ab1') else ler_arquivo_phd)(f)
        if seq: sequencias_fwd.append(seq); qualidades_fwd.append(qual)

    sequencias_rev = []
    qualidades_rev = []
    for f in arquivos_reverse:
        seq_orig, qual_orig = (ler_arquivo_ab1 if f.endswith('.ab1') else ler_arquivo_phd)(f)
        if seq_orig:
            sequencias_rev.append(reverter_complementar(seq_orig))
            qualidades_rev.append(qual_orig[::-1])

    if not sequencias_fwd or not sequencias_rev:
        messagebox.showerror("Erro", "The sequences could not be loaded. Check the files..")
        return

    consenso_fwd, qual_fwd = criar_consenso_de_replicatas(sequencias_fwd, qualidades_fwd, phred_value)
    consenso_rev, qual_rev = criar_consenso_de_replicatas(sequencias_rev, qualidades_rev, phred_value)
    consenso_final, _ = fundir_dois_alinhamentos(consenso_fwd, qual_fwd, consenso_rev, qual_rev, phred_value)
    print("Final consensus generated!")

    ultima_sequencia_consenso = str(consenso_final)
    resultado_final = f">{nome_base}_consensus\n{ultima_sequencia_consenso}"

    # Atualizar GUI e salvar arquivo
    resultado_text.delete(1.0, tk.END)
    resultado_text.insert(tk.END, resultado_final)
    blast_btn.config(state=tk.NORMAL)

    with open(f"{nome_base}_consensus.fasta", "w") as f:
        f.write(resultado_final)

    messagebox.showinfo("Success", f"Processing complete!\nFinal sequence saved in '{nome_base}_consensus.fasta'")


# ---BLAST ---
def executar_blast():
    if not ultima_sequencia_consenso:
        messagebox.showerror("Erro", "No consensus sequence available for BLAST")
        return

    resultado_text.delete(1.0, tk.END)
    resultado_text.insert(tk.END, "Starting the BLAST search.\n")
    root.update_idletasks()

    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", ultima_sequencia_consenso)
        blast_records = NCBIXML.read(result_handle)
        hits = [f"> {aln.title}\nScore: {hsp.score}, E-value: {hsp.expect}\n" for aln in blast_records.alignments for
                hsp in aln.hsps]
        resultado_text.delete(1.0, tk.END)
        resultado_text.insert(tk.END, "\n".join(hits) if hits else "No similar sequence found.")
    except Exception as e:
        resultado_text.delete(1.0, tk.END)
        resultado_text.insert(tk.END, f"Error running BLAST: {e}")


# --- Montagem da Interface Gráfica (GUI) ---
root = tk.Tk()
root.config(bg=cor1)
root.title("Sanger Consensus Assembler - Samanthex v2.0")
root.geometry("1280x720")
root.minsize(800, 600)

main_frame = tk.Frame(root, bg=cor1)
main_frame.pack(fill="both", expand=True, padx=10, pady=10)
main_frame.grid_columnconfigure(1, weight=1)

# Widgets
font_label = ("Calibri", 14)
font_btn = ("Calibri", 12)
font_text = ("Courier New", 10)

selecionar_fwd_btn = tk.Button(main_frame, text="Select Forward sequences", command=selecionar_arquivos_forward,
                               font=font_btn, bg=cor2, fg="white", activebackground=cor1)
selecionar_fwd_btn.grid(row=0, column=0, padx=5, pady=5, sticky="ew")
arquivos_label_fwd = tk.Label(main_frame, text="No forward files selected", font=font_label, bg=cor1,
                              fg="white")
arquivos_label_fwd.grid(row=0, column=1, padx=5, pady=5, sticky="w")

selecionar_rev_btn = tk.Button(main_frame, text="Select Reverse sequences", command=selecionar_arquivos_reverse,
                               font=font_btn, bg=cor2, fg="white", activebackground=cor1)
selecionar_rev_btn.grid(row=1, column=0, padx=5, pady=5, sticky="ew")
arquivos_label_rev = tk.Label(main_frame, text="No reverse files selected", font=font_label, bg=cor1,
                              fg="white")
arquivos_label_rev.grid(row=1, column=1, padx=5, pady=5, sticky="w")

tk.Label(main_frame, text="Final name:", font=font_label, bg=cor1, fg="white").grid(row=2, column=0, padx=5, pady=5,
                                                                                    sticky="w")
nome_sequencia = tk.Entry(main_frame, font=font_label, bg=cor2, fg="white", insertbackground="white", width=40)
nome_sequencia.grid(row=2, column=1, padx=5, pady=5, sticky="w")

tk.Label(main_frame, text="Minimum Phred:", font=font_label, bg=cor1, fg="white").grid(row=3, column=0, padx=5, pady=5,
                                                                                      sticky="w")
phred_minimo = StringVar(root)
phred_minimo.set("20")
opcoes_phred = ["10", "15", "20", "25", "30", "35", "40"]
phred_menu = OptionMenu(main_frame, phred_minimo, *opcoes_phred)
phred_menu.config(bg=cor2, fg="white", font=font_btn, activebackground=cor1, highlightthickness=0)
phred_menu["menu"].config(bg=cor2, fg="white")
phred_menu.grid(row=3, column=1, padx=5, pady=5, sticky="w")

btn_frame = tk.Frame(main_frame, bg=cor1)
btn_frame.grid(row=4, column=0, columnspan=2, pady=10, sticky="w")
iniciar_btn = tk.Button(btn_frame, text="Generate Final Consensus", command=iniciar_processamento, font=font_btn,
                        bg="green", fg="white")
iniciar_btn.pack(side="left", padx=5)
blast_btn = tk.Button(btn_frame, text="BLAST", command=executar_blast, state=tk.DISABLED, font=font_btn, bg="blue",
                      fg="white")
blast_btn.pack(side="left", padx=5)

resultado_text = tk.Text(main_frame, height=15, font=font_text, bg="#1e1e1e", fg="lightgreen", insertbackground="white",
                         wrap="word")
resultado_text.grid(row=5, column=0, columnspan=2, padx=5, pady=5, sticky="nsew")
main_frame.grid_rowconfigure(5, weight=1)

root.mainloop()
