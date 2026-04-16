import numpy as np
from scipy.spatial.distance import pdist, squareform
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog, messagebox


def parse_input(text):
    lines = text.strip().split("\n")
    names = []
    data = []

    for line in lines:
        parts = line.replace(",", " ").split()
        names.append(parts[0])
        values = list(map(float, parts[1:]))
        data.append(values)

    return names, np.array(data)


def load_file():
    filepath = filedialog.askopenfilename(filetypes=[("Text files", "*.txt *.csv")])
    if not filepath:
        return

    with open(filepath, "r") as f:
        content = f.read()

    text_input.delete("1.0", tk.END)
    text_input.insert(tk.END, content)


def generate_tree():
    try:
        text = text_input.get("1.0", tk.END)
        names, data = parse_input(text)

        # liczenie macierzy odległości
        dist_array = squareform(pdist(data, metric='euclidean'))

        # konwersja do BioPython
        matrix = []
        for i in range(len(names)):
            matrix.append(list(dist_array[i, :i+1]))

        distance_matrix = DistanceMatrix(names, matrix)

        # Neighbor Joining
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(distance_matrix)

        # rysowanie
        Phylo.draw(tree)
        plt.show()

        # zapis
        Phylo.write(tree, "tree_output.nwk", "newick")

        messagebox.showinfo("Sukces", "Drzewo wygenerowane i zapisane jako tree_output.nwk")

    except Exception as e:
        messagebox.showerror("Błąd", str(e))


# ---- GUI ----
root = tk.Tk()
root.title("Generator drzewa filogenetycznego (NJ)")

frame = tk.Frame(root)
frame.pack(padx=10, pady=10)

text_input = tk.Text(frame, width=80, height=20)
text_input.pack()

btn_frame = tk.Frame(root)
btn_frame.pack(pady=5)

load_btn = tk.Button(btn_frame, text="Wczytaj plik", command=load_file)
load_btn.pack(side=tk.LEFT, padx=5)

generate_btn = tk.Button(btn_frame, text="Generuj drzewo", command=generate_tree)
generate_btn.pack(side=tk.LEFT, padx=5)

root.mainloop()
