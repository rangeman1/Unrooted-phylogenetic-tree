import tkinter as tk
from tkinter import filedialog, messagebox

from Bio import Phylo
import matplotlib.pyplot as plt
import numpy as np


# 🔢 Liczenie liści
def count_leaves(clade):
    if clade.is_terminal():
        return 1
    return sum(count_leaves(c) for c in clade.clades)


# 🧭 Kąty proporcjonalne
def assign_angles(clade, start_angle, end_angle, angles):
    if clade.is_terminal():
        angles[clade] = (start_angle + end_angle) / 2
        return

    total = sum(count_leaves(c) for c in clade.clades)
    current_angle = start_angle

    for child in clade.clades:
        proportion = count_leaves(child) / total
        span = proportion * (end_angle - start_angle)

        assign_angles(child, current_angle, current_angle + span, angles)
        current_angle += span

    angles[clade] = np.mean([angles[c] for c in clade.clades])


# 📐 Pozycje (z długościami gałęzi)
def compute_positions(clade, angles, depth=0, pos=None):
    if pos is None:
        pos = {}

    angle = angles[clade]
    pos[clade] = (depth, angle)

    for child in clade.clades:
        bl = child.branch_length if child.branch_length else 1
        compute_positions(child, angles, depth + bl, pos)

    return pos


# 🔄 Polar → XY
def to_cartesian(pos):
    xy = {}
    for clade, (r, theta) in pos.items():
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        xy[clade] = (x, y)
    return xy


# 🎨 Rysowanie drzewa ze skalą
def draw_tree(tree):
    try:
        tree.root_at_midpoint()
    except:
        pass

    angles = {}
    assign_angles(tree.root, 0, 2 * np.pi, angles)

    pos_polar = compute_positions(tree.root, angles)
    pos = to_cartesian(pos_polar)

    fig, ax = plt.subplots(figsize=(8, 8))

    # ====== GAŁĘZIE ======
    for clade in pos:
        x, y = pos[clade]
        for child in clade.clades:
            x2, y2 = pos[child]
            ax.plot([x, x2], [y, y2], color="black", linewidth=1)

    # ====== SKALA (OKRĘGI) ======
    max_radius = max(r for r, _ in pos_polar.values())

    # krok skali
    step = max_radius / 5 if max_radius > 0 else 1

    radii = np.arange(step, max_radius + step, step)

    for r in radii:
        circle = plt.Circle(
            (0, 0), r,
            color='gray',
            fill=False,
            linestyle='dotted',
            linewidth=0.8
        )
        ax.add_artist(circle)

        # podpis skali
        ax.text(
            r, 0,
            f"{r:.2f}",
            fontsize=8,
            ha='left',
            va='center',
            color='gray'
        )

    # ====== ETYKIETY ======
    label_offset = 0.4

    for leaf in tree.get_terminals():
        x, y = pos[leaf]
        r, angle = pos_polar[leaf]

        x_text = x + label_offset * np.cos(angle)
        y_text = y + label_offset * np.sin(angle)

        if leaf.name:
            ha = 'left' if np.cos(angle) > 0 else 'right'

            ax.text(
                x_text, y_text, leaf.name,
                fontsize=8,
                ha=ha,
                va='center'
            )

    # ====== WYGLĄD ======
    ax.set_aspect("equal")
    ax.axis("off")
    plt.title("Radial phylogenetic tree (ze skalą)")

    plt.tight_layout()
    plt.show()


# 📂 GUI
def load_file():
    file_path = filedialog.askopenfilename(
        title="Wybierz plik Newick (.nwk)",
        filetypes=[("Newick files", "*.nwk")]
    )

    if not file_path:
        return

    try:
        tree = Phylo.read(file_path, "newick")
        draw_tree(tree)
    except Exception as e:
        messagebox.showerror("Błąd", str(e))


# 🪟 Okno
root = tk.Tk()
root.title("Radial Tree Viewer (ze skalą)")

frame = tk.Frame(root, padx=20, pady=20)
frame.pack()

label = tk.Label(frame, text="Wczytaj plik .nwk i narysuj drzewo radialne")
label.pack(pady=10)

button = tk.Button(frame, text="Wybierz plik", command=load_file)
button.pack(pady=10)

exit_button = tk.Button(frame, text="Zamknij", command=root.quit)
exit_button.pack(pady=10)

root.mainloop()
