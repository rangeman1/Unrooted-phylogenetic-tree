import numpy as np
from scipy.spatial.distance import pdist, squareform
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
import matplotlib.pyplot as plt

# ---- Dane (przecinki zamienione na kropki) ----

A = [0,1.742,1.509,2.669,4.835,0.780,2.161,6.044,1.744,2.689,1.329,1.894,4.729,1.136]
B = [1.742,0,2.258,1.671,3.594,1.828,3.637,6.526,0.496,2.357,1.947,1.068,6.282,1.226]
C = [1.509,2.258,0,3.114,5.605,1.037,1.966,4.642,2.387,2.695,1.262,2.512,4.413,1.330]
D = [2.669,1.671,3.114,0,3.029,2.961,4.577,7.461,1.484,3.631,2.972,2.618,7.121,2.466]
E = [4.835,3.594,5.605,3.029,0,5.259,6.847,9.694,3.341,4.993,5.070,3.972,9.383,4.612]
F = [0.780,1.828,1.037,2.961,5.259,0,1.927,4.780,1.979,2.402,1.191,1.919,4.519,0.898]
G = [2.161,3.637,1.966,4.577,6.847,1.927,0,2.934,3.690,3.340,2.060,3.513,2.673,2.611]
H = [6.044,6.526,4.642,7.461,9.694,4.780,2.934,0,6.579,5.770,4.881,6.258,0.964,5.427]
I = [1.744,0.496,2.387,1.484,3.341,1.979,3.690,6.579,0,2.430,2.016,1.207,6.325,1.321]
J = [2.689,2.357,2.695,3.631,4.993,2.402,3.340,5.770,2.430,0,2.032,2.084,5.614,1.923]
K = [1.329,1.947,1.262,2.972,5.070,1.191,2.060,4.881,2.016,2.032,0,2.067,4.604,1.294]
L = [1.894,1.068,2.512,2.618,3.972,1.919,3.513,6.258,1.207,2.084,2.067,0,6.128,1.394]
M = [4.729,6.282,4.413,7.121,9.383,4.519,2.673,0.964,6.325,5.614,4.604,6.128,0,5.208]
N = [1.136,1.226,1.330,2.466,4.612,0.898,2.611,5.427,1.321,1.923,1.294,1.394,5.208,0]

# ---- Tworzymy macierz (wiersze = taksony) ----
data = np.array([A, B, C, D, E, F, G, H, I, J, K, L, M, N])

names = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]

# ---- Liczymy macierz odległości (Euklidesowa) ----
dist_array = squareform(pdist(data, metric='euclidean'))

# ---- Konwersja do formatu BioPython ----
matrix = []
for i in range(len(names)):
    matrix.append(list(dist_array[i, :i+1]))

distance_matrix = DistanceMatrix(names, matrix)

# ---- Neighbor Joining ----
constructor = DistanceTreeConstructor()
tree = constructor.nj(distance_matrix)

# ---- Rysowanie drzewa nierootowanego ----
Phylo.draw(tree)
plt.show()

# ---- Zapis do Newick ----
Phylo.write(tree, "treeA_B_C_D_E_F_G_H_I_J_K_L_M_N.nwk", "newick")
