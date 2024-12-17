import os
from ete3 import Tree
import itertools
import csv
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import numpy as np

# Example species and methods. Update these paths as needed.
species = "Onil"
methods = {
    "ASTRAL-Pro3": "ASTRAL-pro3/astral-pro3.nwk",
    "IQ-Tree": "IQ-Tree/Onil_nc_iqtree.min4.phy.treefile",
    "TreeQMC": "TreeQMC/TreeQMC.nwk",
    "Weighted_ASTRAL": "Weighted_ASTRAL/wastral.nwk"
}

# Load trees into a dictionary
trees = {}
for method, filepath in methods.items():
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    t = Tree(filepath)
    # Ensure leaves are treated consistently
    for leaf in t:
        leaf.name = leaf.name.strip("'\"")
    trees[method] = t

methods_list = list(trees.keys())

# Function to compute RF distances handling variable returns
def compute_rf(t1, t2):
    # First, just call robinson_foulds and check how many values are returned
    result = t1.robinson_foulds(t2, unrooted_trees=True)
    # result might have 5 or 7 values depending on ete version
    if len(result) == 5:
        rf, max_rf, common_leaves, parts_t1, parts_t2 = result
    elif len(result) == 7:
        # If 7 values are returned, unpack them and ignore the extra name lists
        rf, max_rf, common_leaves, parts_t1, parts_t2, names_t1, names_t2 = result
    else:
        raise ValueError(f"Unexpected number of return values from robinson_foulds: {len(result)}")
    return rf, max_rf

# Compute pairwise RF distances
rf_results = []
for (m1, m2) in itertools.combinations(methods_list, 2):
    t1 = trees[m1]
    t2 = trees[m2]
    rf, max_rf = compute_rf(t1, t2)
    rf_results.append((m1, m2, rf, max_rf))

# Write results to CSV
comparison_csv = f"{species}_ete_rf_comparison.csv"
with open(comparison_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Tree Pair", "RF Distance", "Max RF", "Normalized RF"])
    for (m1, m2, rf, max_rf) in rf_results:
        norm_rf = rf / max_rf if max_rf != 0 else 0
        pair_name = f"{m1} vs {m2}"
        writer.writerow([pair_name, rf, max_rf, norm_rf])

print(f"Pairwise RF distances written to {comparison_csv}")

# Render each tree to an image file
for method, t in trees.items():
    out_img = f"{species}_{method}.png"
    t.render(out_img, w=500, units="px")
    print(f"Rendered {method} tree to {out_img}")

# Optional: Create a dendrogram based on RF distances
n = len(methods_list)
dist_matrix = np.zeros((n, n))
dist_matrix[:] = np.nan

for (m1, m2, rf, max_rf) in rf_results:
    i = methods_list.index(m1)
    j = methods_list.index(m2)
    dist_matrix[i, j] = rf
    dist_matrix[j, i] = rf

# Replace diagonals with 0
for i in range(n):
    dist_matrix[i, i] = 0

condensed = []
for i in range(n):
    for j in range(i+1, n):
        condensed.append(dist_matrix[i,j])
condensed = np.array(condensed)

Z = linkage(condensed, method='average')

plt.figure(figsize=(8,5))
dendrogram(Z, labels=methods_list)
plt.title(f"Hierarchical Clustering of Methods based on RF Distances ({species})")
plt.ylabel("RF Distance")
dendrogram_img = f"{species}_rf_dendrogram.png"
plt.savefig(dendrogram_img, dpi=300)
plt.show()
print(f"Dendrogram saved as {dendrogram_img}")
