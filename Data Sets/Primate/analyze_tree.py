
import os
import itertools
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram

# ETE3 and DendroPy imports
from ete3 import Tree as EteTree
import dendropy
from dendropy.calculate import treecompare, treemeasure

#########################
# CONFIGURATION
#########################
species = "primate"  # Change to your species name
methods = {
    "ASTRAL-Pro3": "ASTRAL-pro3/apro3",
    "TreeQMC": "TreeQMC/treeqmc",
    "Weighted_ASTRAL": "Weighted_ASTRAL/wastral"
}

comparison_output_file = f"tree_comparison_metrics_{species}.csv"
shape_output_file = f"tree_shape_metrics_{species}.csv"
ete_rf_comparison = f"{species}_ete_rf_comparison.csv"
dendrogram_img = f"{species}_rf_dendrogram.png"


#########################
# HELPER FUNCTIONS
#########################

def read_tree_dendropy(file_path, taxon_namespace):
    tree = dendropy.Tree.get(
        path=file_path,
        schema="newick",
        preserve_underscores=True,
        taxon_namespace=taxon_namespace
    )
    # Assign default edge length if None
    for node in tree.postorder_node_iter():
        if node.edge_length is None:
            node.edge_length = 1.0
    tree.is_rooted = False
    return tree

def read_tree_ete(file_path):
    t = EteTree(file_path)
    for leaf in t:
        leaf.name = leaf.name.strip("'\"")
    return t

def compute_ete_rf(t1, t2):
    result = t1.robinson_foulds(t2, unrooted_trees=True)
    # ETE returns either 5 or 7 values:
    if len(result) == 5:
        rf, max_rf, common_leaves, parts_t1, parts_t2 = result
    elif len(result) == 7:
        rf, max_rf, common_leaves, parts_t1, parts_t2, names_t1, names_t2 = result
    else:
        raise ValueError(f"Unexpected number of return values from robinson_foulds: {len(result)}")
    return rf, max_rf, common_leaves, parts_t1, parts_t2

def is_strictly_bifurcating(tree):
    for nd in tree.internal_nodes():
        if len(nd.child_nodes()) != 2:
            return False
    return True

def compute_tree_height(tree):
    # Clone and reroot to find max distance
    tmp_tree = tree.clone()
    tmp_tree.reroot_at_midpoint(update_bipartitions=False)
    distances = [leaf.distance_from_root() for leaf in tmp_tree.leaf_node_iter()]
    return max(distances) if distances else "N/A"

def compute_total_branch_length(tree):
    return sum(e.length for e in tree.postorder_edge_iter() if e.length is not None)

def count_cherries(tree):
    # A cherry is an internal node that has exactly two leaves as children
    # We'll count how many internal nodes have exactly 2 descendants that are leaves.
    cherry_count = 0
    for nd in tree.internal_nodes():
        children = nd.child_nodes()
        if len(children) == 2 and all(child.is_leaf() for child in children):
            cherry_count += 1
    return cherry_count

def compute_tree_shape_metrics(tree, label):
    tmp_tree = tree.clone()
    tmp_tree.resolve_polytomies(update_bipartitions=True)
    sackin = treemeasure.sackin_index(tmp_tree)
    num_leaves = len(list(tmp_tree.leaf_node_iter()))
    tree_height = compute_tree_height(tmp_tree)
    total_bl = compute_total_branch_length(tmp_tree)
    cherries = count_cherries(tmp_tree)
    return {
        "Tree Label": label,
        "Sackin Index": sackin,
        "Number of Leaves": num_leaves,
        "Tree Height": tree_height,
        "Total Branch Length": total_bl,
        "Cherry Count": cherries
    }

def compare_trees_dendropy(t1, t2, label1, label2):
    rf = treecompare.symmetric_difference(t1, t2)
    # Weighted RF only if one of them is IQ-Tree (assuming IQ-Tree = Weighted)
    if "IQ-Tree" in [label1, label2]:
        try:
            wrf = treecompare.weighted_robinson_foulds_distance(t1, t2)
        except ValueError:
            wrf = "N/A"
    else:
        wrf = "N/A"
    t1_splits = set(str(bp.split_bitmask) for bp in t1.encode_bipartitions())
    t2_splits = set(str(bp.split_bitmask) for bp in t2.encode_bipartitions())
    conflicting_splits = len((t1_splits - t2_splits) | (t2_splits - t1_splits))
    return {
        "Tree Pair": f"{label1} vs {label2}",
        "Robinson-Foulds Distance": rf,
        "Weighted RF Distance": wrf,
        "Conflicting Splits": conflicting_splits
    }


#########################
# MAIN SCRIPT
#########################

# Define a single TaxonNamespace for dendropy
tns = dendropy.TaxonNamespace()

# Read trees for both DendroPy and ETE3
dendro_trees = {}
ete_trees = {}

for method, filepath in methods.items():
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    dendro_trees[method] = read_tree_dendropy(filepath, tns)
    ete_trees[method] = read_tree_ete(filepath)

methods_list = list(methods.keys())

# Compute pairwise ETE RF distances + additional metrics
rf_results_ete = []
for (m1, m2) in itertools.combinations(methods_list, 2):
    t1 = ete_trees[m1]
    t2 = ete_trees[m2]
    rf, max_rf, common_leaves, parts_t1, parts_t2 = compute_ete_rf(t1, t2)
    norm_rf = rf / max_rf if max_rf != 0 else 0
    # parts_t1, parts_t2 are sets of partitions unique to each tree
    unique_t1 = len(parts_t1)
    unique_t2 = len(parts_t2)
    rf_results_ete.append((m1, m2, rf, max_rf, norm_rf, common_leaves, unique_t1, unique_t2))

# Write ETE-based RF results to CSV
with open(ete_rf_comparison, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Tree Pair", "RF Distance", "Max RF", "Normalized RF", "Common Leaves", "Unique Partitions T1", "Unique Partitions T2"])
    for (m1, m2, rf, max_rf, norm_rf, common_leaves, unique_t1, unique_t2) in rf_results_ete:
        pair_name = f"{m1} vs {m2}"
        writer.writerow([pair_name, rf, max_rf, norm_rf, common_leaves, unique_t1, unique_t2])

print(f"ETE-based pairwise RF distances and extra metrics written to {ete_rf_comparison}")

# Render each ETE tree to an image file
for method, t in ete_trees.items():
    out_img = f"{species}_{method}.png"
    t.render(out_img, w=500, units="px")
    print(f"Rendered {method} tree to {out_img}")

# Create a dendrogram based on ETE RF distances
n = len(methods_list)
dist_matrix = np.zeros((n, n))
dist_matrix[:] = np.nan

for (m1, m2, rf, max_rf, norm_rf, common_leaves, unique_t1, unique_t2) in rf_results_ete:
    i = methods_list.index(m1)
    j = methods_list.index(m2)
    dist_matrix[i, j] = rf
    dist_matrix[j, i] = rf

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
plt.savefig(dendrogram_img, dpi=300)
plt.show()
print(f"Dendrogram saved as {dendrogram_img}")


#########################
# DENDROPY-BASED METRICS AND SHAPE ANALYSES
#########################
# Compute shape metrics and write them to a CSV
shape_results = []
for method, t in dendro_trees.items():
    shape_res = compute_tree_shape_metrics(t, method)
    shape_results.append(shape_res)

with open(shape_output_file, "w", newline="") as outf:
    headers = ["Tree Label", "Sackin Index", "Number of Leaves", "Tree Height", "Total Branch Length", "Cherry Count"]
    writer = csv.DictWriter(outf, fieldnames=headers)
    writer.writeheader()
    for res in shape_results:
        writer.writerow(res)

print("Tree shape metrics written to:", shape_output_file)

# Compute pairwise DendroPy-based comparison metrics (RF, Weighted RF [if IQ-Tree involved], Conflicts)
comparison_results = []
for (m1, m2) in itertools.combinations(methods_list, 2):
    t1 = dendro_trees[m1]
    t2 = dendro_trees[m2]
    comp_res = compare_trees_dendropy(t1, t2, m1, m2)
    comparison_results.append(comp_res)

with open(comparison_output_file, "w", newline="") as outf:
    headers = ["Tree Pair", "Robinson-Foulds Distance", "Weighted RF Distance", "Conflicting Splits"]
    writer = csv.DictWriter(outf, fieldnames=headers)
    writer.writeheader()
    for res in comparison_results:
        writer.writerow(res)

print("DendroPy-based comparison metrics written to:", comparison_output_file)

print("All analyses completed.")
