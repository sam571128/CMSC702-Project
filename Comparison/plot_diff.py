import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.manifold import MDS
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage
import itertools

sns.set_style("whitegrid")

# List of species
species_list = ["mzeb", "onil", "salicaceae", "crataegus", "primate"]  # Adjust as needed

shape_data = []
comparison_data = []

for sp in species_list:
    shape_file = f"tree_shape_metrics_{sp}.csv"
    comparison_file = f"tree_comparison_metrics_{sp}.csv"
    
    # Check if files exist
    try:
        shape_df = pd.read_csv(shape_file)
        shape_df["Species"] = sp
        shape_data.append(shape_df)
        
        comp_df = pd.read_csv(comparison_file)
        comp_df["Species"] = sp
        comparison_data.append(comp_df)
    except FileNotFoundError:
        print(f"Missing data for {sp}. Please run analyze_tree.py for this species or remove it from species_list.")
        continue

# Concatenate all species data
shape_all = pd.concat(shape_data, ignore_index=True)
comparison_all = pd.concat(comparison_data, ignore_index=True)

### Basic Plots ###

# 1. Sackin Index Comparison Across Species
plt.figure(figsize=(8,6))
sns.barplot(data=shape_all, x="Species", y="Sackin Index", hue="Tree Label", ci=None)
plt.title("Sackin Index Across Species and Methods")
plt.ylabel("Sackin Index")
plt.legend(title="Method", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("sackin_index_comparison.png", dpi=300)
plt.show()

# 2. Robinson-Foulds Distance Comparison Across Species
plt.figure(figsize=(10,6))
sns.barplot(data=comparison_all, x="Species", y="Robinson-Foulds Distance", hue="Tree Pair", ci=None)
plt.title("Robinson-Foulds Distance Across Species and Tree Pairs")
plt.ylabel("Robinson-Foulds Distance")
plt.legend(title="Tree Pair", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("rf_distance_comparison.png", dpi=300)
plt.show()

### Additional Plots ###

# Boxplots for Tree Shape Metrics (Sackin Index)
shape_melted = shape_all.melt(id_vars=["Species", "Tree Label"], 
                              value_vars=["Sackin Index"],
                              var_name="Metric", value_name="Value")

plt.figure(figsize=(12,8))
sns.boxplot(data=shape_melted, x="Metric", y="Value", hue="Tree Label")
plt.title("Tree Shape Metrics (Sackin Index) Across Methods")
plt.ylabel("Metric Value")
plt.legend(title="Method", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("tree_shape_metrics_boxplot.png", dpi=300)
plt.show()

# Plot Cherry Count Comparison
plt.figure(figsize=(8,6))
sns.barplot(data=shape_all, x="Species", y="Cherry Count", hue="Tree Label", ci=None)
plt.title("Cherry Count Across Species and Methods")
plt.ylabel("Cherry Count")
plt.legend(title="Method", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("cherry_count_comparison.png", dpi=300)
plt.show()

# Plot Conflicting Splits Comparison
plt.figure(figsize=(10,6))
sns.barplot(data=comparison_all, x="Species", y="Conflicting Splits", hue="Tree Pair", ci=None)
plt.title("Conflicting Splits Across Species and Tree Pairs")
plt.ylabel("Conflicting Splits")
plt.legend(title="Tree Pair", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("conflicting_splits_comparison.png", dpi=300)
plt.show()

# Hierarchical Clustering Heatmap of RF Distances
comparison_all[["Method1","Method2"]] = comparison_all["Tree Pair"].str.split(" vs ", expand=True)
method_pair_mean = comparison_all.groupby(["Method1", "Method2"])["Robinson-Foulds Distance"].mean().reset_index()
all_methods = sorted(set(method_pair_mean["Method1"]).union(set(method_pair_mean["Method2"])))
dist_matrix = pd.DataFrame(index=all_methods, columns=all_methods, data=0.0)

for _, row in method_pair_mean.iterrows():
    m1, m2, dist = row["Method1"], row["Method2"], row["Robinson-Foulds Distance"]
    dist_matrix.loc[m1, m2] = dist
    dist_matrix.loc[m2, m1] = dist

dist_matrix = dist_matrix.fillna(dist_matrix.max().max()).astype(float)

plt.figure(figsize=(10,8))
sns.heatmap(dist_matrix, annot=True, fmt=".2f", cmap="viridis")
plt.title("Average Robinson-Foulds Distance Between Methods Across Species")
plt.xlabel("Method")
plt.ylabel("Method")
plt.tight_layout()
plt.savefig("rf_distance_heatmap.png", dpi=300)
plt.show()

# Check within-method comparisons
same_method_data = comparison_all[comparison_all["Method1"] == comparison_all["Method2"]]
if not same_method_data.empty:
    consistency = same_method_data.groupby("Method1")["Robinson-Foulds Distance"].mean().reset_index()
    consistency = consistency.rename(columns={"Method1":"Method", "Robinson-Foulds Distance":"Mean Within-Method RF"})
    
    plt.figure(figsize=(8,6))
    sns.barplot(data=consistency, x="Method", y="Mean Within-Method RF", palette="pastel")
    plt.title("Method Consistency (Average Within-Method RF Distance)")
    plt.ylabel("Mean Within-Method RF Distance")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("method_consistency.png", dpi=300)
    plt.show()
else:
    print("No within-method comparisons found for consistency analysis.")

# Species-Specific Method Disagreement
species_disagreement = comparison_all.groupby("Species")["Robinson-Foulds Distance"].mean().reset_index()

plt.figure(figsize=(8,6))
sns.barplot(data=species_disagreement, x="Species", y="Robinson-Foulds Distance", palette="magma")
plt.title("Average RF Distance by Species (Across All Method Pairs)")
plt.ylabel("Mean RF Distance")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("species_disagreement.png", dpi=300)
plt.show()

# Species-Method Disagreement Heatmap
species_method_dist = (
    comparison_all
    .groupby(["Species", "Method1", "Method2"])["Robinson-Foulds Distance"]
    .mean()
    .unstack(level=[1,2])
)

plt.figure(figsize=(12,8))
sns.heatmap(species_method_dist, annot=True, fmt=".2f", cmap="coolwarm")
plt.title("RF Distance Between Methods per Species")
plt.xlabel("Method Pairs")
plt.ylabel("Species")
plt.tight_layout()
plt.savefig("species_method_disagreement_heatmap.png", dpi=300)
plt.show()

# Violin Plots for Robinson-Foulds Distances
plt.figure(figsize=(10,6))
sns.violinplot(data=comparison_all, x="Species", y="Robinson-Foulds Distance", hue="Tree Pair", split=True)
plt.title("Distribution of Robinson-Foulds Distances Across Species and Tree Pairs")
plt.ylabel("Robinson-Foulds Distance")
plt.legend(title="Tree Pair", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("rf_distance_violinplot.png", dpi=300)
plt.show()

# MDS Visualization of All Trees
comparison_all["Tree1_ID"] = comparison_all["Species"] + "_" + comparison_all["Method1"]
comparison_all["Tree2_ID"] = comparison_all["Species"] + "_" + comparison_all["Method2"]

unique_trees = pd.unique(comparison_all[["Tree1_ID", "Tree2_ID"]].values.ravel())
tree_dist_matrix = pd.DataFrame(index=unique_trees, columns=unique_trees, data=0.0)

for _, row in comparison_all.iterrows():
    t1, t2, dist = row["Tree1_ID"], row["Tree2_ID"], row["Robinson-Foulds Distance"]
    tree_dist_matrix.loc[t1, t2] = dist
    tree_dist_matrix.loc[t2, t1] = dist

tree_dist_matrix = tree_dist_matrix.fillna(tree_dist_matrix.max().max()).astype(float)

mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
tree_coords = mds.fit_transform(tree_dist_matrix.values)

tree_coords_df = pd.DataFrame(tree_coords, columns=["MDS1", "MDS2"], index=tree_dist_matrix.index)
tree_coords_df["Species"] = tree_coords_df.index.str.split("_").str[0]

tree_coords_df["Method"] = (
    tree_coords_df.index
    .to_series()
    .str.split("_")
    .apply(lambda parts: "_".join(parts[1:]))
)

plt.figure(figsize=(12,8))
sns.scatterplot(data=tree_coords_df, x="MDS1", y="MDS2", hue="Method", style="Species", s=100)
plt.title("MDS Projection of All Trees Based on Robinson-Foulds Distances")
plt.xlabel("MDS Dimension 1")
plt.ylabel("MDS Dimension 2")
plt.legend(title="Method", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("all_trees_mds.png", dpi=300)
plt.show()

# Average Tree Shape Metrics (Sackin Index)
average_shape = shape_all.groupby("Tree Label")[["Sackin Index"]].mean().reset_index()
average_shape_melted = average_shape.melt(id_vars="Tree Label", var_name="Metric", value_name="Average Value")

plt.figure(figsize=(10,6))
sns.barplot(data=average_shape_melted, x="Metric", y="Average Value", hue="Tree Label", ci=None)
plt.title("Average Tree Shape Metric (Sackin Index) Across Methods")
plt.ylabel("Average Metric Value")
plt.legend(title="Method", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("average_tree_shape_metrics.png", dpi=300)
plt.show()
