import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram

# Define the names of the isolates
isolate_names = ["MgBl", "MgBn", "MgC21", "MgC25", "MgL1", "MgP", 
                 "MgVN6", "MgVN27", "MgJv2", "MgBr", "MgL2", "MgVN11", "MgVN18"]

# Read the pairwise overlap matrix
overlap_df = pd.read_csv("pairwise_overlap_matrix.csv", index_col=0)

# Perform hierarchical clustering
Z = linkage(overlap_df, method='ward')

# Plot and save a dendrogram
plt.figure(figsize=(10, 7))
dendrogram(Z, labels=isolate_names, leaf_rotation=90)
plt.title('Hierarchical Clustering of Isolates based on LoH Regions')
plt.savefig("dendrogram.png", format="png", dpi=300)  # Save as PNG with 300 dpi
plt.show()

# Create a heatmap and save it
plt.figure(figsize=(10, 8))
sns.heatmap(overlap_df, annot=True, cmap='coolwarm', fmt='.0f', linewidths=.5, 
            xticklabels=isolate_names, yticklabels=isolate_names)

plt.title('Heatmap of Pairwise Overlaps between Isolates')
plt.savefig("heatmap.png", format="png", dpi=300)  # Save as PNG with 300 dpi
plt.show()
