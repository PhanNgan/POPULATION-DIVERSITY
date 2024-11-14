import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram

# Read the pairwise overlap matrix
overlap_df = pd.read_csv("pairwise_overlap_matrix.csv", index_col=0)

# Perform hierarchical clustering
Z = linkage(overlap_df, method='ward')

# Plot a dendrogram
plt.figure(figsize=(10, 7))
dendrogram(Z, labels=overlap_df.index, leaf_rotation=90)
plt.title('Hierarchical Clustering of Isolates based on LoH Regions')
plt.savefig("dendrogram.png", format="png", dpi=300)  # Save as PNG with 300 dpi
plt.show()

# Create a heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(overlap_df, annot=True, cmap='coolwarm', linewidths=.5)
plt.title('Heatmap of Pairwise Overlaps between Isolates')
plt.savefig("heatmap.png", format="png", dpi=300)
plt.show()

