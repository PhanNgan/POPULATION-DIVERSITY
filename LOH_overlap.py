import numpy as np
import pandas as pd

# Load the bedtools multiinter output
file = "common_isolate_loh_2.bed"  # Replace with your actual file name
data = pd.read_csv(file, sep='\t', header=None)

# Extract isolate presence/absence data (columns 6 onward, excluding first 5 columns)
isolate_data = data.iloc[:, 5:]  # Starting from the 6th column to include only isolate data

# Check if the number of columns matches 13 (your number of isolates)
num_isolates = isolate_data.shape[1]
if num_isolates != 13:
    print(f"Warning: Expected 13 isolates, but found {num_isolates} columns in the data.")

# Initialize a matrix to store the overlap counts between each pair of isolates
overlap_matrix = np.zeros((num_isolates, num_isolates))

# Iterate over each row and calculate pairwise overlaps
for index, row in isolate_data.iterrows():
    overlapping_isolates = np.where(row != 0)[0]  # Get indices of overlapping isolates
    for i in range(len(overlapping_isolates)):
        for j in range(i, len(overlapping_isolates)):
            overlap_matrix[overlapping_isolates[i], overlapping_isolates[j]] += 1
            if i != j:
                overlap_matrix[overlapping_isolates[j], overlapping_isolates[i]] += 1  # Symmetric matrix

# Convert to a DataFrame and normalize by total possible overlaps
overlap_df = pd.DataFrame(overlap_matrix, columns=[f"Isolate{i+1}" for i in range(num_isolates)],
                          index=[f"Isolate{i+1}" for i in range(num_isolates)])

# Save the overlap matrix to a CSV file
overlap_df.to_csv("pairwise_overlap_matrix.csv")

print(overlap_df)
