#!/bin/bash

# List of input files
files=("MgBl_loh.merged.bed" "MgBn_loh.merged.bed" "MgBr_loh.merged.bed" 
       "MgC21_loh.merged.bed" "MgC25_loh.merged.bed" "MgJv2_loh.merged.bed" 
       "MgL1_loh.merged.bed" "MgL2_loh.merged.bed" "MgP_loh.merged.bed" 
       "MgVN11_loh.merged.bed" "MgVN27_loh.merged.bed" "MgVN6_loh.merged.bed")

# Loop through each file and apply the awk command
for file in "${files[@]}"; do
    # Create an output filename
    output_file="${file%.bed}_50_loh.merged.bed"
    
    # Run the awk command
    awk -F'\t' '$4 >= 50 {print}' "$file" > "$output_file"
    
    echo "Processed $file -> $output_file"
done
