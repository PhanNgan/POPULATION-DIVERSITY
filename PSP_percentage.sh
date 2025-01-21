#!/bin/bash

# Create output table header

echo -e "Isolate\tNPSP_LoH\tNgene_LoH\tPSP_LoH_ratio\tNPSP_rest\tNgene_rest\tPSP_rest_ratio" > output_table_PSP_LOH.tsv

# List of isolate prefixes
isolates=("MgBl" "MgBn" "MgBr" "MgC21" "MgC25" "MgJv2" "MgL1" "MgL2" "MgP" "MgVN11" "MgVN18" "MgVN27" "MgVN6")

# Loop through each isolate
for isolate in "${isolates[@]}"; do
    # Define file name
    bed_file="${isolate}_loh.merged_50_loh.merged_concatenated.bed"
    
    # Calculate values
    NPSP_LoH=$(bedtools intersect -wa  -a PSP_gene.bed -b $bed_file | wc -l)
    Ngene_LoH=$(bedtools intersect -wa -a 15520_genes.bed -b $bed_file | wc -l)
    if [ $Ngene_LoH -ne 0 ]; then
        PSP_LoH_ratio=$(echo "scale=4; $NPSP_LoH * 100 / $Ngene_LoH" | bc)
    else
        PSP_LoH_ratio="NA"
    fi

    NPSP_rest=$(bedtools intersect -v -a PSP_gene.bed  -b $bed_file | wc -l)
    Ngene_rest=$(bedtools intersect -v -a 15520_genes.bed -b $bed_file | wc -l)
    if [ $Ngene_rest -ne 0 ]; then
        PSP_rest_ratio=$(echo "scale=4; $NPSP_rest * 100 / $Ngene_rest" | bc)
    else
        PSP_rest_ratio="NA"
    fi

    # Append results to the output table
    echo -e "$isolate\t$NPSP_LoH\t$Ngene_LoH\t$PSP_LoH_ratio\t$NPSP_rest\t$Ngene_rest\t$PSP_rest_ratio" >>  output_table_PSP_LOH.tsv
done

echo "Calculation complete. Results saved in output_table_PSP_LOH.tsv "
