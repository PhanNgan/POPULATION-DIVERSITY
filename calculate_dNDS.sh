#!/bin/bash

# Output file for the results
output_file="dn_ds_ratios.tsv"
echo -e "Isolate\tdN_Mg_PSP_LoH/dS_Mg_PSP_LoH\tdN_Mg_PSP_rest/dS_Mg_PSP_rest" > $output_file

# List of isolate prefixes
isolates=(MgBl MgBn MgBr MgC21 MgC25 MgJv2 MgL1 MgL2 MgP MgVN11 MgVN18 MgVN27 MgVN6)

# Loop through each isolate
for isolate in "${isolates[@]}"; do
    # Define input files
    loh_file="${isolate}_loh.merged_50_loh.merged_concatenated.bed"
    psp_loh_file="${isolate}_PSP_LoH.bed"
    psp_rest_file="${isolate}_PSP-rest.bed"

    # Intersect to create PSP_LoH regions
    bedtools intersect -wb -a "$loh_file" -b PSP_gene.bed > "$psp_loh_file"

    # Filter VCF for PSP_LoH
    bedtools intersect -wa -a synonymous_missense_variants.vcf -b "$psp_loh_file" > "synonymous_missense_variants_PSP_${isolate}_LoH.vcf"

    # Calculate dN and dS for PSP_LoH
    dS_PSP_LoH=$(grep -c "synonymous" "synonymous_missense_variants_PSP_${isolate}_LoH.vcf")
    dN_PSP_LoH=$(grep -c "missense" "synonymous_missense_variants_PSP_${isolate}_LoH.vcf")

    if [ $dS_PSP_LoH -ne 0 ]; then
        ratio_PSP_LoH=$(echo "scale=4; $dN_PSP_LoH / $dS_PSP_LoH" | bc)
    else
        ratio_PSP_LoH="NA"
    fi

    # Intersect to create PSP-rest regions
    bedtools intersect -v -a PSP_gene.bed -b "$loh_file" > "$psp_rest_file"

    # Filter VCF for PSP-rest
    bedtools intersect -wa -a synonymous_missense_variants.vcf -b "$psp_rest_file" > "synonymous_missense_variants_PSP_${isolate}_rest.vcf"

    # Calculate dN and dS for PSP-rest
    dS_PSP_rest=$(grep -c "synonymous" "synonymous_missense_variants_PSP_${isolate}_rest.vcf")
    dN_PSP_rest=$(grep -c "missense" "synonymous_missense_variants_PSP_${isolate}_rest.vcf")

    if [ $dS_PSP_rest -ne 0 ]; then
        ratio_PSP_rest=$(echo "scale=4; $dN_PSP_rest / $dS_PSP_rest" | bc)
    else
        ratio_PSP_rest="NA"
    fi

    # Append results to the output file
    echo -e "${isolate}\t${ratio_PSP_LoH}\t${ratio_PSP_rest}" >> $output_file
done

echo "Calculation complete. Results saved in $output_file."
