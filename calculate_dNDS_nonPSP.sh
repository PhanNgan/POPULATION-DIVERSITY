#!/bin/bash

# Output file for the results
output_file="dn_ds_ratios_nonPSP.tsv"
echo -e "Isolate\dS_nonPSP_LoH\dN_nonPSP_LoH\ratio_nonPSP_LoH\dS_nonPSP_rest\dN_nonPSP_rest\ratio_nonPSP_rest" > $output_file

# List of isolate prefixes
isolates=(MgBl MgBn MgBr MgC21 MgC25 MgJv2 MgL1 MgL2 MgP MgVN11 MgVN18 MgVN27 MgVN6)

# Loop through each isolate
for isolate in "${isolates[@]}"; do
    # Define input files
    loh_file="${isolate}_loh.merged_50_loh.merged_concatenated.bed"
    nonpsp_loh_file="${isolate}_nonPSP_LoH.bed"
    nonpsp_rest_file="${isolate}_nonPSP-rest.bed"

    # Intersect to create PSP_LoH regions
    bedtools intersect -wb -a "$loh_file" -b non_PSP_14098gene.bed > "$nonpsp_loh_file"

    # Filter VCF for PSP_LoH
    bedtools intersect -wa -a synonymous_missense_variants.vcf -b "$nonpsp_loh_file" > "synonymous_missense_variants_nonPSP_${isolate}_LoH.vcf"

    # Calculate dN and dS for PSP_LoH
    dS_nonPSP_LoH=$(grep -c "synonymous" "synonymous_missense_variants_nonPSP_${isolate}_LoH.vcf")
    dN_nonPSP_LoH=$(grep -c "missense" "synonymous_missense_variants_nonPSP_${isolate}_LoH.vcf")

    if [ $dS_nonPSP_LoH -ne 0 ]; then
        ratio_nonPSP_LoH=$(echo "scale=4; $dN_nonPSP_LoH / $dS_nonPSP_LoH" | bc)
    else
        ratio_nonPSP_LoH="NA"
    fi

    # Intersect to create PSP-rest regions
    bedtools intersect -v -a non_PSP_14098gene.bed -b "$loh_file" > "$nonpsp_rest_file"

    # Filter VCF for PSP-rest
    bedtools intersect -wa -a synonymous_missense_variants.vcf -b "$nonpsp_rest_file" > "synonymous_missense_variants_nonPSP_${isolate}_rest.vcf"

    # Calculate dN and dS for PSP-rest
    dS_nonPSP_rest=$(grep -c "synonymous" "synonymous_missense_variants_nonPSP_${isolate}_rest.vcf")
    dN_nonPSP_rest=$(grep -c "missense" "synonymous_missense_variants_nonPSP_${isolate}_rest.vcf")

    if [ $dS_nonPSP_rest -ne 0 ]; then
        ratio_nonPSP_rest=$(echo "scale=4; $dN_nonPSP_rest / $dS_nonPSP_rest" | bc)
    else
        ratio_nonPSP_rest="NA"
    fi

    # Append results to the output file
    echo -e "${isolate}\t${dS_nonPSP_LoH}\t${dN_nonPSP_LoH}\t${ratio_nonPSP_LoH}\t${dS_nonPSP_rest}\t${dN_nonPSP_rest}\t${ratio_nonPSP_rest}" >> $output_file
done

echo "Calculation complete. Results saved in $output_file."
