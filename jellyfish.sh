#!/bin/sh
## Give a name to  your job
#SBATCH --job-name=jelly_fish_mg
## precise the logfile for your job
#SBATCH --output=jellyfish_mg.out
## precise the error file for your job
#SBATCH --error=jellyfish_mg.er
# Precise the partion you want to use
#SBATCH --partition=normal
# precise when you receive the email
#SBATCH --mail-type=end
# precise to  which address you have to send the mail to
#SBATCH --mail-user=thi-ngan.phan@ird.fr
# number of cpu you want to use on you node
#SBATCH --cpus-per-task=2

############################################################

path_to_dir="/projects/medium/graminicola_evolution/ISOLATES/3.cleaned_FASTQ/FASTQ";
path_to_tmp="/scratch/mg_kmer1";
path_to_dest="/projects/medium/graminicola_evolution/ISOLATES";

############# chargement du module load toggle
mkdir $path_to_tmp
cd $path_to_tmp

module load jellyfish/2.3.0

scp $path_to_dir/*fastq.gz $path_to_tmp
gunzip *.gz

# Define the list of genome names
genomes=("MgBl" "MgBr" "MgC25" "MgL1" "MgP" "MgVN18" "MgVN6" "MgBn" "MgC21" "MgJv2" "MgL2" "MgVN11" "MgVN27")

# Loop through each genome
for genome in "${genomes[@]}"; do
    echo "Processing genome: $genome"
    
    # Loop through k-mer sizes (17, 21, 27, 47)
    for kmer in 17 21 27 47; do
        echo "Counting kmers for $genome with k-mer size $kmer"
        
        # Count kmers using Jellyfish for the given k-mer size
        jellyfish count -C -m $kmer -s 1000000000 -t 10 ${genome}_R*.fastq -o ${genome}_k${kmer}_reads.jf
        
        # Export the k-mer count histogram
        jellyfish histo -h 1000000 -t 10 ${genome}_k${kmer}_reads.jf > ${genome}_k${kmer}_reads.histo
    done
done

echo "K-mer counting completed for all genomes."
rm *fastq*
scp -r $path_to_tmp $path_to_dest
rm -rp $path_to_tmp
