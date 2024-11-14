#!/bin/sh
## Give a name to  your job
#SBATCH --job-name=bwa_lumpy
## precise the logfile for your job
#SBATCH --output=bwa_lumpy.out
## precise the error file for your job
#SBATCH --error=bwa_lumpy_err.out
# Precise the partion you want to use
#SBATCH --partition=normal
# precise when you receive the email
#SBATCH --mail-type=end
# precise to  which address you have to send the mail to
#SBATCH --mail-user=ngan.phan-thi@ird.fr
# number of cpu you want to use on you node
#SBATCH --cpus-per-task=2

############################################################

path_to_dir="/data3/projects/graminicola_evolution/ISOLATES";
path_to_tmp="/scratch/bwa_lumpy";
path_to_dest="/data3/projects/graminicola_evolution/ISOLATES/CNV";

############# chargement du module load CNVnator
module load bioinfo/lumpy/0.2.13
module load bioinfo/samtools/1.10
module load bioinfo/samblaster/0.1.20
module load bioinfo/bwa/0.7.17

###### Creation du repertoire temporaire sur  la partition /scratch du noeud
mkdir $path_to_tmp

####### copie du repertoire de données  vers la partition /scratch du noeud
scp nas3:$path_to_dir/CANU_283scf_renamed.fasta $path_to_tmp/


echo "tranfert donnees master -> noeud"
cd $path_to_tmp/

###### Execution du programme

scp nas3:$path_to_dir/FASTQ/MgBl* $path_to_tmp/
bwa index CANU_283scf_renamed.fasta

for f in *_R1.fastq.gz; do sample=${f%%_*};

bwa mem -R "@RG\tID:id\tSM:MgBl\tLB:lib" CANU_283scf_renamed.fasta ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 | samtools view -S -b -o ${sample}.bam ;
# Extract the discordant paired-end alignments.
samtools view -b -F 1294 ${sample}.bam -o ${sample}.discordants.unsorted.bam ;
# Extract the split-read alignments
samtools view -h ${sample}.bam | /usr/local/lumpy-0.2.13/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb -o ${sample}.splitters.unsorted.bam ;

# Sort both alignments
samtools sort ${sample}.discordants.unsorted.bam -o ${sample}.discordants ;
samtools sort ${sample}.splitters.unsorted.bam -o ${sample}.splitters ;

done; 






rm *fastq.gz





###nsfert des données du noeud vers master

scp -rp $path_to_tmp  nas3:$path_to_dest/;
echo "Transfert donnees node -> master";

#### Suppression du repertoire tmp noeud
rm -rf $path_to_tmp;
echo "Suppression des donnees sur le noeud";

