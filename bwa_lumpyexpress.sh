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

path_to_dir="/data3/projects/graminicola_evolution/ISOLATES"
path_to_tmp="/scratch/bwa_lumpy"
path_to_dest="/data3/projects/graminicola_evolution/ISOLATES/CNV"

############# chargement du module load CNVnator
module load bioinfo/lumpy/0.2.13
module load bioinfo/samtools/1.10
module load bioinfo/samblaster/0.1.20
module load bioinfo/bwa/0.7.17
module load bioinfo/lumpyexpress/0.2.13

###### Creation du repertoire temporaire sur  la partition /scratch du noeud
mkdir $path_to_tmp

####### copie du repertoire de données  vers la partition /scratch du noeud
scp nas3:$path_to_dir/CANU_283scf_renamed.fasta $path_to_tmp/
scp nas3:$path_to_dir/FASTQ/*fastq.gz $path_to_tmp/
echo "tranfert donnees master -> noeud"
cd $path_to_tmp/

###### Execution du programme

bwa index CANU_283scf_renamed.fasta

for f in *_R1.fastq.gz; do sample=`basename ${f%%_*}`; 
  if [ ! -f ${sample}.splitters.bam ]; then 
     echo "processing sample ${sample}"; 

	bwa mem -R "@RG\tID:id\tSM:${sample}\tLB:lib" CANU_283scf_renamed.fasta ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 | samtools view -S -b -o ${sample}.bam ;
   # Extract the discordant paired-end alignments.
	samtools view -b -F 1294 ${sample}.bam -o ${sample}.discordants.unsorted.bam ;
   # Extract the split-read alignments
	samtools view -h ${sample}.bam | /usr/local/lumpy-0.2.13/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb -o ${sample}.splitters.unsorted.bam ;

   # Sort both alignments
	samtools sort ${sample}.discordants.unsorted.bam -o ${sample}.discordants.bam ;
	samtools sort ${sample}.splitters.unsorted.bam -o ${sample}.splitters.bam ;  
  fi;
done


lumpyexpress -B MgBl.bam,MgBn.bam,MgBr.bam,MgC21.bam,MgC25.bam,MgJv2.bam,MgL1.bam,MgL2.bam,MgP.bam,MgVN11.bam,MgVN18.bam,MgVN27.bam,MgVN6.bam -S MgBl.splitters.bam,MgBn.splitters.bam,MgBr.splitters.bam,MgC21.splitters.bam,MgC25.splitters.bam,MgJv2.splitters.bam,MgL1.splitters.bam,MgL2.splitters.bam,MgP.splitters.bam,MgVN11.splitters.bam,MgVN18.splitters.bam,MgVN27.splitters.bam,MgVN6.splitters.bam -D MgBl.discordants.bam,MgBn.discordants.bam,MgBr.discordants.bam,MgC21.discordants.bam,MgC25.discordants.bam,MgJv2.discordants.bam,MgL1.discordants.bam,MgL2.discordants.bam,MgP.discordants.bam,MgVN11.discordants.bam,MgVN18.discordants.bam,MgVN27.discordants.bam,MgVN6.discordants.bam -o lumpyexpress_multi_sample.vcf


rm *fastq.gz
rm *fasta
###nsfert des donnes du noeud vers master

scp -rp $path_to_tmp  nas3:$path_to_dest/
echo "Transfert donnees node -> master"

#### Suppression du repertoire tmp noeud
#rm -rf $path_to_tmp
echo "Suppression des donnees sur le noeud"




