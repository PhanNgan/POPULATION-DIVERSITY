#!/bin/sh
## Give a name to  your job
#SBATCH --job-name=CNVnator_ex
## precise the logfile for your job
#SBATCH --output=CNVnator_ex.out
## precise the error file for your job
#SBATCH --error=CNVnator_err.out
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
path_to_tmp="/scratch/CNVnator";
path_to_dest="/data3/projects/graminicola_evolution/ISOLATES/CNV";

############# chargement du module load CNVnator
module load bioinfo/CNVnator/0.4.1 

###### Creation du repertoire temporaire sur  la partition /scratch du noeud
mkdir $path_to_tmp

####### copie du repertoire de données  vers la partition /scratch du noeud
scp -r nas3:$path_to_dir/isolates_canu/Isolates_canu/output/MgBl/8_gatkIndelRealigner/MgBl.GATKINDELREALIGNER.bam $path_to_tmp/





echo "tranfert donnees master -> noeud";
cd $path_to_tmp/;

###### Execution du programme
cnvnator -root MgBl.root -tree MgBl.GATKINDELREALIGNER.bam 
cnvnator -root MgBl.root -his 100 -f .
cnvnator -root MgBl.root -stat 100
cnvnator -root MgBl.root -partition 100 
cnvnator -root MgBl.root -call 100




##### Transfert des données du noeud vers master
#rm *bam
scp -rp $path_to_tmp  nas3:$path_to_dest/;
echo "Transfert donnees node -> master";

#### Suppression du repertoire tmp noeud
#rm -rf $path_to_tmp;
echo "Suppression des donnees sur le noeud";

