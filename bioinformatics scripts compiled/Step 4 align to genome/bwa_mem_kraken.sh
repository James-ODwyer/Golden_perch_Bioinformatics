#!/bin/bash
#SBATCH 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000
#SBATCH --partition=compute
#SBATCH --time=149:00:00
#SBATCH --mail-user=18088076@students.latrobe.edu.au
#SBATCH --ntasks=1


module load bwa-gcc/0.7.17 

echo "Starting at: $(date)"
 
cd /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/bwamem_GP/
 
for file in *unclassified.kraken

do

file2=${file%unclassified.kraken}

    bwa mem -t 16 \
    /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/GP_genome/GOP001_1.1_HiC.fasta \
     $file > /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/bwamem_GP/$file2.sam


echo "Finished ind $file "

echo " $file2 "

#rm $file

done

 
echo "Finished at: $(date)"
