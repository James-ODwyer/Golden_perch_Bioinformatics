#!/bin/bash
#SBATCH 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --time=100:00:00
#SBATCH --mail-user=18088076@students.latrobe.edu.au
#SBATCH --ntasks=1


echo "Starting at: $(date)"

cd /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/bwamem_GP

module load samtools-gcc/1.9

for f in *.sam
do
    samtools view -uhq 20 -@ 16 \
    /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/bwamem_GP/$f | \
    samtools sort -m 1G -@ 16 \
   -o \
    /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/samtools_viewsort_GP/${f%%.sam}.bam

rm $f

done
echo "Finishing at: $(date)"
