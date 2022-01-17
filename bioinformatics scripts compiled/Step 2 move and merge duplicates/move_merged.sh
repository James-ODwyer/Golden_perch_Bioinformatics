#!/bin/bash
#SBATCH 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --time=1:00:00
#SBATCH --mail-user=18088076@students.latrobe.edu.au
#SBATCH --ntasks=1


cd /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/process_radtags
  
mkdir duplicated

for f in *.1.fq.gz
do
base=${f%%.1.fq.gz}
mv ${base}*.gz duplicated/
done

cd /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/process_radtags/merged

for f in *
do
mv $f ../${f%%_merged.fq.gz}.fq.gz

done