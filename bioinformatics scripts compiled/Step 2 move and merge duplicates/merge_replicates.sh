#!/bin/bash
#SBATCH 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --time=3:00:00
#SBATCH --mail-user=18088076@students.latrobe.edu.au
#SBATCH --ntasks=1

echo "Starting at: $(date)"

cd /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/process_radtags

mkdir -p merged
 
for f in *.1.fq.gz
do
base=${f%%.1.fq.gz}
cat ${base}*.gz > merged/${base}_merged.fq.gz
done

echo "Finished at: $(date)"
