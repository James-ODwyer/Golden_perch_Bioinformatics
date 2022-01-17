#!/bin/bash
#SBATCH 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --time=16:00:00
#SBATCH --mail-user=18088076@students.latrobe.edu.au
#SBATCH --ntasks=1


echo "Starting at: $(date)"

cd /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/Kraken

for file in *output.output

do
echo " $file "

rm $file


done

echo "Finished at: $(date)"
