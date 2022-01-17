#!/bin/bash
#SBATCH 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --time=48:00:00
#SBATCH --mail-user=18088076@students.latrobe.edu.au
#SBATCH --ntasks=1

echo "Starting at: $(date)"
cd /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/Kraken/

module load kraken-gcc/0.10.4

#kraken-build --download-taxonomy --db bacteria 

kraken-build --download-library bacteria --db bacteria 



#--threads 16 --max-db-size 300