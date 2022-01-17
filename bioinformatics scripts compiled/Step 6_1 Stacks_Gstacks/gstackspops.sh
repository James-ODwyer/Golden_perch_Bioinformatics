#!/bin/bash
#SBATCH 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --time=140:00:00
#SBATCH --mail-user=18088076@students.latrobe.edu.au
#SBATCH --ntasks=1


echo "Starting at: $(date)"


module load stacks-gcc7/2.41

cd /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/stacks_gstacks_populations/
 
gstacks \
-I /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/samtools_viewsort_GP/ \
-M /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/stacks_gstacks_populations/poplist_GP.txt \
-O /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/stacks_gstacks_populations/ \
-t 16 --min-mapq 30

echo "Finishing at: $(date)"
