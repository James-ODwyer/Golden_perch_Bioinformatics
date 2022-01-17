#!/bin/bash
#SBATCH 
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --time=10:00:00
#SBATCH --mail-user=18088076@students.latrobe.edu.au
#SBATCH --ntasks=1


cd /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/stacks_gstacks_populations/

echo "Starting at: $(date)"

mkdir ./p20r7

module load stacks-gcc7/2.41

populations \
-P /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/stacks_gstacks_populations/ \
-M /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/stacks_gstacks_populations/poplist_GP.txt \
-O /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/stacks_gstacks_populations/p20r7 \
-t 16 -p 20 -r 0.7 --write_single_snp --hwe --fstats --ordered_export --fasta_samples --fasta_loci --plink --vcf --vcf_haplotypes --genepop --structure --phylip --phylip_var \
--fasta_samples_raw

echo "Finishing at: $(date)"
