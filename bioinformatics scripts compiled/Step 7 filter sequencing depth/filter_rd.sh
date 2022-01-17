#!/bin/bash
#SBATCH 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --partition=compute
#SBATCH --time=2:00:00
#SBATCH --mail-user=18088076@students.latrobe.edu.au
#SBATCH --ntasks=1


module load vcftools-gcc/0.1.13

echo "Starting at: $(date)"


vcftools --vcf populations.snps.vcf --min-meanDP 10 --minDP 5 --recode --recode-INFO-all --out vcf_filtered_MAFDP10minDP5


echo "Finishing at: $(date)"