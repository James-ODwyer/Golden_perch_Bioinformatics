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


vcftools --vcf GP_filtered_step3_call_loc60_ind50_maf_singleton.vcf --thin 16000 --recode --recode-INFO-all --out vcf_GP_filtered_lax_complete

echo "Finishing at: $(date)"