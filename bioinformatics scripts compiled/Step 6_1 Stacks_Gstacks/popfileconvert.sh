#!/bin/bash
#SBATCH 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --time=1:00:00
#SBATCH --mail-user=18088076@students.latrobe.edu.au
#SBATCH --ntasks=1



sed 's/,/\t/g' poplist_GP.csv > poplist_GP.txt