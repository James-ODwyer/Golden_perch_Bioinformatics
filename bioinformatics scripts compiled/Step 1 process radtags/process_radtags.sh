#!/bin/bash
#SBATCH 
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --time=48:00:00
#SBATCH --mail-user=18088076@students.latrobe.edu.au
#SBATCH --ntasks=1

echo "Starting at: $(date)"
cd /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/fastp_a/

module load stacks-gcc7/2.41
 
for file in *.FASTQ.gz
do
    	id=${file%.FASTQ.gz}
    	process_radtags -f  $file \
    	-b \
/data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/fastp_a/single_barcodes/$id.txt \
    	-o /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/process_radtags \
-e pstI --retain_header \
    	-c -q -r
    	cp \
/data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/process_radtags/process_radtags.log /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/process_radtags/process_radtags_$id.log
done
 
echo "Finished at: $(date)"
