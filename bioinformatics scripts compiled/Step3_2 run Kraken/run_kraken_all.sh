#!/bin/bash
#SBATCH 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --time=149:00:00
#SBATCH --mail-user=18088076@students.latrobe.edu.au
#SBATCH --ntasks=1


echo "Starting at: $(date)"

module load kraken-gcc/0.10.4

cd /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/process_radtags/


for file in *.gz

do
echo "$file"

cp "$file" --target-directory=/data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/Kraken

a="classified.kraken"
b="unclassified.kraken"
c="output.output"
rep="taxonreport"
file2=${file%.fq.gz}
d=$file2$a
e=$file2$b
f=$file2$c
g=$file2$rep

#echo " $file2 $d   $e   $f   $g "


cd /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/Kraken/

kraken --db bacteria $file --threads 16 --classified-out $d --unclassified-out $e > $f

echo " "$file2" kraken run complete"

kraken-report --db bacteria $f > $g.txt

echo " "$file2" kraken report complete"

rm $file

cd /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/process_radtags/

done