#!/bin/bash
#SBATCH 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --time=120:00:00
#SBATCH --mail-user=18088076@students.latrobe.edu.au
#SBATCH --ntasks=1


echo "Starting at: $(date)"
cd /data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/Kraken/
module load kraken-gcc/0.10.4

TAXONOMY_DIR="$/data/group/murphylab/project/James/GP_bioinformatics/DMacq18-3638/Kraken/bacteria/taxonomy"
NCBI_SERVER="ftp.ncbi.nih.gov"
FTP_SERVER="ftp://$NCBI_SERVER"
THIS_DIR=$PWD

cd "$TAXONOMY_DIR"

  wget $FTP_SERVER/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
  wget $FTP_SERVER/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
  touch accmap.dlflag
  echo "Downloaded accession to taxon map(s)"

wget $FTP_SERVER/pub/taxonomy/taxdump.tar.gz
  touch taxdump.dlflag
  echo "Downloaded taxonomy tree data"
  
if ls | grep -q 'accession2taxid\.gz$'
then
  echo -n "Uncompressing taxonomy data... "
  gunzip *accession2taxid.gz
  echo "done."
fi


gunzip gi_taxid_nucl.dmp.gz
  touch gimap.flag
  echo "Uncompressed GI to taxon map"

tar zxf taxdump.tar.gz
  touch taxdump.flag
  echo "Uncompressed taxonomy tree data"

