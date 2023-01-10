#! /bin/bash
#
# bash trainPipeline.sh [/PATH/TO/GenomicModelCreator/] <nthread>
#
# This script will run the entire prediction pipeline to train the 
# model.  It will only train the first model unless the break 
# statement is commented out
#

##
# Setup and read parameters
##
SDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
CDIR=$(echo $1 | sed 's/\/$//g')

NTHREAD=1
if [ $# -gt 2 ]; then
	NTHREAD=$2
fi

##
# downloadFTP.sh
##

echo "Downloading FTP"
bash $SDIR/downloadFTP.sh

##
# parseFTP.py
##

echo "Parsing FTP"
python $SDIR/parseFTP.py -f ftp/genomes -o all.100plf

##
# runKMC.sh
##

echo "Running KMC"
bash $SDIR/runKMC.sh all.100plf.con.fasta/ all.100plf.con.kmc/

##
# Step: run KMC on all genomes
##

echo "Catting fastas and running KMC"
cat all.100plf.con.fasta/* > allFasta.fasta
kmc.sh 7 allFasta.fasta allFasta.fasta . 1

##
# getClusters.py
##

echo "Getting clusters"
python $SDIR/getClusters.py all.100plf.con.kmc/ allFasta.fasta.7.kmrs all.100plf

##
# getSubsample.sh
##

echo "Subsampling"
bash $SDIR/getSubsample.sh all.100plf.con.fasta/ all.100plf.con.kmc all.100plf.500.clusts

##
# buildModel.py
##

echo "Training"
mkdir all.100plf.500.models

for i in all.100plf.acc.tabs/*.tab; do
	python $CDIR/buildModel.py -f all.100plf.con.fasta.500/ -t $i -T temp/ -o all.100plf.500.models/$(basename $i) -n $NTHREAD -d 16 -k 7 -K all.100plf.con.kmc.500/ -c True -w True -W False
	break
done
