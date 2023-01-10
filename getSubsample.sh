#! /bin/bash
#
# bash getSubsample.sh [plf.con.fasta dir] [plf.con.kmc dir] 
# [clust file]
# 
# reads a subsampled list of genome IDs and creates directories for
# that subsample and copies files over for subsample.  

# clean up fasta and kmc directory names
fIn=$(echo $1 | sed 's/\/$//g')
kIn=$(echo $2 | sed 's/\/$//g')

# get the number of genomes
n=$(basename $3 | cut -f3 -d'.')
# get the output fasta and kmc directories
fOut=$fIn.$n
kOut=$kIn.$n

# if output directories don't exist, make them
if [ ! -d $fOut ]; then
	mkdir $fOut
fi
if [ ! -d $kOut ]; then
	mkdir $kOut
fi

# for each gid in the cluster file
#   cp the fasta into the output fasta output directory
#   cp the KMC file into the KMC output directory
for i in $(cut -f1 -d$'\t' $3); do
	echo $i
	cp $fIn/$i.fasta $fOut
	cp $kIn/$i.fasta.*.kmrs $kOut
done