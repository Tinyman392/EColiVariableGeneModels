#!/bin/bash
#
# Downloads the relevant ftp data from the PATRIC FTP.  Places a 
# directory named "ftp" in the working directory.  
#
# IF "ftp" DIRECTORY EXISTS IN THE WORKING DIRECTORY, IT WILL BE 
# REMOVED AND REMADE!!!
#
# While downloading genomes, the script will only attempt to 
# download each file once upon failure.  In case, at a later date 
# the genomes are removed for any reason, no infinite loop will 
# exist.

# get script directory
SDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# make ftp directory, if it doesn't exist, make it
# otherwise clear it out then make it
if [ ! -d ftp/ ]; then
	mkdir ftp/
else
	rm -r ftp/
	mkdir ftp/
fi

# make directories for the RELEASE_NOTES and genomes directories
mkdir ftp/RELEASE_NOTES ftp/genomes

# copy the genome summary into the RELEASE_NOTES directory
curl ftp:/ftp.bvbrc.org/RELEASE_NOTES/genome_summary > ftp/RELEASE_NOTES/genome_summary 2> /dev/null

# for each genome ID in the genomes.gid.lst file
#   make a directory for the genome
#   try to download the PATRIC.faa, PATRIC.features.tab, PATRIC.ffn, 
#     and fna files.  Upon failure, attempt redownload once.  
for i in $(cat $SDIR/genomes.gid.lst); do
	echo $i
	mkdir ftp/genomes/$i/
	for j in PATRIC.faa PATRIC.features.tab PATRIC.ffn fna; do
		curl ftp:/ftp.bvbrc.org/genomes/$i/$i.$j > ftp/genomes/$i/$i.$j 2> /dev/null
		if [ ! -f ftp/genomes/$i/$i.$j ] || [ -s ftp/genomes/$i/$i.$j ]; then
			curl ftp:/ftp.bvbrc.org/genomes/$i/$i.$j > ftp/genomes/$i/$i.$j 2> /dev/null
		fi
	done
done