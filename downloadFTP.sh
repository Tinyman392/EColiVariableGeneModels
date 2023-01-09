#!/bin/bash
#
# Downloads the relevant ftp data from the PATRIC FTP.  Places a directory named "ftp" in the working directory.  
#
# IF "ftp" DIRECTORY EXISTS IN THE WORKING DIRECTORY, IT WILL BE REMOVED AND REMADE!!!

SDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [ ! -d ftp/ ]; then
	mkdir ftp/
else
	rm -r ftp/
	mkdir ftp/
fi

mkdir ftp/RELEASE_NOTES ftp/genomes

curl ftp:/ftp.bvbrc.org/RELEASE_NOTES/genome_summary > ftp/RELEASE_NOTES/genome_summary 2> /dev/null

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