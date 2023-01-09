#! /bin/bash
#
# bash getSubsample.sh [plf.con.fasta dir] [plf.con.kmc dir] [clust file]

fIn=$(echo $1 | sed 's/\/$//g')
kIn=$(echo $2 | sed 's/\/$//g')

n=$(basename $3 | cut -f3 -d'.')
fOut=$fIn.$n
kOut=$kIn.$n

if [ ! -d $fOut ]; then
	mkdir $fOut
fi
if [ ! -d $kOut ]; then
	mkdir $kOut
fi

for i in $(cut -f1 -d$'\t' $3); do
	echo $i
	cp $fIn/$i.fasta $fOut
	cp $kIn/$i.fasta.*.kmrs $kOut
done