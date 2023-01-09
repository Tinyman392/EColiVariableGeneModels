#/bin/bash
#
# runKMC.sh [fasta dir] [out dir]

FDIR=$(echo $1 | sed 's/\/$//g')
ODIR=$(echo $2 | sed 's/\/$//g')

if [ ! -d $ODIR ]; then
	mkdir $ODIR
fi

for i in $FDIR/*.fasta; do
	echo $i
	oFil=$ODIR/$(basename $i)
	kmc.sh 7 $i $oFil ./ 12 > /dev/null
	rm $oFil.kmc_pre $oFil.kmc_suf
done
