#/bin/bash
#
# runKMC.sh [fasta dir] [out dir]
#
# runs KMC for all fasta files in the fasta directory
# writes output to the output directory

# clean Fasta and output directory names
FDIR=$(echo $1 | sed 's/\/$//g')
ODIR=$(echo $2 | sed 's/\/$//g')

# if output directory doesn't exist, make it
if [ ! -d $ODIR ]; then
	mkdir $ODIR
fi

# for each fasta
#   set the output file name
#   run KMC
#   remove the intermediate KMC files
for i in $FDIR/*.fasta; do
	echo $i
	oFil=$ODIR/$(basename $i)
	kmc.sh 7 $i $oFil ./ 12 > /dev/null
	rm $oFil.kmc_pre $oFil.kmc_suf
done
