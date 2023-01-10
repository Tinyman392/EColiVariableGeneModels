'''
python parseFTP.py [-opt val]

Options:
-f | --ftp : required: location of a copied FTP directory to use.  
-g | --gid : optional: file containing gids to filter from
-o | --out_pref : optional: output prefix for directory/file 
  output.  Defaults to "out"
-n | --n_plf : optional: specify the number of top PLFs to grab per 
  genome.  Defaults to 100
'''

from sys import stderr
import sys
import os
from optparse import OptionParser
from glob import glob

# set PLF conservation hi and low thresholds for accessory
# currently set to get all genes that occur in
#   < 90% of genomes (no conserved in accessory)
#   > 10% of genomes (avoid "one-offs")
PLFTHRESHLO = 0.1
PLFTHRESHHI = 0.9

# output s to stderr
def err(s):
	stderr.write(s)

# takes a directory name and cleans it
def cleanDirNm(dNm):
	if dNm[-1] != '/':
		dNm += '/'

	return dNm

# prases options
# returns set of parsed options
def getOptions():
	parser = OptionParser()

	parser.add_option('-f', '--ftp', help="Required: Location of the FTP genomes directory", metavar="DIR", default='', dest="ftpDir")
	parser.add_option('-g', '--gid', help="Optional: File containing genome ID list to filter off of", metavar='FILE', default='', dest='gidFNm')
	parser.add_option('-o', '--out_pref', help="Optional: Output directory/file prefix to use.  Output will be [out_pref].acc.plf/, [out_pref].con.fasta/, [out_pref].cnts.tab", metavar="STR", default='out', dest='outPref')
	parser.add_option('-n', '--n_plf', help="Optional: specify the number of top PLF to get per genome", metavar="INT", type=int, default=100, dest='nPLF')

	options,args = parser.parse_args()
	options.ftpDir = cleanDirNm(options.ftpDir)

	return options, parser

# given a file name for a file containing genome IDs
# returns the list of genome IDs
# if file name is empty, returns emtpy hash.
def getGIDLst(fNm):
	if fNm == '':
		return {}
	else:
		f = open(fNm)

		gids = {}
		for i in f:
			i = i.strip('\n')
			gids[i] = 0

		f.close()

		return gids

# given options and list of genome IDs
# returns a list of directories for each genome ID
# if genome ID list empty, returns list of all directories
def getFLst(options, gids):
	if len(gids) == 0:
		return glob(options.ftpDir + '*/')
	else:
		dLst = []
		for i in gids:
			dLst.append(options.ftpDir + i + '/')
		return dLst

# given a file stream parses header
# returns a hash that maps header string to index
def getHeader(f):
	ln = f.readline().strip('\n').split('\t')
	headHsh = {}
	for i in range(0,len(ln)):
		headHsh[ln[i]] = i

	return headHsh

# given a directory for a genome, a genome ID, a hash that maps FIG 
# IDs to PLF, hash that maps PLF to length, hash that maps FIG to 
# length, and map that maps GID to a list of FIG
# Returns a PLF hash of counts
def parseGenome(dNm, gid, figPLFHsh, plfLenHsh, figLenHsh, gidFigHsh):
	# set the GID Fig hash for the given GID
	gidFigHsh[gid] = {}

	# open features tabular for the GID
	f = open(dNm + gid + '.PATRIC.features.tab')

	# parse header
	headHsh = getHeader(f)

	# init PLF hash
	pHsh = {}
	for i in f:
		# strip and split line
		i = i.strip('\n').split('\t')
		# get PLF, FIG, and len
		plf = i[headHsh['plfam_id']]
		fig = i[headHsh['patric_id']]
		gLn = int(i[headHsh['na_length']])

		# if no PLF, continue to next line
		if plf == '':
			continue

		# if PLF not in hash, add it and init value to 0
		if plf not in pHsh:
			pHsh[plf] = 0
		# if PLF not in len hash, add it and init empty list
		if plf not in plfLenHsh:
			plfLenHsh[plf] = []
		# set the FIG PLF hash
		figPLFHsh[fig] = plf
		# set the FIG len hash
		figLenHsh[fig] = gLn
		# append len to the PLF len hash
		plfLenHsh[plf].append(gLn)
		# increment count for PLF by 1
		pHsh[plf] += 1

		# set the gid fig hash value
		gidFigHsh[gid][fig] = 0

	f.close()

	# return plf hash
	return pHsh

# given set of options and genome ID list, parses FTP directory
# returns a bunch of hashes
#   PLF hash maps GID to a hash of PLFs
#   fig PLF hash maps fig to PLF
#   PLF len hash maps PLF to len
#   fig len hash maps fig to len
#   gid fig hash maps a gid to figs it has
def parseFTP(options, gids):
	# get list of directories
	dLst = getFLst(options, gids)#[:100]

	# init
	# PLF hash
	# fig PLF hash
	# PLF len hash
	# fig len hash
	# gid fig hash
	plfHsh = {}
	figPLFHsh = {}
	plfLenHsh = {}
	figLenHsh = {}
	gidFigHsh = {}

	# progress bar stuff
	cnt = 0
	inc = len(dLst) / 50.
	err("Parsing feature tabs...\n\t")
	# for each directory
	#   get the GID
	#   parse the genome
	#   set the PLF hash
	for i in dLst:
		# progress bar stuff
		if cnt >= inc:
			err('=')
			cnt = 0
		cnt += 1

		gid = i.split('/')[-2]
		pHsh = parseGenome(i, gid, figPLFHsh, plfLenHsh, figLenHsh, gidFigHsh)
		plfHsh[gid] = pHsh
	err('\n')

	# for each PLF in the len hash
	#   sort the lengths
	#   get the median and average
	#   set the value to median and average
	for i in plfLenHsh:
		plfLenHsh[i].sort()
		med = float(plfLenHsh[i][len(plfLenHsh[i])/2])
		avg = float(sum(plfLenHsh[i])) / len(plfLenHsh[i])

		plfLenHsh[i] = [med, avg]

	# return all hashes
	return plfHsh, figPLFHsh, plfLenHsh, gidFigHsh, figLenHsh

# given a PLF hash
# returns a hash of stats per PLF
def getPLFStaHsh(plfHsh):
	# init the stat hash
	plfStaHsh = {}

	# progress bar stuff
	err("Getting unique PLF list...\n\t")
	cnt = 0
	inc = len(plfHsh) / 50.
	# for each gid in the hash
	#   for each PLF in the gid's hash
	#     if the PLF not i nthe stat hash, init the stat hash
	for i in plfHsh:
		# progress bar stuff
		if cnt >= inc:
			err('=')
			cnt = 0
		cnt += 1

		for j in plfHsh[i]:
			if j not in plfStaHsh:
				plfStaHsh[j] = []
	err('\n')

	# progress bar stuff
	err("Getting stats...\n\t")
	cnt = 0
	inc = len(plfStaHsh) / 50.
	# for each PLF in the stat hash
	#   init counts array
	#   init presence/absence array
	#   for each genome in the hash
	#     check if PLF exists
	#     if so, append 1 for pr/ab hash and count for cnts
	#     otherwise, append 0
	#   compute percent genomes with PLF
	#   sort counts
	#   compute average and median lengths
	#   append stats to stat hash
	for i in plfStaHsh:
		if cnt >= inc:
			err('=')
			cnt = 0
		cnt += 1

		cnts = []
		prab = []
		for j in plfHsh:
			if i in plfHsh[j]:
				cnts.append(plfHsh[j][i])
				prab.append(1)
			else:
				cnts.append(0)
				prab.append(0)

		pGn = float(sum(prab)) / len(prab)
		cnts.sort()
		avg = float(sum(cnts)) / len(cnts)
		med = float(cnts[len(cnts)/2])

		plfStaHsh[i] = [pGn, avg, med]
	err('\n')

	return plfStaHsh

# given options and stat hash, prints the stat hash
def printStats(options, plfStaHsh):
	# open file to print to
	f = open(options.outPref + '.cnts.tab', 'w')

	# for each element in the stat hash sorted by the conservation 
	# of the PLF
	#   init an output array containing PLF and stats
	#   convert array to strings
	#   print to file
	for i in sorted(plfStaHsh, key = lambda x: plfStaHsh[x][0], reverse = True):
		arr = [i] + plfStaHsh[i]
		for j in range(0,len(arr)):
			arr[j] = str(arr[j])
		f.write('\t'.join(arr) + '\n')

	f.close()

# given options and hashes from parseFTP function
# returns a hash of good genomes/GIDs/figs and the top X PLFs
def getGoodGenomes(options, plfStaHsh, gidFigHsh, figPLFHsh, figLenHsh, plfLenHsh):
	# init the top PLFs
	topPLFs = {}
	# for each PLF in the stat hash sorted by % conservation
	#   get the average count per genome
	#   if avg is within 0.01 of 1, then it occurs once per genome
	#     add to top PLF hash
	#   if len of hash > number of top PLFs to get
	#     break
	for i in sorted(plfStaHsh, key = lambda x: plfStaHsh[x][0], reverse = True):
		avg = plfStaHsh[i][1]
		if abs(1 - avg) > 0.01:
			continue
		topPLFs[i] = 0
		if len(topPLFs) >= options.nPLF:
			break

	# initialize the good genome hash
	gGenFighsh = {}

	# progress bar stuff
	err("Getting good GIDs\n\t")
	cnt = 0
	inc = len(gidFigHsh) / 50.
	# for each genome
	#   get the gid
	#   create a genome FIG hash
	#   for each fig in the PLF
	#     get the fig ID
	#     get the PLF for fig
	#     init len to 0
	#     if plf in the top
	#       set len to len of fig
	#     get median length of plf
	#     if len of fig < 0.5x len or > 2x the len
	#       continue (missing or has duplicate PLF)
	#     set the fig for the genome
	#   if genome has all figs in its hash passing the len
	#   cutoffs, add it to the good hash
	for i in gidFigHsh:
		# progress bar stuff
		if cnt >= inc:
			err('=')
			cnt = 0
		cnt += 1

		gid = i
		gFig = {}
		for j in gidFigHsh[i]:
			fig = j
			plf = figPLFHsh[fig]

			l = 0
			if plf in topPLFs:
				l = figLenHsh[fig]

			med = plfLenHsh[plf][0]

			if l < 0.5*med or l > 2.0*med:
				continue

			gFig[fig] = 0

		if len(gFig) == len(topPLFs):
			gGenFighsh[gid] = gFig

	err('\n')

	return gGenFighsh, topPLFs

# given a fasta hash, a gene name and a gene sequence
# adds the gene to the hash if non-empty
def addGeneToHsh(fHsh, cGen, cSeq):
	if cGen != '' and cSeq != '':
		fHsh[cGen] = cSeq

# given options, a gid, and good fig hash, parses the fasta file 
# returns a fasta hash of good conserved genes
def parseFasta(options, gid, gFig):
	# open the fasta file
	f = open(options.ftpDir + gid + '/' + gid + '.PATRIC.ffn')

	# init hash, current gene name, and current gene sequence
	fHsh = {}
	cGen = ''
	cSeq = ''
	# for each line in file
	#   if len of line is 0, continue
	#   if first character is '>', new gene
	#     add the current gene to the hash
	#     set current gene to fig ID
	#     empty out current sequence
	#     if current gene not in good fig IDs
	#       empty out current gene name so it won't be added to hash
	#       on go around
	#     continue to next line
	#   append line to current sequence
	for i in f:
		i = i.strip()
		if len(i) == 0:
			continue
		if i[0] == '>':
			addGeneToHsh(fHsh, cGen, cSeq)
			cGen = '|'.join(i.split('>')[1].split()[0].split('|')[:2])
			cSeq = ''

			# print [cGen]

			if cGen not in gFig:
				cGen = ''
			continue
		cSeq += i

	# add last gene to sequence.  
	addGeneToHsh(fHsh, cGen, cSeq)

	# close file
	f.close()

	return fHsh

# given an output directory, genome ID, and fasta hash
# outputs fasta hash to a fasta formatted file
def writeFasta(oDir, gid, fHsh):
	# open output file
	f = open(oDir + gid + '.fasta', 'w')

	# for each gene in fasta hash
	#   write fasta header
	#   write fasta sequence
	for i in fHsh:
		f.write('>' + i + '\n')
		f.write(fHsh[i] + '\n')

	# close file
	f.close()

# given options and good genome figh hash, writes all fasta
# files for conserved genes
def writeFastas(options, gGenFighsh):
	# set output directory
	# if directory doesn't exist, make it
	oDir = options.outPref + '.con.fasta/'
	if not os.path.exists(oDir):
		os.mkdir(oDir)

	# progress bar stuff
	err("Writing fastas...\n\t")
	cnt = 0
	inc = len(gGenFighsh) / 50.
	# for each good genome
	#   parse the fasta file for the genome
	#   write the fasta file for the genome
	for i in gGenFighsh:
		# progress bar stuff
		if cnt >= inc:
			err('=')
			cnt = 0
		cnt += 1

		fHsh = parseFasta(options, i, gGenFighsh[i])
		writeFasta(oDir, i, fHsh)
	err('\n')

# given a plf stat hash, good genome fig hahs, and plf hash
# table of PLFs for training
def makePLFTabHsh(plfStaHsh, gGenFighsh, plfHsh):
	# init stat table for plfs
	plfTabHsh = {}

	# progress bar stuff
	err("Getting accessory PLF labels...\n\t")
	inc = len(plfStaHsh) / 50.
	cnt = 0
	# for each PLF in stat hash
	#   if PLF outside threshold, continue
	#   if PLF not in hash, add it with additional label 0
	#   else add it with additional label 1
	for i in plfStaHsh:
		if cnt >= inc:
			cnt = 0
			err('=')
		cnt += 1

		if plfStaHsh[i][0] > PLFTHRESHHI or plfStaHsh[i][0] < PLFTHRESHLO:
			continue
		plfTabHsh[i] = []
		for j in gGenFighsh:
			if i not in plfHsh[j]:
				plfTabHsh[i].append([j, '0'])
			else:
				plfTabHsh[i].append([j, '1'])
	err('\n')

	# dLst = []
	# for i in plfTabHsh:
	# 	cnt = [0,0]
	# 	for j in plfTabHsh[i]:
	# 		if j[1] == '1':
	# 			cnt[1] += 1
	# 		if j[1] == '0':
	# 			cnt[0] += 1
	# 	if cnt[0] == 0 and cnt[1] == 0:
	# 		dLst.append(plfTabHsh)

	# for i in dLst:
	# 	del plfTabHsh[i]

	return plfTabHsh

# given output name and table prints table to file
def printPLFTab(oFNm, tab):
	f = open(oFNm, 'w')

	for i in tab:
		f.write('\t'.join(i) + '\n')

	f.close()

# given options and plfTabHsh, outputs presence/absence table per PLF
def printPLFTabHsh(options, plfTabHsh):
	# set output directory
	# if not exists, make directory
	oDir = options.outPref + '.acc.tabs/'
	if not os.path.exists(oDir):
		os.mkdir(oDir)

	# progress bar stuff
	err("Writing accessory PLF tabs...\n\t")
	inc = len(plfTabHsh) / 50.
	cnt = 0
	# for each PLF in table hash
	#   set output file name
	#   print table to file
	for i in plfTabHsh:
		if cnt >= inc:
			cnt = 0
			err('=')
		cnt += 1

		oFNm = oDir + i + '.tab'
		printPLFTab(oFNm, plfTabHsh[i])
	err('\n')

# given options and a list of top PLFs, prints them to a file
def printTopPLFs(options, topPLFs):
	oFNm = options.outPref + '.plf.con.lst'

	f = open(oFNm, 'w')
	for i in topPLFs:
		f.write(i + '\n')
	f.close()

# main driver program
# 
# parse options
# get GIDs
#
# parse the FTP to ge ta bunch of hashes
#
# get the stats for PLFs
# print the stats to file
# 
# get good genomes
# 
# print top PLFs to file
#
# write conserved fasta files
#
# make the PLF tab hash for accessories to train on
# print the tabular files
def main():
	options, parser = getOptions()
	gids = getGIDLst(options.gidFNm)

	plfHsh,figPLFHsh,plfLenHsh,gidFigHsh, figLenHsh = parseFTP(options, gids)

	plfStaHsh = getPLFStaHsh(plfHsh)
	printStats(options, plfStaHsh)

	gGenFighsh, topPLFs = getGoodGenomes(options, plfStaHsh, gidFigHsh, figPLFHsh, figLenHsh, plfLenHsh)

	printTopPLFs(options, topPLFs)

	writeFastas(options, gGenFighsh)

	plfTabHsh = makePLFTabHsh(plfStaHsh, gGenFighsh, plfHsh)
	printPLFTabHsh(options, plfTabHsh)
	
if __name__ == '__main__':
	main()

