'''
python parseFTP.py [-opt val]
'''

from sys import stderr
import sys
import os
from optparse import OptionParser
from glob import glob

PLFTHRESHLO = 0.1
PLFTHRESHHI = 0.9

def err(s):
	stderr.write(s)

def cleanDirNm(dNm):
	if dNm[-1] != '/':
		dNm += '/'

	return dNm

def getOptions():
	parser = OptionParser()

	parser.add_option('-f', '--ftp', help="Required: Location of the FTP genomes directory", metavar="DIR", default='', dest="ftpDir")
	parser.add_option('-g', '--gid', help="Optional: File containing genome ID list to filter off of", metavar='FILE', default='', dest='gidFNm')
	parser.add_option('-o', '--out_pref', help="Optional: Output directory/file prefix to use.  Output will be [out_pref].acc.plf/, [out_pref].con.fasta/, [out_pref].cnts.tab", metavar="STR", default='out', dest='outPref')
	parser.add_option('-n', '--n_plf', help="Optional: specify the number of top PLF to get per genome", metavar="INT", type=int, default=100, dest='nPLF')

	options,args = parser.parse_args()
	options.ftpDir = cleanDirNm(options.ftpDir)

	return options, parser

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

def getFLst(options, gids):
	if len(gids) == 0:
		return glob(options.ftpDir + '*/')
	else:
		dLst = []
		for i in gids:
			dLst.append(options.ftpDir + i + '/')
		return dLst

def getHeader(f):
	ln = f.readline().strip('\n').split('\t')
	headHsh = {}
	for i in range(0,len(ln)):
		headHsh[ln[i]] = i

	return headHsh

def parseGenome(dNm, gid, figPLFHsh, plfLenHsh, figLenHsh, gidFigHsh):
	gidFigHsh[gid] = {}

	f = open(dNm + gid + '.PATRIC.features.tab')

	headHsh = getHeader(f)

	pHsh = {}
	for i in f:
		i = i.strip('\n').split('\t')
		plf = i[headHsh['plfam_id']]
		fig = i[headHsh['patric_id']]
		gLn = int(i[headHsh['na_length']])

		if plf == '':
			continue

		if plf not in pHsh:
			pHsh[plf] = 0
		if plf not in plfLenHsh:
			plfLenHsh[plf] = []
		figPLFHsh[fig] = plf
		figLenHsh[fig] = gLn
		plfLenHsh[plf].append(gLn)
		pHsh[plf] += 1

		gidFigHsh[gid][fig] = 0

	f.close()

	return pHsh

def parseFTP(options, gids):
	dLst = getFLst(options, gids)#[:100]

	plfHsh = {}
	figPLFHsh = {}
	plfLenHsh = {}
	figLenHsh = {}
	gidFigHsh = {}

	cnt = 0
	inc = len(dLst) / 50.
	err("Parsing feature tabs...\n\t")
	for i in dLst:
		if cnt >= inc:
			err('=')
			cnt = 0
		cnt += 1

		gid = i.split('/')[-2]
		pHsh = parseGenome(i, gid, figPLFHsh, plfLenHsh, figLenHsh, gidFigHsh)
		plfHsh[gid] = pHsh
	err('\n')

	for i in plfLenHsh:
		plfLenHsh[i].sort()
		med = float(plfLenHsh[i][len(plfLenHsh[i])/2])
		avg = float(sum(plfLenHsh[i])) / len(plfLenHsh[i])

		plfLenHsh[i] = [med, avg]

	return plfHsh, figPLFHsh, plfLenHsh, gidFigHsh, figLenHsh

def getPLFStaHsh(plfHsh):
	plfStaHsh = {}
	err("Getting unique PLF list...\n\t")
	cnt = 0
	inc = len(plfHsh) / 50.
	for i in plfHsh:
		if cnt >= inc:
			err('=')
			cnt = 0
		cnt += 1

		for j in plfHsh[i]:
			if j not in plfStaHsh:
				plfStaHsh[j] = []
	err('\n')

	err("Getting stats...\n\t")
	cnt = 0
	inc = len(plfStaHsh) / 50.
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

def printStats(options, plfStaHsh):
	f = open(options.outPref + '.cnts.tab', 'w')

	for i in sorted(plfStaHsh, key = lambda x: plfStaHsh[x][0], reverse = True):
		arr = [i] + plfStaHsh[i]
		for j in range(0,len(arr)):
			arr[j] = str(arr[j])
		f.write('\t'.join(arr) + '\n')

	f.close()

def getGoodGenomes(options, plfStaHsh, gidFigHsh, figPLFHsh, figLenHsh, plfLenHsh):
	topPLFs = {}
	for i in sorted(plfStaHsh, key = lambda x: plfStaHsh[x][0], reverse = True):
		avg = plfStaHsh[i][1]
		if abs(1 - avg) > 0.01:
			continue
		topPLFs[i] = 0
		if len(topPLFs) >= options.nPLF:
			break

	gGenFighsh = {}

	err("Getting good GIDs\n\t")
	cnt = 0
	inc = len(gidFigHsh) / 50.
	for i in gidFigHsh:
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

def addGeneToHsh(fHsh, cGen, cSeq):
	if cGen != '' and cSeq != '':
		fHsh[cGen] = cSeq

def parseFasta(options, gid, gFig):
	f = open(options.ftpDir + gid + '/' + gid + '.PATRIC.ffn')

	fHsh = {}
	cGen = ''
	cSeq = ''
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

	addGeneToHsh(fHsh, cGen, cSeq)

	f.close()

	return fHsh

def writeFasta(oDir, gid, fHsh):
	f = open(oDir + gid + '.fasta', 'w')

	for i in fHsh:
		f.write('>' + i + '\n')
		f.write(fHsh[i] + '\n')

	f.close()

def writeFastas(options, gGenFighsh):
	oDir = options.outPref + '.con.fasta/'
	if not os.path.exists(oDir):
		os.mkdir(oDir)

	err("Writing fastas...\n\t")
	cnt = 0
	inc = len(gGenFighsh) / 50.
	for i in gGenFighsh:
		if cnt >= inc:
			err('=')
			cnt = 0
		cnt += 1

		fHsh = parseFasta(options, i, gGenFighsh[i])
		writeFasta(oDir, i, fHsh)
	err('\n')

def makePLFTabHsh(plfStaHsh, gGenFighsh, plfHsh):
	plfTabHsh = {}
	err("Getting accessory PLF labels...\n\t")
	inc = len(plfStaHsh) / 50.
	cnt = 0
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

def printPLFTab(oFNm, tab):
	f = open(oFNm, 'w')

	for i in tab:
		f.write('\t'.join(i) + '\n')

	f.close()

def printPLFTabHsh(options, plfTabHsh):
	oDir = options.outPref + '.acc.tabs/'
	if not os.path.exists(oDir):
		os.mkdir(oDir)

	err("Writing accessory PLF tabs...\n\t")
	inc = len(plfTabHsh) / 50.
	cnt = 0
	for i in plfTabHsh:
		if cnt >= inc:
			cnt = 0
			err('=')
		cnt += 1

		oFNm = oDir + i + '.tab'
		printPLFTab(oFNm, plfTabHsh[i])
	err('\n')

def printTopPLFs(options, topPLFs):
	oFNm = options.outPref + '.plf.con.lst'

	f = open(oFNm, 'w')
	for i in topPLFs:
		f.write(i + '\n')
	f.close()

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

