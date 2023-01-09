'''
python getClusters.py [kmc dir] [all fasta kmc] [out pref]
'''

from sys import argv,stderr
import os
from random import shuffle
from glob import glob
import numpy as np
from sklearn.cluster import AgglomerativeClustering

if argv[1][-1] != '/':
	argv[1] += '/'

clustsToGet = [500, 1000, 2000, 4000]

def err(s):
	stderr.write(s)

def parseKMC(fNm):
	f = open(fNm)

	kHsh = {}
	for i in f:
		i = i.strip('\n').split('\t')
		kHsh[i[0]] = int(i[1])

	f.close()

	return kHsh

def getKInd(fNm):
	kInd = parseKMC(fNm)

	cnt = 0
	for i in sorted(kInd):
		kInd[i] = cnt
		cnt += 1

	return kInd

def getKArr(fNm, kInd):
	kHsh = parseKMC(fNm)
	arr = [float('nan')] * len(kInd)
	for i in kHsh:
		ind = kInd[i]
		arr[ind] = kHsh[i]

	return np.asarray(arr)

def makeMatrix(dNm, fNm):
	kInd = getKInd(fNm)

	fLst = glob(dNm + '*.kmrs')

	mat = []
	gidOrd = []
	err("Parsing KMC...\n\t")
	inc = len(fLst) / 50.
	cnt = 0
	for i in fLst:
		if cnt >= inc:
			cnt = 0
			err('=')
		cnt += 1

		gid = '.'.join(os.path.basename(i).split('.')[:2])
		arr = getKArr(fNm, kInd)

		mat.append(arr)
		gidOrd.append(gid)
	err('\n')

	return gidOrd, np.asarray(mat)

def getClust(mat, nClust):
	mod = AgglomerativeClustering(n_clusters=nClust, affinity='l1', linkage='average')
	pred = mod.fit_predict(mat)

	return pred

def getClusts(mat, gidOrd, oPref):
	err("Getting clusters...\n")
	for i in clustsToGet:
		err('\t' + str(i) + '\n')
		pred = getClust(mat, i)

		cHsh = {}
		for j in range(0,len(pred)):
			cNum = pred[j]
			gid = gidOrd[j]

			if cNum not in cHsh:
				cHsh[cNum] = []
			cHsh[cNum].append(gid)

		for j in cHsh:
			shuffle(cHsh[j])

		f = open(oPref + '.' + str(i) + '.clusts', 'w')
		for j in cHsh:
			f.write('\t'.join(cHsh[j]) + '\n')
		f.close()

def main():
	gidOrd, mat = makeMatrix(argv[1], argv[2])
	getClusts(mat, gidOrd, argv[3])

if __name__ == '__main__':
	main()