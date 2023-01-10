'''
python getClusters.py [kmc dir] [all fasta kmc] [out pref]

Given a directory containing KMC output, a KMC output for all the fasta files concatenated together, and a file output prefix, this script will cluster all the KMC files into clusters of 500, 1000, 2000, and 4000.  From each of these, the script will select one "representative" geonome from each cluster.  These subsampled genomes can then be used as a diverse set.  
'''

from sys import argv,stderr
import os
from random import shuffle
from glob import glob
import numpy as np
from sklearn.cluster import AgglomerativeClustering

# clean up first parameter directory
if argv[1][-1] != '/':
	argv[1] += '/'

# array of clusters to get.  Change it if you want others
clustsToGet = [500, 1000, 2000, 4000]

# output s to stderr
def err(s):
	stderr.write(s)

# this takes a kmc file name and parses it
# returns a hash that maps k-mer to count
def parseKMC(fNm):
	# open file
	f = open(fNm)

	# hash to return
	kHsh = {}
	# for each line in file
	#   split by tab
	#   set k-mer hash
	for i in f:
		i = i.strip('\n').split('\t')
		kHsh[i[0]] = int(i[1])

	f.close()

	return kHsh

# this takes in a file name to a file which contains a list of all
# k-mers.  
# returns a hash that maps a k-mer to an array index
def getKInd(fNm):
	# parse the KMC file
	kInd = parseKMC(fNm)

	# init the index count to 0
	cnt = 0
	# for each element in the k-mer hash
	#   set the value of the k-mer to the count
	#   increment the count by 1
	for i in sorted(kInd):
		kInd[i] = cnt
		cnt += 1

	return kInd

# takes a KMC output file name and k-mer index hash (from getKInd)
# returns an array of k-mer counts
def getKArr(fNm, kInd):
	# parse the KMC file
	kHsh = parseKMC(fNm)
	# init the array of 'nan' values
	arr = [float('nan')] * len(kInd)

	# for each element in the hash
	#   get the k-mer index
	#   set respective array value in array
	for i in kHsh:
		ind = kInd[i]
		arr[ind] = kHsh[i]

	return np.asarray(arr)

# takes in a directory name for a directory which contains all
# the KMC files and a file name which corresponds to the filing 
# containing all k-mers.  
# returns a list containing genome order and a matrix which rows are 
# genomes and columns are k-mer counts
def makeMatrix(dNm, fNm):
	# get the k-mer index hash
	kInd = getKInd(fNm)

	# get list of KMC files
	fLst = glob(dNm + '*.kmrs')

	# init the matrix and genome order 
	mat = []
	gidOrd = []
	# progress bar stuff
	err("Parsing KMC...\n\t")
	inc = len(fLst) / 50.
	cnt = 0
	# for each file in the list
	#   get the gid
	#   get the k-mer count array
	#   append array to matrix
	#   append gid to gid list
	for i in fLst:
		# progress bar stuff
		if cnt >= inc:
			cnt = 0
			err('=')
		cnt += 1

		gid = '.'.join(os.path.basename(i).split('.')[:2])
		arr = getKArr(fNm, kInd)

		mat.append(arr)
		gidOrd.append(gid)
	# progress bar stuff
	err('\n')

	return gidOrd, np.asarray(mat)

# takes in a matrix and number of clusters
# returns cluster predictions for each row in the matrix
def getClust(mat, nClust):
	# init the agglomerative clustering model
	mod = AgglomerativeClustering(n_clusters=nClust, affinity='l1', linkage='average')
	# fit and make predictions
	pred = mod.fit_predict(mat)

	return pred

# given a matrix, genome order list, and output file prefix this gets
# clusters and writes them out to files
def getClusts(mat, gidOrd, oPref):
	# progress stuff
	err("Getting clusters...\n")
	# for each cluster to get
	#   get predictions for the matrix
	#   initialize the cluster hash
	#   for each element in predictions
	#     get the cluster number
	#     get the gid
	#     if the cluster not in cluster hash
	#       add it
	#     add gid to the cluster hash
	#   shuffle the cluster hash for each cluster
	#   open a file and write the first gid in each cluster to file
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

# main driver function
def main():
	# get the GID order and K-mer matrix
	gidOrd, mat = makeMatrix(argv[1], argv[2])
	# get clusters
	getClusts(mat, gidOrd, argv[3])

if __name__ == '__main__':
	main()