'''
python predict.py [-opt val] | [--opt val]

-f | --feature_tab : feature tabular file from an annotated genome on BVBRC
-p | --cons_plfs : list of conserved PLFs used to train the model
-m | --models : directory containing all the models which will be predicted on
-t | --temp_dir : Temporary directory used to hold fasta files, KMC output, etc.  This directory may get cleared when running the script.  Defaults to "temp/"
'''

from sys import stderr
import sys
from optparse import OptionParser
import os
from ast import literal_eval
import shutil
from glob import glob
import xgboost as xgb

# outputs s to stderr
def err(s):
	stderr.write(s)

# given a directory name and a remove flag creates a directory 
# if it doesn't exist and clears it if rm is set to True
def makeDir(d, rm = True):
	# flag to check if successful
	flag = True
	# if directory isn't a directory
	#   if path exists it's not a directory
	#     output error
	#     and set flag to False
	#   else make the directory
	# else it is a directory
	#   if flag is true, remove the directory and recreate it
	if not os.path.isdir(d):
		if os.path.exists(d):
			print d, 'is not a directory'
			flag = False
		else:
			os.mkdir(d)
	else:
		if rm:
			shutil.rmtree(d)
			os.mkdir(d)

	return flag

# given a directory checks if it exists
def checkDir(d):
	if os.path.isdir(d):
		return True
	else:
		print 'directory does not exist:', d
		return False

# given a directory name, cleans it so it ends with a '/'
def cleanDir(d):
	if d == '':
		return d

	if d[-1] != '/':
		d += '/'

	return d

# grab options for the script and returns them
def getOptions():
	parser = OptionParser()

	parser.add_option('-f', '--feature_tab', help="Feature tabular file from the annotated genome to predict on", metavar="FILE", default='', dest="featFile")
	parser.add_option('-p', '--cons_plfs', help="List of conserved PLFs that the model is based off of.", metavar='FILE', default='', dest='conPLFFile')
	parser.add_option('-m', '--models_dir', help="Directory name containing all the models to preedict with.  Each model would be a directory named PLF_XXX_XXXXXXXX.tab", metavar="DIR", default='', dest='modelsDir')
	parser.add_option('-t', '--temp_dir', help="Directory name to be used as a temporary directory for temporary files and other things", metavar="DIR", default='temp/', dest='tempDir')

	options,args = parser.parse_args()

	options.modelsDir = cleanDir(options.modelsDir)
	options.tempDir = cleanDir(options.tempDir)

	options.fastaFile = os.path.basename(options.featFile).replace('.txt','.fasta', 1)

	makeDir(options.tempDir, True)

	return options, parser

# given options, gets the list of conserved genes
def getConPLFs(options):
	# open the conserved genes file
	f = open(options.conPLFFile)

	# init the hash to hold conserved gene names
	cPLFHsh = {}
	# for each line, add it to the hash
	for i in f:
		i = i.strip('\n')
		cPLFHsh[i] = ''

	f.close

	return cPLFHsh

# given file stream parses header
# returns hash that maps header name to index
def parseHeader(f):
	# split line
	ln = f.readline().strip('\n').split('\t')
	# init hash
	headHsh = {}
	# for each element in line, set hash
	for i in range(0,len(ln)):
		headHsh[ln[i]] = i

	return headHsh

# given options and conserved genes parses the feature file
# for the genome and extracts the conserved PLFs and nucleotide
# sequences for the given genome.
# returns a hash that maps conserved figs to their sequences
def parseFeatures(options, cPLFHsh):
	# open the feature file
	f = open(options.featFile)

	# get the header setup
	headHsh = parseHeader(f)
	# init the hash
	gFigHsh = {}
	# for each line in the file
	#   split the line
	#   get the fig
	#   get the PLF
	#   get the nucl seq
	#   if the plf is one of the conserved ones
	#     add to hash
	for i in f:
		i = i.strip('\n').split('\t')
		figID = i[headHsh['feature_id']]
		plfam = i[headHsh['plfam']]
		dnaSe = i[headHsh['nucleotide_sequence']]

		if plfam in cPLFHsh:
			gFigHsh[figID] = dnaSe

	f.close()

	return gFigHsh

# given options, creates a fasta of conserved genes (which was 
# trained on) and writes it to a file
def makeConsFastaFile(options):
	# get the cosnerved genes
	cPLFHsh = getConPLFs(options)
	# parse the feature tabular file for the cosnerved genes and
	# sequences
	gFigHsh = parseFeatures(options, cPLFHsh)

	# open file
	f = open(options.tempDir + os.path.basename(options.fastaFile), 'w')

	# for each element in the hash, write it in fasta file
	for i in gFigHsh:
		f.write('>' + i + '\n')
		f.write(gFigHsh[i] + '\n')

	f.close()

# given options, finds the k-mer size used to train the models
def getK(options):
	# get the first model trained
	dNm = glob(options.modelsDir + '*.tab/')[0]

	# open params model
	# reads line
	# closes file
	f = open(dNm + 'model.params')
	ln = f.readline()
	f.close()

	# the params are a hash in string form, so we do a literal_eval
	hsh = literal_eval(ln)

	# return the k-mer size
	return hsh['kmerSize']

# given options and k-mer size, runs KMC
def runKMC(options, k):
	# create command to run KMC
	cmdArr = ['kmc.sh', str(k), options.tempDir + os.path.basename(options.fastaFile), options.tempDir + os.path.basename(options.fastaFile), options.tempDir,  '1', '> /dev/null']
	cmd = ' '.join(cmdArr)

	# run KMC
	os.system(cmd)

	# create command to delete intermediate KMC output
	cmdArr = ['rm', options.tempDir + os.path.basename(options.fastaFile) + '.kmc_*']
	cmd = ' '.join(cmdArr)

	# remove intermediate KMC output
	os.system(cmd)

# given a file name, parses the KMC file
# returns hash that maps k-mer to count
def parseKMC(fNm):
	# open the file
	f = open(fNm)

	# init the k-mer hash
	kHsh = {}
	# for each line in the file
	#   split the line
	#   add k-mer and count to hash
	for i in f:
		i = i.strip('\n').split('\t')
		kHsh[i[0]] = int(i[1])

	f.close()

	return kHsh

# given options, gets the order of features for train/prediction
# matrix.  
# returns hash that maps feature to index
def getAttrOrder(options):
	# get the first model trained
	dNm = glob(options.modelsDir + '*.tab/')[0]

	# open attribute order file
	f = open(dNm + 'model.attrOrder')

	# init the hash
	attrOrder = {}
	# for each line, split then add to hash
	for i in f:
		i = i.strip('\n').split('\t')
		attrOrder[i[0]] = int(i[1])

	f.close()

	return attrOrder

# given options and k-mer size, creates a matrix to predict with
# returns a DMatrix to predict XGB models
def makeMatrix(options, k):
	# parse the KMC file
	kHsh = parseKMC(options.tempDir + os.path.basename(options.fastaFile) + '.' + str(k) + '.kmrs')
	# get attribute order
	attrOrder = getAttrOrder(options)

	# init the array
	arr = [0] * len(attrOrder)

	# for each k-mer in the hash
	#   if the k-mer is not in the features list
	#     continue
	#   get the index number for the feature
	#   set array appropriately
	for i in kHsh:
		if i not in attrOrder:
			continue
		ind = attrOrder[i]
		arr[ind] = kHsh[i]

	# convert the array to a DMatrix
	dMat = xgb.DMatrix([arr])

	return dMat

# given a DMatrix and directory name for a model predicts the
# presence/absence of the gene from the model for the genome
# returns a set of predictions.
def predPLF(dMat, dNm):
	# get list of folds to predict with
	mLst = glob(dNm + 'all/model*pkl')

	# array to hold predictions
	preds = []
	# for each fold
	#   load the model
	#   set threads to 1
	#   make prediction
	#   append to array
	for i in mLst:
		mod = xgb.Booster(model_file=i)
		mod.set_param('njobs', 1)
		pred = mod.predict(dMat)
		preds.append(pred[0])

	return preds

# given an array of predictions converts to prediction string
#   N = strong absence
#   N* = weak absence
#   Y = strong presence
#   Y* = weak presence
# returns prediction string
def convPredArr(pArr):
	s = sum(pArr) / float(len(pArr))

	sLab = ''
	if s < 0.33:
		sLab = 'N'
	elif s < 0.5:
		sLab = 'N*'
	elif s < 0.66:
		sLab = 'Y*'
	elif s > 0.66:
		sLab = 'Y'

	return sLab

# given options and a DMatrix, predicts on all models
# returns a table of predictions
def predAll(options, dMat):
	# get list of models
	dLst = glob(options.modelsDir + '*.tab/')

	# init predictions
	preds = []
	# progress bar stuff
	err("Making PLF predictions...\n\t")
	inc = len(dLst) / 50
	cnt = 0
	# for each model
	#   get PLF ID
	#   make predictions for model
	#   get the string label for the fold
	#   append to table of predictions
	for i in dLst:
		if cnt >= inc:
			err('=')
			cnt = 0
		cnt += 1

		plf = os.path.basename(i[:-1]).replace('.tab', '', 1)
		predArr = predPLF(dMat, i)
		sLab = convPredArr(predArr)
		preds.append([plf, sLab, str(sum(predArr))])
	err('\n')

	return preds

# given table of predictions, prints them
def printPreds(preds):
	# create and output header
	arr = ['PLFam', 'Prediction']
	print '\t'.join(arr)
	# for each prediction, print it
	for i in preds:
		print '\t'.join(i)

# main driver program
# get options
# get k-mer size
# make conserved fasta file
# run KMC on said fasta file
# create DMatrix
# make predictions
# print predictions
def main():
	options, parser = getOptions()
	k = getK(options)
	makeConsFastaFile(options)
	runKMC(options, k)
	dMat = makeMatrix(options, k)
	preds = predAll(options, dMat)
	printPreds(preds)

if __name__ == '__main__':
	main()
