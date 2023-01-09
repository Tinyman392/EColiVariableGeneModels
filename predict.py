'''
python predict.py
'''

from sys import stderr
import sys
from optparse import OptionParser
import os
from ast import literal_eval
import shutil
from glob import glob
import xgboost as xgb

def err(s):
	stderr.write(s)

def makeDir(d, rm = True):
	flag = True
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

def checkDir(d):
	if os.path.isdir(d):
		return True
	else:
		print 'directory does not exist:', d
		return False

def cleanDir(d):
	if d == '':
		return d

	if d[-1] != '/':
		d += '/'

	return d

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

def getConPLFs(options):
	f = open(options.conPLFFile)

	cPLFHsh = {}
	for i in f:
		i = i.strip('\n')
		cPLFHsh[i] = ''

	f.close

	return cPLFHsh

def parseHeader(f):
	ln = f.readline().strip('\n').split('\t')
	headHsh = {}
	for i in range(0,len(ln)):
		headHsh[ln[i]] = i

	return headHsh

def parseFeatures(options, cPLFHsh):
	f = open(options.featFile)

	headHsh = parseHeader(f)
	gFigHsh = {}
	for i in f:
		i = i.strip('\n').split('\t')
		figID = i[headHsh['feature_id']]
		plfam = i[headHsh['plfam']]
		dnaSe = i[headHsh['nucleotide_sequence']]

		if plfam in cPLFHsh:
			gFigHsh[figID] = dnaSe

	f.close()

	return gFigHsh

def makeConsFastaFile(options):
	cPLFHsh = getConPLFs(options)
	gFigHsh = parseFeatures(options, cPLFHsh)

	f = open(options.tempDir + os.path.basename(options.fastaFile), 'w')

	for i in gFigHsh:
		f.write('>' + i + '\n')
		f.write(gFigHsh[i] + '\n')

	f.close()

def getK(options):
	dNm = glob(options.modelsDir + '*.tab/')[0]

	f = open(dNm + 'model.params')
	ln = f.readline()
	f.close()

	hsh = literal_eval(ln)

	return hsh['kmerSize']

def runKMC(options, k):
	cmdArr = ['kmc.sh', str(k), options.tempDir + os.path.basename(options.fastaFile), options.tempDir + os.path.basename(options.fastaFile), options.tempDir,  '1', '> /dev/null']
	cmd = ' '.join(cmdArr)

	os.system(cmd)

	cmdArr = ['rm', options.tempDir + os.path.basename(options.fastaFile) + '.kmc_*']
	cmd = ' '.join(cmdArr)

	os.system(cmd)

def parseKMC(fNm):
	f = open(fNm)

	kHsh = {}
	for i in f:
		i = i.strip('\n').split('\t')
		kHsh[i[0]] = int(i[1])

	f.close()

	return kHsh

def getAttrOrder(options):
	dNm = glob(options.modelsDir + '*.tab/')[0]

	f = open(dNm + 'model.attrOrder')

	attrOrder = {}
	for i in f:
		i = i.strip('\n').split('\t')
		attrOrder[i[0]] = int(i[1])

	f.close()

	return attrOrder

def makeMatrix(options, k):
	kHsh = parseKMC(options.tempDir + os.path.basename(options.fastaFile) + '.' + str(k) + '.kmrs')
	attrOrder = getAttrOrder(options)

	arr = [0] * len(attrOrder)

	for i in kHsh:
		if i not in attrOrder:
			continue
		ind = attrOrder[i]
		arr[ind] = kHsh[i]

	dMat = xgb.DMatrix([arr])

	return dMat

def predPLF(dMat, dNm):
	mLst = glob(dNm + 'all/model*pkl')

	preds = []
	for i in mLst:
		mod = xgb.Booster(model_file=i)
		mod.set_param('njobs', 1)
		pred = mod.predict(dMat)
		preds.append(pred[0])

	return preds

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

def predAll(options, dMat):
	dLst = glob(options.modelsDir + '*.tab/')

	preds = []
	err("Making PLF predictions...\n\t")
	inc = len(dLst) / 50
	cnt = 0
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

def printPreds(preds):
	arr = ['PLFam', 'Prediction']
	print '\t'.join(arr)
	for i in preds:
		print '\t'.join(i)

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
