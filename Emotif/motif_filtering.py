from __future__ import division
import os
import sys
import argparse
import shutil
import re
from utils import *
import warnings
warnings.filterwarnings("ignore")
import sys
import numpy as np
from scipy.io.arff import loadarff as la
from sklearn import svm
from sklearn import datasets
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.metrics import precision_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import recall_score
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.tree import DecisionTreeClassifier as DTC




def RF_gini_filter(jid,confDict,arff,fileList):
	data = la(arff)
	### get all features
	features=list(data[1])[:-1]
	top = int(confDict['RF_gini_filter']['top'])
	total_ranking_file = jid + "_gini_total_ranking.tsv"
	total_ranking = open(total_ranking_file,"wb")
	Y=np.array(data[0]["Class"])
	X=np.array(map(lambda x:list(x),data[0][features].tolist()))
	a=RFC(criterion='gini')
	a.fit(X,Y)
	count = 0
	hash={}
	for i in a.feature_importances_:
		hash[count]=[i,features[count],i*i]
		count+=1
	rank_list = sorted(hash,key=lambda x:hash[x][2],reverse=True)
	print >>total_ranking,len(rank_list),"feature  \t square of weight value"
	for pos in rank_list:
		print >>total_ranking, hash[pos][1],"\t",hash[pos][2]
	top_list_file = jid + "_gini_top.csv"
	top_list = open(top_list_file,"wb")
	top_motifs = map(lambda x:hash[x][1],rank_list[0:top])
	print >>top_list,",".join(top_motifs)
	fileList.append(top_list_file)
	fileList.append(total_ranking_file)
	return top_motifs
		
def RF_entropy_filter(jid,confDict,arff,fileList):
	data = la(arff)
	### get all features
	features=list(data[1])[:-1]
	top = int(confDict['RF_entropy_filter']['top'])
	total_ranking_file = jid + "_entropy_total_ranking.tsv"
	total_ranking = open(total_ranking_file,"wb")
	Y=np.array(data[0]["Class"])
	X=np.array(map(lambda x:list(x),data[0][features].tolist()))
	a=RFC(criterion='entropy')
	a.fit(X,Y)
	count = 0
	hash={}
	for i in a.feature_importances_:
		hash[count]=[i,features[count],i*i]
		count+=1
	rank_list = sorted(hash,key=lambda x:hash[x][2],reverse=True)
	print >>total_ranking,len(rank_list),"feature  \t square of weight value"
	for pos in rank_list:
		print >>total_ranking, hash[pos][1],"\t",hash[pos][2]
	top_list_file = jid + "_entropy_top.csv"
	top_list = open(top_list_file,"wb")
	top_motifs = map(lambda x:hash[x][1],rank_list[0:top])
	print >>top_list,",".join(top_motifs)
	fileList.append(top_list_file)
	fileList.append(total_ranking_file)
	return top_motifs	





	
		
def prefilter(jid, confDict):

	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	output_folder = jid + '_coverage_filter'
	pos_mhit = confDict['input']['pos_mhit']
	neg_mhit = confDict['input']['neg_mhit']
	maxNegCov = float(confDict['coverage_filter']['maxnegcov'])
	minPosCov = float(confDict['coverage_filter']['minposcov'])
	
	
	os.makedirs(output_folder)
	fileList = []
	
	totalNumPosSeqs = general_utils.findNumSeqs(pos_seq)
	totalNumNegSeqs = general_utils.findNumSeqs(neg_seq)
	print 'numPosSeqs:', totalNumPosSeqs,'numNegSeqs:', totalNumNegSeqs
	#dict between motif names and their seq lists
	posMotifDict = readHitFile(pos_mhit)
	negMotifDict = readHitFile(neg_mhit)
	# for k in negMotifDict.keys():
		# print "Neg Check, ",k,len(negMotifDict[k])
	
	print 'maxNegCov:', maxNegCov,'minPosCov:', minPosCov
	filterPosHitFile = jid + '_pos_filter_hitFile_maxNegCov_' + str(int(maxNegCov*100)) + '.mhit'
	filterNegHitFile = jid + '_neg_filter_hitFile_maxNegCov_' + str(int(maxNegCov*100)) + '.mhit'
	fileList.append(filterPosHitFile)
	fileList.append(filterNegHitFile)
	filtMotifList = filterHits(posMotifDict, negMotifDict, maxNegCov, minPosCov, totalNumPosSeqs, totalNumNegSeqs)
	
	#write hit files
	writeFiltFiles(filtMotifList, filterPosHitFile, filterNegHitFile, posMotifDict, negMotifDict)
	
	#write the pos and neg cov of the filtered motifs
	covInfoFileName = jid + '_cov_info.csv'
	fileList.append(covInfoFileName)
	writeCovFiltFiles(filtMotifList, posMotifDict, negMotifDict, totalNumPosSeqs, totalNumNegSeqs, covInfoFileName)
	
	#posFiltFileName = jid + '_pos_filtered_hits_' + str(minNegCov)
	#fileList.append(posFiltFileName)
	#negFiltFileName = jid + '_neg_filtered_hits_' + str(minNegCov)
	#fileList.append(negFiltFileName)
	
	
	#move files to results folder
	for outFile in fileList:
		shutil.move(outFile, output_folder)
	
	return output_folder,filterPosHitFile,filterNegHitFile

def readHitFile(hitFile):
	"""Read the hit file and make a dict between motifnames and their seqs """
	motifDict = {}
	motifName = ''
	with open(hitFile, 'rb') as handler:
		for line in handler:
			line = line.strip()
			if re.search(r'>', line):
				motifName = line[1:]
				if motifName not in motifDict:
					motifDict[motifName] = []
				continue
			motifDict[motifName].append(line)
	
	# print 'len of motifDict:', len(motifDict)
	return motifDict
	
	
def hit2hit(in_pos_hit,in_neg_hit,out_pos_hit,out_neg_hit,mName):
	posMotifDict = readHitFile(in_pos_hit)
	negMotifDict = readHitFile(in_neg_hit)
	pos_hit = open(out_pos_hit,"wb")	
	neg_hit = open(out_neg_hit,"wb")	
	for i in mName:	
		print >>pos_hit,(">"+i)
		print >>pos_hit,"\n".join(posMotifDict[i])
		if negMotifDict.has_key(i):
			print >>neg_hit,(">"+i)
			print >>neg_hit,"\n".join(negMotifDict[i])
	

def subPwmFile(pwmFileName, selected_motifs_name, allPwmFile):
	pwmFile = open(pwmFileName, 'wb')
	pwmFile.write('MEME version 4.4\nALPHABET= ACGT\nstrands: + -\nBackground letter frequencies (from web form):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n')
	
	
	#read the allPWM file and get the sol PWMs and write them to a file
	flag = 0
	with open(allPwmFile, 'rb') as handler:
		
		for line in handler:
			# print "This is the line in allpwm",line
			line = line.strip()
			if re.search(r'MOTIF', line):
				# print "asd"
				split = line.split()
				motifName  = split[1]
				if motifName in selected_motifs_name:
					# print "In :",motifName
					flag = 1
					pwmFile.write(line + '\n')
					continue
				else:
					flag = 0
					
			if flag == 1:
				pwmFile.write(line + '\n')	
	#close the file
	pwmFile.close()

	