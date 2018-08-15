from __future__ import division
import os
import sys
import argparse
import shutil
import re
from algo import greedy
from algo import ILP
from algo import branch_cut
from algo import required_cover
from utils import *
from utils import Tomtom
from utils import general_utils

def run_rami_welch_clustering(jid, confDict):

	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	output_folder = jid + '_rami_welch_motif_clustering'
	pos_mhit = confDict['input']['pos_mhit']
	neg_mhit = confDict['input']['neg_mhit']
	allPwmFile = confDict['input']['pwm_file']
	fileList = []
	
	#parameters
	#max negative coverage like 20% or 10%
	maxNeg = confDict['prefilter']['maxnegcov']
	maxNeg = maxNeg * 100
	#min positive coverage like 20% or 10%
	minPos = confDict['prefilter']['minposcov']
	minPos = maxNeg * 100
	# minPos = 20
	#max negative coverage per super motif
	maxNegSuper = confDict['rami_welch']['maxnegcov']
	maxNegSuper = maxNegSuper * 100
	# maxNegSuper = 100
	#min positive coverage per super motif
	minPosSuper = confDict['rami_welch']['minposcov']
	minPosSuper = minPosSuper * 100
	# minPosSuper = 0
	#tomtom match e-value
	tomEValue = confDict['rami_welch']['evalue']
	tomStrand = confDict['rami_welch']['strand']
	# tomEValue = 0.0005
	# tomStrand = '+'
	#name of job for pca plot
	jobDesc = 'test'
	
	#find number of seqs in the fasta file
	foreTotalNumSeqs = general_utils.findNumSeqs(pos_seq) 
	backTotalNumSeqs = general_utils.findNumSeqs(neg_seq)
	print 'num fore seqs:', foreTotalNumSeqs, 'num back seqs:', backTotalNumSeqs
	
	#read the hit files and get a dict of motif names and hits
	foreMotifDict = readHitFile(pos_mhit)
	backMotifDict = readHitFile(neg_mhit)
	
	#write to a file the coverage of the selected motifs,
	allMotifInfoFileName = jid + '_all_motifs_info.csv'
	writeMotifInfo(foreMotifDict, backMotifDict, foreTotalNumSeqs, backTotalNumSeqs, allMotifInfoFileName)
	fileList.append(allMotifInfoFileName)
	
	#filter the motifs based on their hits and make a new filtered PWM file
	filtMotifList = filterMotifs(foreMotifDict, backMotifDict, maxNeg, minPos, foreTotalNumSeqs, backTotalNumSeqs)
	print 'len filtmotif list:', len(filtMotifList)
	
	#create a filtered PWM file
	filtPwmFileName = jid + '_filt_PWM.pwm'
	writeFiltPwmFile(inputPwmFile, filtPwmFileName, filtMotifList)
	fileList.append(filtPwmFileName)
	
	#create the motif logos for all the filtered motifs
	allMotifLogoFolderName = jid + '_all_logos'
	general_utils.makeMotifLogoFromPwm(filtPwmFileName, allMotifLogoFolderName)
	fileList.append(allMotifLogoFolderName)
	
	#run tomtom on the list of filtered motifs
	tomtomDict = {}
	#dict that stores parameters for the tomtom tool
	tomParaDict = {}
	tomParaDict['evalue'] = tomEValue
	#tomtom output directory name
	tomtomOutDir = jid + '_Tomtom'
	#we are comparing the motifs against each other, so use the same PWM file
	Tomtom.callTomtom(filtPwmFileName, filtPwmFileName, tomParaDict, tomtomOutDir)
	fileList.append(tomtomOutDir)
	#parse the tomtom output file and return a dict between motif Name and list[] of motifs that match it
	tomtomDict = Tomtom.parse(tomtomOutDir+'/tomtom.txt')
	#2 dimensional dict for e value
	tomEvalDict = parseTomEval(tomtomOutDir+'/tomtom.txt', tomStrand)
	
	#write the set of motif clusters info to a file, sort the motifs on fore cove and the top motifs is considered a seed motif
	setsFileName = jid + '_motif_sets_1.csv'
	setsFile = open(setsFileName, 'wb')
	selectList = []
	#dict between group ID and motif group object
	groupObjDict = {}
	groupId  = 1
	for key, value in sorted(foreMotifDict.iteritems(), key = lambda (k,v):len(v), reverse=True):
		motifName = key
		#check if motif passed the filtering or not
		if motifName not in filtMotifList:
			continue
		#check if motif already was selected as a first motif or not
		if motifName in selectList:
			continue
		seqList = foreMotifDict[motifName]
		#if motifName not in foreMotifDict:
			#continue
		setsFile.write('\n\nfirst motif:' + motifName + '\n')
		matchList = []
		selectList.append(motifName)
		if motifName in tomtomDict:
			for val in tomtomDict[motifName]:
				if val not in foreMotifDict:
					continue
				if val in selectList:
					continue
				setsFile.write(val + '\n')
				matchList.append(val)
				selectList.append(val)
		
		groupObjDict[groupId] = MotifGroup(groupId, motifName, matchList)
			
		groupId += 1
	setsFile.close()
	fileList.append(setsFileName)
	print 'num of motif groups:', len(groupObjDict)	
	
	#file that stores and comapres the stats
	statsFileName = jid + '_group_stats.csv'
	fileList.append(statsFileName)
	statsFile = open(statsFileName, 'wb')
	statsFile.write('Motif,Group_num,clusterSize,reducedCluster,foreground_cov(%),fore_numSeqs,background_cov(%),back_numSeqs\n')
	#create the logo group folders for each motif group
	#folder list for the group logos
	allGroupLogoFolder = jid + '_groups_logos'
	os.makedirs(allGroupLogoFolder)
	allGroupLogoFileList = []
	#create a folder to store the groups hit files and coverage results
	allGroupCovFolder  = jid + '_groups_coverage'
	os.makedirs(allGroupCovFolder)
	allGroupCovList = []
	#file for motif match e-values
	simEvalueFileName = jid + '_similarity_evalues.csv'
	fileList.append(simEvalueFileName)
	simEvalueFile = open(simEvalueFileName, 'wb')
	simEvalueFile.write('seedMotif,similarMotif,E-value\n')
	#file for motif hits of the combined cov selected motifs; super motifs
	foreSuperMotifHitFileName = jid + '_fore_superMotif_hits'
	foreSuperMotifHitFile = open(foreSuperMotifHitFileName, 'wb')
	backSuperMotifHitFileName = jid + '_back_superMotif_hits'
	backSuperMotifHitFile = open(backSuperMotifHitFileName, 'wb')
	fileList.append(foreSuperMotifHitFileName)
	fileList.append(backSuperMotifHitFileName)
	
	
	#loop thru the groups
	for groupId in sorted(groupObjDict.iterkeys()):
		groupObj = groupObjDict[groupId]
		#if len(groupObj.simList) == 0:
			#continue
		print 'group Id:', groupId,'seedMotif:', groupObj.seedMotif,'num matches:', len(groupObj.simList)
		groupLogoFolder = jid + '_group_' + str(groupId)
		#create the folder
		os.makedirs(groupLogoFolder)
		#copy the first motif logo to this folder
		try:
			#MEME logos changes any - to _ , so replace that character, also any '.' in the name is changed to '_' as well
			tmpName = groupObj.seedMotif
			if re.search(r'-', tmpName):
				tmpName = tmpName.replace('-', '_')
			if tmpName.count('.') > 0:
				tmpName = tmpName.replace('.', '_')
			cmd = 'cp ./' + allMotifLogoFolderName + '/logo' + tmpName + '.png ' +  groupLogoFolder
			os.system(cmd)
		except:
			print 'copy logos failed'
			exit()
		#copy the similar motif logos
		for simMotif in groupObj.simList:
			tmpName = simMotif
			if re.search(r'-', tmpName):
				tmpName = tmpName.replace('-', '_')
			if tmpName.count('.') > 0:
				tmpName = tmpName.replace('.', '_')
			try:
				cmd = 'cp ' + allMotifLogoFolderName + '/logo' + tmpName + '.png ' +  groupLogoFolder
				os.system(cmd)
			except:
				print 'copy similar logos failed'
				exit()
		allGroupLogoFileList.append(groupLogoFolder)
		
		###write the match e-values
		simEvalueFile.write(groupObj.seedMotif + '\n')
		for simMotif in groupObj.simList:
			eValue = tomEvalDict[groupObj.seedMotif][simMotif]
			simEvalueFile.write(',' + simMotif + ',' + str(eValue) + '\n')
		
		####################################################################
		#make a hits file per group to run motif selection on it
		#create the folder
		####################################################################
		groupCovFolder = jid + '_group_' + str(groupId) + '_cov'
		os.makedirs(groupCovFolder)
		covFolderList = []
		foreGroupHitFileName = jid + '_group_' + str(groupId) + '_fore_hits'
		backGroupHitFileName = jid + '_group_' + str(groupId) + '_back_hits'
		writeGroupHits(foreGroupHitFileName, backGroupHitFileName, groupObj.seedMotif, groupObj, foreMotifDict, backMotifDict)
		covFolderList.append(foreGroupHitFileName)
		covFolderList.append(backGroupHitFileName)
		
		
		
		#********************* motif selection on motif similarity groups ********************# 
		#call the greedy algthm on the group hit files
		jobID = jid + '_Greedy_group_' + str(groupId)
		groupConfFileName = jid + '_confFile_group_' + str(groupId) + '.conf'
		
		
		
		callGreedySelect_per_group(groupObj.seedMotif, groupObj.simList, foreGroupHitFileName, backGroupHitFileName, inputPwmFile, groupConfFileName, motifSelectDefaultConfFile, foreFastaFile, backFastaFile, jobID)
		
		
		covFolderList.append(groupConfFileName)
		covFolderList.append(jobID)
		#read the depth results file to get how many motifs produced and their coverage
		covForeCumCov, covBackCumCov, greedyMotifList = readGreedyResults(jobID, statsFile, foreGroupHitFileName, groupId) 
		print 'group ID:', groupId, 'greedy list:', greedyMotifList
		#read the greedyMotifList and make a new hit file
		#create a new hit file based on these requirements
		foreGroupSelectHitFileName = jid + '_pos.mhit'
		foreGroupSelectHitFile = open(foreGroupSelectHitFileName, 'wb')
		covFolderList.append(foreGroupSelectHitFileName)
		backGroupSelectHitFileName = jid + '_neg.mhit'
		backGroupSelectHitFile = open(backGroupSelectHitFileName, 'wb')
		covFolderList.append(backGroupSelectHitFileName)
		makeHitFile(foreGroupSelectHitFile, backGroupSelectHitFile, foreMotifDict, backMotifDict, greedyMotifList)
		foreGroupSelectHitFile.close()
		backGroupSelectHitFile.close()
		#read the group cov selected hit file and add it to the super motif all hit file. This will be the hit file of the super motifs
		#check how much was the fore and back cum coverage and add if pass the check
		if covForeCumCov >= minPosSuper and covBackCumCov <= maxNegSuper:
			writeHitsToSuperMotifFile(foreGroupSelectHitFileName, foreSuperMotifHitFile)
			writeHitsToSuperMotifFile(backGroupSelectHitFileName, backSuperMotifHitFile)
		
		#move the group cov files and folders
		for outFile in covFolderList:
			shutil.move(outFile, groupCovFolder)
		allGroupCovList.append(groupCovFolder)
		##break outer loop
		#break
	
	#close the sim e value file
	simEvalueFile.close()
	statsFile.close()
	#move the group logo folders
	for outFile in allGroupLogoFileList:
		shutil.move(outFile, allGroupLogoFolder)
	fileList.append(allGroupLogoFolder)
	#move the cov folders
	for outFile in allGroupCovList:
		shutil.move(outFile, allGroupCovFolder)
	fileList.append(allGroupCovFolder)
	
	
	foreSuperMotifHitFile.close()
	backSuperMotifHitFile.close()

	
	
	#move files to results folder
	for outFile in fileList:
		shutil.move(outFile, resultsDirName)
	exit()
	
	return output_folder

class MotifGroup:
	"""general word (site) class to store info for a specific word
	
	"""
	def __init__(self, groupId, seedMotif, simList):
		self.groupId = groupId     #the ID of this group which is based on the ranking of the seed motif
		self.seedMotif = seedMotif  #the seed motif in the group which has the highest fore coverage
		self.simList = simList    #the list of similar motifs



def parseTomEval(inFile, strand):
	"""parse the tomtotm.txt file and return a dict of e-values between motifs """
	#dict between query and target Ids
	matchDict = {}
	#read the tomtom.txt file and see if any known are found
	with open(inFile, 'rb') as handler:
		for line in handler:
			line = line.strip()
			if re.search(r'Query', line):
				continue
			#skip if the target is reverse compliment 
			lineSplit = line.split('\t')
			#when making the dict, the targets are usually motif IDs as digits so check them, otherwise they are strings/names
			if lineSplit[0].isdigit():
				query = int(lineSplit[0])
			else:
				query = lineSplit[0]
			if lineSplit[1].isdigit():
				target = int(lineSplit[1])
			else:
				target = lineSplit[1]
			
			eVal = float(lineSplit[4])
			overlap = int(lineSplit[6])
			strand = lineSplit[9]
			#if strand != '+':
				#continue
			if query not in matchDict:
				matchDict[query] = {}
			if target not in matchDict[query] and target != query:
				matchDict[query][target]=eVal
	
	
	#for query in matchDict:
		#print 'query:', query
		#for target in matchDict[query]:
			#print '\ttarget:', target,'eval:', matchDict[query][target]
	
	return matchDict

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
	
	print 'len of motifDict:', len(motifDict)
	return motifDict


def writeMotifInfo(foreMotifDict, backMotifDict, totalNumPosSeqs, totalNumNegSeqs, outFileName):
	outFile = open(outFileName, 'wb')
	outFile.write('Motif,PosCov,PosNumSeqs,NegCov,NegNumSeqs\n')
	for motifName in foreMotifDict:
		posNumSeqs = len(foreMotifDict[motifName])
		posCov = round( 100*(posNumSeqs/totalNumPosSeqs), 1)
		if motifName in backMotifDict:
			negNumSeqs = len(backMotifDict[motifName])
			negCov = round( 100*(negNumSeqs/totalNumNegSeqs), 1)
		else:
			negNumSeqs = 0
			negCov= 0
		lineList = [motifName, str(posCov), str(posNumSeqs), str(negCov), str(negNumSeqs)]
		lineStr = ','.join(lineList)
		outFile.write(lineStr + '\n')
		
	outFile.close()


def filterMotifs(foreMotifDict, backMotifDict, maxNeg, minPos, foreTotalNumSeqs, backTotalNumSeqs):
	"""filter the motifs based on fore and back coverage """
	filtMotifList = []
	for motifName in foreMotifDict:
		posNumSeqs = len(foreMotifDict[motifName])
		posCov = round( 100*(posNumSeqs/foreTotalNumSeqs), 1)
		if motifName in backMotifDict:
			negNumSeqs = len(backMotifDict[motifName])
			negCov = round( 100*(negNumSeqs/backTotalNumSeqs), 1)
		else:
			negCov = 0
		#check the cov
		if posCov >= minPos and negCov < maxNeg:
			filtMotifList.append(motifName)
	
	print 'maxNeg:', maxNeg, 'minPos:', minPos, 'num of filtered motifs:', len(filtMotifList)		
	return filtMotifList 
		

def writeFiltPwmFile(inputPwmFile, filtPwmFileName, filtMotifList):
	filtPwmFile = open(filtPwmFileName, 'wb')
	filtPwmFile.write('MEME version 4.4\nALPHABET= ACGT\nstrands: + -\nBackground letter frequencies (from web form):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n')
	
	flag = 0
	count = 0
	with open(inputPwmFile, 'rb') as handler:
		for line in handler:
			line = line.strip()
			if re.search(r'MOTIF', line):
				split = line.split()
				motifName = split[1]
				if motifName not in filtMotifList:
					flag = 0
					continue
				else:
					filtPwmFile.write(line + '\n')
					flag = 1
					count += 1
					continue
			
			if flag == 1:
				filtPwmFile.write(line + '\n') 
			
	print 'count:', count
	filtPwmFile.close()


def writeGroupHits(foreGroupHitFileName, backGroupHitFileName, seedMotif, groupObj, foreMotifDict, backMotifDict):
	"""write the motif hits per group """
	foreGroupHitFile = open(foreGroupHitFileName, 'wb')
	backGroupHitFile = open(backGroupHitFileName, 'wb')
	
	#foreSeqList = []
	#backSeqList = []
	
	#write the hits for the seed motif
	foreGroupHitFile.write('>' + seedMotif + '\n')
	if seedMotif in backMotifDict:
		backGroupHitFile.write('>' + seedMotif + '\n')
	for seqName in foreMotifDict[seedMotif]:
		foreGroupHitFile.write(seqName + '\n')
		#foreSeqList.append(seqName)
	if seedMotif in backMotifDict:
		for seqName in backMotifDict[seedMotif]: 
			backGroupHitFile.write(seqName + '\n')
			#backSeqList.append(seqName)
	
	#write the hits for the similar motifs
	for simMotif in groupObj.simList:
		foreGroupHitFile.write('>' + simMotif + '\n')
		for seqName in foreMotifDict[simMotif]:
			#if seqName not in foreSeqList:
			foreGroupHitFile.write(seqName + '\n')
			#foreSeqList.append(seqName)
		if simMotif in backMotifDict:
			backGroupHitFile.write('>' + simMotif + '\n')
			for seqName in backMotifDict[simMotif]:
				#if seqName not in backSeqList:
				backGroupHitFile.write(seqName + '\n')
				#backSeqList.append(seqName) 
	
	foreGroupHitFile.close()
	backGroupHitFile.close()



def writeGreedySelectConfFile(foreHitFileName, backHitFileName, pwmFileName, groupConfFileName, defaultConfFileName, foreFastaFile, backFastaFile):
	"""
	Read the conf file and find which elements to change
	"""
	outFile = open(groupConfFileName, 'wb')
	with open(defaultConfFileName, 'rb') as handler:
		for line in handler:
			line = line.strip()
			if re.search(r'fore_testing_file=', line):
				outFile.write('fore_testing_file=' + foreFastaFile + '\n')
				continue
			if re.search(r'back_testing_file=', line):
				outFile.write('back_testing_file=' + backFastaFile + '\n')
				continue
			if re.search(r'pwm_file=', line):
				outFile.write('pwm_file=' + pwmFileName + '\n')
				continue
			if re.search(r'foremotifhits_file=', line):
				outFile.write('foremotifhits_file=' + foreHitFileName + '\n')
				continue
			if re.search(r'backmotifhits_file=', line):
				outFile.write('backmotifhits_file=' + backHitFileName + '\n')
				continue
			outFile.write(line + '\n')
	
	outFile.close()


def callGreedySelect_per_group(seedMotif, simList, foreGroupHitFileName, backGroupHitFileName, pwmFileName, groupConfFileName, defaultConfFileName, foreFastaFile, backFastaFile, jobID):
	""" 
	seedMotif: the seed motif of the group
	simList: list of motifs similar to the seed motif
	"""
	#make the conf file
	writeGreedySelectConfFile(foreGroupHitFileName, backGroupHitFileName, pwmFileName, groupConfFileName, defaultConfFileName, foreFastaFile, backFastaFile)
	flag = 1
	cmd = 'python /home/rami/Documents/seq_cov/main.py -jid ' + jobID + ' -confFile ' +  groupConfFileName
	print 'CMD:', cmd
	os.system(cmd)
	


def readGreedyResults(jobID, statsFile, groupForeFileName, groupCount):
	""" read the depth results file to get how many motifs produced and their coverage """
	
	#find how many motifs in the combined list
	totalNumMotifs = 0
	with open(groupForeFileName, 'rb') as handler:
		for line in handler:
			if re.search(r'>', line):
				totalNumMotifs += 1
		
	motifNameList = []
	cumForeNumSeqs = 0
	cumBackNumSeqs = 0
	for fileName in os.listdir(jobID):
		if re.search(r'_depth_results', fileName):
			covLineCount = 0
			flag = 1
			with open(jobID + '/' + fileName, 'rb') as handler:
				for line in handler:
					if not line.strip() or re.search(r'foreground', line) or re.search(r'Depth', line):
						continue
					line = line.strip()
					covLineCount += 1
					split = line.split(',')
					motifNameList.append(split[0])
					cumForeNumSeqs += int(split[3])
					cumBackNumSeqs += int(split[7])
					#to get the seed motif name
					if flag == 1:
						mName = split[0]
						flag=0
			#process the last line
			covForeCumCov = float(split[4])
			covBackCumCov = float(split[8])
			#print 'cov numMotifs:', covLineCount, 'cumCov:', covForeCumCov, 'backCov:', covBackCumCov
			#print info to stats file
			statsFile.write(mName + ',' + str(groupCount) + ',' + str(totalNumMotifs) + ',' +  str(covLineCount) + ',' +  str(covForeCumCov) + ',' + str(cumForeNumSeqs) + ',' 
			+ str(covBackCumCov) + ',' + str(cumBackNumSeqs) +'\n')
			
	return covForeCumCov, covBackCumCov, motifNameList



def makeHitFile(foreFinalHitFile, backFinalHitFile, foreMotifDict, backMotifDict, greedyMotifList):
	""" create a new hit file based on the reduced combined list of motifs """
	count = 0
	foreSeqList = []
	backSeqList = []
	#the motifs already sorted on fore cov 
	for motifName in greedyMotifList:
		if count == 0:
			foreFinalHitFile.write('>' + motifName + '\n')
			if motifName in backMotifDict:
				backFinalHitFile.write('>' +  motifName + '\n')
			count += 1
		
		for seqName in foreMotifDict[motifName]:
			if seqName not in foreSeqList:
				foreFinalHitFile.write(seqName + '\n')
				foreSeqList.append(seqName)
		
		if motifName in backMotifDict:
			for seqName in backMotifDict[motifName]:
				if seqName not in backSeqList:
					backFinalHitFile.write(seqName + '\n')
					backSeqList.append(seqName)
			
def writeHitsToSuperMotifFile(inHitFile, foreSuperMotifHitFile):
	with open(inHitFile, 'rb') as handler:
		for line in handler:
			foreSuperMotifHitFile.write(line)

def callGreedySelect_super(foreSuperMotifHitFileName, backSuperMotifHitFileName, inputPwmFile, superConfFileName, defaultConfFileName, foreFastaFile, backFastaFile, jobID):
	#make the conf file
	writeGreedySelectConfFile(foreSuperMotifHitFileName, backSuperMotifHitFileName, inputPwmFile, superConfFileName, defaultConfFileName, foreFastaFile, backFastaFile)
	cmd = 'python /home/rami/Documents/seq_cov/main.py -jid ' + jobID + ' -confFile ' +  superConfFileName
	print 'CMD:', cmd
	os.system(cmd)




















		