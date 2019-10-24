# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
from __future__ import division
"""
General utils files for general purpose functions

"""
import re
import os
#Biopython
from Bio import SeqIO
import motif_utils

def read_fasta(file):
	from Bio import SeqIO as io
	seq_hash = {}
	f=open(file,"rU")
	record = io.parse(f,"fasta")
	for r in record:
		seq_hash[str(r.id)]=str(r.seq).upper()
	return seq_hash

def getSeqMotifDict(fimoDict):
	"""
	Make a dict between the seq names and list of motifs that occur in it
	Args:
		fimoDict: dict between motif names and the seqs it hits
	Returns:
		a dict between seq names and a list of motif IDs that hit it
	"""
	seqMotifDict = {}
	tmpCount = 0
	for motifId, seqList in fimoDict.items():
		for seqName in seqList:
			if seqName not in seqMotifDict:
				seqMotifDict[seqName] = []
			if motifId not in seqMotifDict[seqName]:
				seqMotifDict[seqName].append(motifId)
				
	
	return seqMotifDict
	
def makeArffFile(jid, fileList, foreFimoDict, backFimoDict, motifDict, foreFileSeqList, backFileSeqList, arffFileName):
	#open the file for writing
	arffFile = open(arffFileName, 'wb')
	#read the fore and backfimo dict and make dict between the seq name and the list of motif hits
	foreSeqMotifDict = getSeqMotifDict(foreFimoDict)
	#
	backSeqMotifDict = getSeqMotifDict(backFimoDict)
	
	#make a dict between the motif names and motif number
	countNameDict = {}
	nameCountDict = {}
	motifCount = 1
	for motifId in motifDict.iterkeys():
		#this is the filtering step applied to remove motifs that have high coverage in the neg data set
		if motifId not in foreFimoDict:
			continue
		countNameDict[motifCount] = motifId
		nameCountDict[motifId] = motifCount
		motifCount += 1
	
	#write the header
	arffFile.write('@RELATION ' + jid + '\n')
	for motifName in motifDict.iterkeys():
		#this is the filtering step applied to remove motifs that have high coverage in the neg data set
		if motifName not in foreFimoDict:
			continue
		arffFile.write('@ATTRIBUTE ' + motifName + ' NUMERIC' + '\n')
	arffFile.write('@ATTRIBUTE ' + 'class ' + '{positive,negative}' + '\n\n')
	arffFile.write('@DATA\n\n')	
	
	##write the first part of the  file (foreground)
	foreLabel = 'positive'
	writePartArff(foreLabel, foreFileSeqList, motifDict, foreSeqMotifDict, arffFile, countNameDict, nameCountDict)
	
	##write the second part of the ARFF file (background)
	backLabel = 'negative'
	writePartArff(backLabel, backFileSeqList, motifDict, backSeqMotifDict, arffFile, countNameDict, nameCountDict)
	
	#close the file
	arffFile.close()
		
def makeMotifDictFromPwmFile(pwmFileName):
	"""
	read a PWM file that we would like to scan and make a motif dict between motif ID and motif object
	Args:
		PWM file in MEME format
	"""
	#print 'PP:', pwmFileName
	#dict between motif IDs and motif objects
	motifDict = {}
	motifName = ''
	motifId = 0
	#dicts between motif IDs and motif names
	with open(pwmFileName, 'rb') as handler:
		for line in handler:
			line = line.strip()
			if not line.strip() or re.search(r'letter', line) or re.search(r'version', line) or re.search(r'ALPHABET',line) or re.search(r'strands', line)\
			or re.search(r'A 0.25000', line):
				continue
			#print 'l:',line
			if re.search(r'MOTIF', line):
				split = line.split()
				print split
				line = '_'.join(split[1:])
				motifName = line
				print "motifName",motifName
				motifId += 1
				motifObj = motif_utils.MyMotif()
				motifObj.Id = motifId
				motifObj.regExp = motifName #notice here it is actually the motif name not a reg. expression
				motifDict[motifName] = motifObj
				continue
			#add the PWM lines to the motif
			motifDict[motifName].pwmLines.append(line)
	
	##check the hash
	#for motifName in motifDict.iterkeys():
		#print motifName
		#motifObj = motifDict[motifName]
		#print '\t',motifObj.regExp
		#for line in motifObj.pwmLines:
			#print '\t', line
			
	return motifDict


def PWMtoPFM(pwmlines,name):
	# only for calculating AUC using DiMO
	out = name + '.pfm'
	f=open(out,"wb")
	line0 = ">motif " + name
	line1 = "A | "
	line2 = "C | "
	line3 = "G | "
	line4 = "T | "
	for line in pwmlines:
		line = line.split()
		line1 = line1 + line[0] + "\t"
		line2 = line2 + line[1] + "\t"
		line3 = line3 + line[2] + "\t"
		line4 = line4 + line[3] + "\t"
	print >>f,line0
	print >>f,line1
	print >>f,line2
	print >>f,line3
	print >>f,line4
	return out
	
def processMotifHitFile_1(inFile):
	""" Process a motif hit/occurrence file and return a universe set (U) and a dict of set of sets
	The inFile should look like:
	>motif_name_1
	sequence name
	sequence name
	>motif_name_2 
	sequence name
	sequence name
	
	You can check the wiki page about set cover problem to get an idea about Universe (U) set and set S of subsets
	
	Args:
		inFile: name of file of motif hits in the format explained above
	
	
	
	Returns:
		motifIdDict: A dict between motif names and motif numbers (IDs)
		idMotifDict: A dict between seq names and seq numbers (IDs)
		seqIdDict: A dict between seq numbers (IDs) and seq names
		idSeqDict: A dict between seq numbers (IDs) and seq names and 
		Uset : universe set of all seq IDs
		Sdict: dict between motif IDs and the sequence IDs the motif hits
	"""
	#a dict between motif names in the file and motif numbers (IDs)
	motifIdDict = {}
	#a dict between motif numbers (IDs) and motif names
	idMotifDict = {}
	#a dict between seq names and seq numbers (IDs)
	seqIdDict = {}
	#a dict between seq numbers and seq names
	idSeqDict = {}
	#all numbers/IDs start from 1 inside the dicts/sets
	motifId = 0
	seqId = 0
	#Universe set which is the set of all the sequences, the U set
	Uset = set()
	#A dict between motif names and the set of seqs that occur in them, the S set
	Sdict = {}
	print 'in process motif hit'
	#start reading the file
	with open(inFile, 'rb') as handler:
		for line in handler:
			line = line.strip()
			#skip empty spaces
			if not line.strip():
				continue
			if re.search(r'>', line):#found a motif
				motifId += 1
				motifName = line[1:]#take off the > character
				if motifName not in motifIdDict:
					motifIdDict[motifName] = motifId
				else:
					print 'There are motif names duplicates. the file should have unique motif names'
					exit()
				idMotifDict[motifId] = motifName
				#initialize the element in the Sdict
				Sdict[motifId] = set()
				continue
			#all the other lines are sequence names, check if the seq name has been checked and inserted or not
			seqName = line
			if seqName not in seqIdDict:
				seqId += 1
				seqIdDict[seqName] = seqId
				idSeqDict[seqId] = seqName
				tmpSeqId = seqId
			else:
				tmpSeqId = seqIdDict[seqName]
			
			if tmpSeqId not in Sdict[motifId]:
				Sdict[motifId].add(tmpSeqId)
			#add the seqs to the U universe dict
			if tmpSeqId not in Uset:
				Uset.add(seqId)
				
	#check the U dict
	#print Uset
	#check the S dict
	#for motifId, seqSet in Sdict.items():
		#print 'motifId:', motifId, 'set:', seqSet
		
	return motifIdDict, idMotifDict, seqIdDict, idSeqDict, Uset, Sdict
	

def processMotifHitFile(foreMotifFile, backMotifFile):
	"""Read the fore and back motif hit files and make motif dict and dict of motifs with their hits
	The input motf hit file should look like:
	>motif_name_1
	sequence name
	sequence name
	>motif_name_2 
	sequence name
	sequence name
	"""
	motifId = 0
	#dict between motif names and motif object
	motifObjDict = {}
	#dict between the motif name and list of seqs it occurs in
	foreHitDict = {}
	#list of sequences in the foreground file
	foreSeqList = []
	motifName = ''
	#start reading the foreground file
	with open(foreMotifFile, 'rb') as handler:
		for line in handler:
			line = line.strip()
			#skip empty spaces
			if not line.strip():
				continue
			if re.search(r'>', line):#found a motif
				motifId += 1
				motifName = line[1:]#take off the > character
				#make a motif object
				motifObj = motif_utils.MyMotif()
				motifObj.Id = motifId
				motifObj.regExp = motifName #notice here it is actually the motif name not a reg. expression
				if motifName not in motifObjDict:
					motifObjDict[motifName] = motifObj
				else:
					print 'motif names should be unique. Exiting'
					exit()
				continue
			#add the seqs to the motif
			if motifName not in foreHitDict:
				foreHitDict[motifName] = []
			foreHitDict[motifName].append(line)
			
			if line not in foreSeqList:
				foreSeqList.append(line)
	
	##check the motif dict
	#for motifName in motifObjDict:
		#print 'mName:', motifName
	#for motifName in foreHitDict:
		#print 'mName:', motifName,'seqs:', foreHitDict[motifName]
	#print 'whole seq list:', foreSeqList
	#process the background
	#dict between the motif name and list of seqs it occurs in the background
	backHitDict = {}
	#list of sequences in the background file
	backSeqList = []
	#check if there is a background file as well and make data for it
	backMotifName = ''
	if backMotifFile != 'none':
		with open(backMotifFile, 'rb') as handler:
			for line in handler:
				line = line.strip()
				#skip empty spaces
				if not line.strip():
					continue
				if re.search(r'>', line):#found a motif
					backMotifName = line[1:]#take off the > character
					continue
				#add the seqs to the motif
				if backMotifName not in backHitDict:
					backHitDict[backMotifName] = []
				backHitDict[backMotifName].append(line)
				
				if line not in backSeqList:
					backSeqList.append(line)
	
	
	#for motifName in backHitDict:
	#	print 'mName:', motifName,'seqs:', backHitDict[motifName]
	#print 'whole seq list:', backSeqList
	
	return motifObjDict, foreHitDict, backHitDict, foreSeqList, backSeqList

def writePWMFromMotifs(motifDict, pwmFileName,prefix):
	"""Go thru the motif dict og BioPython objects and write them to a file in MEME PWM format
	
	Args:
		motif dictionary between motif IDs and  MyMotif objects
	Returns:
		name of PWM motif file
	"""
	pwmFile = open(pwmFileName, 'wb') 
	pwmFile.write('MEME version 4.4\nALPHABET= ACGT\nstrands: + -\nBackground letter frequencies (from web form):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n')
		
	alphaList = ['A', 'C', 'G', 'T']
	#loop thru the motifs
	for motifId in motifDict.iterkeys():
		motifObj = motifDict[motifId].bioMotifObj
		#get motif length
		motifLength = len(motifObj)
		pwmFile.write('\nMOTIF ' + prefix + "_" + str(motifId) + '\n')
		pwmFile.write('letter-probability matrix: alength= 4 w= '+ str(motifLength) + ' nsites= 100 E= 0.001\n') 
		pwm = motifObj.counts.normalize(pseudocounts=0)
		for i in range(motifLength):
			for alpha in alphaList:
				pwmFile.write(' ' + str(pwm[alpha][i]))
			pwmFile.write('\n')
			
	#close the file	
	pwmFile.close()


def writePWMFromMotifsSelected(motifDict, pwmFileName, selectIdList, idMotifDict):
	"""Write the PWMs to a file, only write the IDs chosen in the selectIdList
	
	
	"""
	pwmFile = open(pwmFileName, 'wb') 
	pwmFile.write('MEME version 4.4\nALPHABET= ACGT\nstrands: + -\nBackground letter frequencies (from web form):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n')
		
	alphaList = ['A', 'C', 'G', 'T']
	#loop thru the motifs
	for motifId in selectIdList:
		if motifId not in idMotifDict:
			motifName = motifId
		else:
			motifName = idMotifDict[motifId]
		#print 'motifId:', motifId,motifName
		motifObj = motifDict[motifName].bioMotifObj
		#get motif length
		motifLength = len(motifObj)
		pwmFile.write('\nMOTIF ' + str(motifName) + '\n')
		pwmFile.write('letter-probability matrix:\n')
		pwm = motifObj.counts.normalize(pseudocounts=0)
		for i in range(motifLength):
			for alpha in alphaList:
				pwmFile.write(' ' + str(pwm[alpha][i]))
			pwmFile.write('\n')
				
	#close the file	
	pwmFile.close()
	


def writePWMFromMotifsSelectedScan(motifDict, pwmFileName, selectIdList, idMotifDict, finalSelectMotifList):
	"""Write the PWMs to a file, only write the IDs chosen in the selectIdList
	This is for the motif PWM scan problem where the PWm is already privided as input
	"""
	
	pwmFile = open(pwmFileName, 'wb') 
	pwmFile.write('MEME version 4.4\nALPHABET= ACGT\nstrands: + -\nBackground letter frequencies (from web form):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n')
	
	
	for motifId in selectIdList:
		if motifId not in idMotifDict:
			motifName = motifId
		else:
			motifName = idMotifDict[motifId]
		if motifName not in finalSelectMotifList:
			continue
		motifObj = motifDict[motifName]
		#motifName = motifObj.regExp
		pwmFile.write('\nMOTIF ' + motifName + '\n')
		pwmFile.write('letter-probability matrix:\n')
		for line in motifObj.pwmLines:
			pwmFile.write(line+'\n')
	#close the file	
	pwmFile.close()
	
				
	
def writeMotifSeqFile(motifDict, outFileName):
	"""take a dict between motif IDs and seqs it occurs in and write to a file
	The format is:
	>motif_name
	seq_1
	seq_2
	Args:
		motifDict: dict between motif IDs/names and a list of sequences the motif occurs in
		outFileName" name of file to write to
	"""
	outFile = open(outFileName, 'wb')
	for motifId in sorted(motifDict.iterkeys()):
		seqList = motifDict[motifId]
		outFile.write('>' + str(motifId) + '\n')
		for seqName in seqList:
			outFile.write(seqName + '\n')
		
	outFile.close()
	

def writeMotifSeqFileFromMotifDict(motifDict, outFileName):
	"""Take a motif dict beween motif Ids and MyMotif objects and write the seq of each motif to a file
	use the seq coverage file format
	>motif_name
	seq_1
	seq_2
	Args:
		motifDIct: dict between motif IDs and MyMotif objects
		outFileName" name of file to write to
	"""
	outFile = open(outFileName, 'wb')
	for motifId in motifDict.iterkeys():
		motifObj = motifDict[motifId]
		seqList = motifObj.foreSeqList
		outFile.write('>' + str(motifId) + '\n')
		for seqName in seqList:
			outFile.write(seqName + '\n')
		
	outFile.close()

def findNumSeqs(fastaFile):
	"""Use Biopython to find the number of fasta seqeunces in a fasta file
	
	Args:
		fasta file
	Returns:
		number of fasta sequences
	"""
	handle = open(fastaFile, "rU")
	numSeqs = 0
	for record in SeqIO.parse(handle, "fasta"):
		seqId = str(record.id)
		numSeqs += 1
	handle.close()
	
	return numSeqs


def findSeqList(fastaFile):
	"""
	Read the fasta file and return a list of sequence names in the fasta file 
	"""	
	handle = open(fastaFile, "rU")
	numSeqs = 0
	seqList = []
	for record in SeqIO.parse(handle, "fasta"):
		seqId = str(record.id)
		seqList.append(seqId)
		numSeqs += 1
	handle.close()
	
	return seqList
	
def getSeqList(fastaFile):
	"""Use Biopython to return a list of fasta seq names
	
	Args:
		fasta file
	Returns:
		list of seqNames
	"""
	handle = open(fastaFile, "rU")
	numSeqs = 0
	seqList = []
	for record in SeqIO.parse(handle, "fasta"):
		seqId = str(record.id)
		seqList.append(seqId)
		numSeqs += 1
	handle.close()
	
	return seqList

def makeMotifLogo(motifIdList, motifDict, outDirName, pwmFileName, idMotifDict, op, finalSelectMotifList):
	"""Make motif logos for a set of motifs as specified in the motifIdList
	Add a motif logo ID (logo file name) to each motif object
	Args:
		motifIdList: list of motif IDs to make logos for
		motifDict: a dict between motif ID and a MyMotif object
		idMotifDict: dict from IDs to motif names
		op: operation; is it motif disocvery or a cov 
		finalSelectMotifList: list of motif names , final filtered list
		
	Return:
	
	""" 
	if op == 'disc':	
		writePWMFromMotifsSelected(motifDict, pwmFileName, motifIdList, idMotifDict)
	if op== 'cov':
		writePWMFromMotifsSelectedScan(motifDict, pwmFileName, motifIdList, idMotifDict, finalSelectMotifList)
	#get the path of the module	
	utils_path = os.path.realpath(__file__)
	#remove the name of the module from the end of it. Replace won't work since w ehave DME.py and DME.pyc
	split = utils_path.split('/')
	utils_path = '/'.join(split[:len(split)-1])
	
	inPath = os.path.realpath(__file__)
	split = inPath.split('/')
	inPath = '/'.join(split[:len(split)-2])
	print 'inPath:', inPath
	command =  inPath + '/' + 'meme2images -png ' + pwmFileName + ' ' + outDirName
	#print 'meme2image command:', command
	try:
		os.system(command)
	except:
		print 'meme2image execution failed. Exiting'
		exit()
	
	
def makeMotifLogoFromPwm(pwmFileName, outDirName):
	"""
	Given the pwm file in MEME format make the logos
	"""	
	#get the path of the module	
	utils_path = os.path.realpath(__file__)
	#remove the name of the module from the end of it. Replace won't work since w ehave DME.py and DME.pyc
	split = utils_path.split('/')
	utils_path = '/'.join(split[:len(split)-1])
	
	inPath = os.path.realpath(__file__)
	split = inPath.split('/')
	inPath = '/'.join(split[:len(split)-2])
	command =  inPath + '/' + 'meme2images -png ' + pwmFileName + ' ' + outDirName
	print 'meme2image command:', command
	try:
		os.system(command)
	except:
		print 'meme2image execution failed. Exiting'
		exit()
	


def writeMotifBasicStat(motifDict, foreFimoDict, backFimoDict, foreNumSeqs, backNumSeqs, outFileName, tomtomDict):
	"""
	Write a file with basic info about the motif stats like coverage
	Args:
		motifDict: dict between motif names and motif objects
	"""
	outFile = open(outFileName, 'wb')
	headerList = ['#MotifName','Num_Fg_seqs','Fg_cov','Num_Bg_seqs', 'Bg_cov','FG/BG','similar_motifs']
	headerLine = ','.join(headerList)
	outFile.write(headerLine + '\n')
	#go thru the motif dict
	for motifName in motifDict.iterkeys():
		motifObj = motifDict[motifName]
		if motifName in foreFimoDict:
			foreSeqs = len(foreFimoDict[motifName])
			foreCov = 100*(foreSeqs/foreNumSeqs)
		else:
			foreSeqs = 0
			foreCov = 0
		if motifName in backFimoDict:
			backSeqs = len(backFimoDict[motifName])
			backCov = 100*(backSeqs/backNumSeqs)
			fg_over_bg = foreCov/backCov
		else:
			backSeqs = 0
			backCov = 0
			fg_over_bg = foreCov/1
		
		#find similar motifs to this motif
		simList = []
		simStr = ''
		if motifName in tomtomDict:
			simList = tomtomDict[motifName]
			simStr = '//'.join(simList)
		lineList = [motifName, str(foreSeqs), str(foreCov), str(backSeqs), str(backCov), str(fg_over_bg), simStr]
		lineStr = ','.join(lineList)
		outFile.write(lineStr + '\n')
		
		
		
		
	
	outFile.close()
	

def main(args):
    pass

##
if( __name__ == "__main__" ):
    main(sys.argv)
#--eof--#
