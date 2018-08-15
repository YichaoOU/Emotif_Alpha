# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
from __future__ import division
import sys
import os
import re



def processRequiredoutput(inFileName):
	""" Read the Required output file and get the motif IDs found 
	Args:
		inFileName: name of the Required results file
	"""
	motifIdList = []
	with open(inFileName, 'rb') as handler:
		flag = 0
		for line in handler:
			line = line.strip()
			if re.search(r'>', line):
				lineSplit = line.split(':')
				motifId = lineSplit[2]
				motifId = motifId.strip()
				motifId = motifId[1:]
				motifIdList.append(motifId)
				
	return motifIdList
	


def callRequiredDepth(jid, motifDict, tomtomDict, depth, fileList):
	"""Call the Required algorithm to select motifs
	To run the Required algthm. Usage:rc [-d <coverage_depth>] <input_file> 
	Args:
		motifDict: dict between motif IDs and seqs it occurs in
		tomtomDict: dict between motif Ids and a list of motif Ids that are similar to it
	Returns:
		depthDict: dictionary between the depth and a tuple of (motifIdList, list of sets of seqs per motif, newSeqsAdded)
	"""
	print 'Required sequence coverage with depth:', depth
	#get the path of the module	
	Required_path = os.path.realpath(__file__)
	#remove the name of the module from the end of it. 
	split = Required_path.split('/')
	Required_path = '/'.join(split[:len(split)-1])
	#this list has the motif IDs that have been already written to a file and been selected, used when selecting motifs for depth 2
	chosenMotifId = []
	#dict between a depth value and a tuple of motif Ids, seq sets, and newly added seqs by each motif
	depthDict = {} 
	for i in range(1,depth+1):
		print 'depth:', i
		#prepare the motif seq file
		motifSeqFileName = jid + '_fore_motif_hits_in_seqs_depth_' + str(depth)
		fileList.append(motifSeqFileName)
		outFile = open(motifSeqFileName, 'wb')
		for motifId in sorted(motifDict.iterkeys()):
			if motifId in chosenMotifId:
				continue
			seqList = motifDict[motifId]
			outFile.write('>' + str(motifId) + '\n')
			for seqName in seqList:
				outFile.write(seqName + '\n')
		outFile.close()
		#run the algthm
		Required_path = Required_path + '/rc -d ' + str(1) + ' ' + motifSeqFileName
		RequiredOutFile = jid + '_Required_out_depth_' + str(depth)
		fileList.append(RequiredOutFile)
		cmd = Required_path + ' > ' + RequiredOutFile
		print 'Required cmd:', cmd 
		try:
			os.system(cmd)
		except:
			print 'Required execution failed. Exiting'
		#read and process the Required file and return a list of motif IDs ordered as output by the algthm
		motifIdList = processRequiredoutput(RequiredOutFile)
		print 'motifIdList len:', len(motifIdList)
		#Loop thru the IDs and see how much is the coverage and how many new seqs added , and add these motifs to list of chosen motif IDs
		seqsAddedList = []
		seqSetList = []
		numAddedList = []
		for motifId in motifIdList:
			chosenMotifId.append(motifId)
			seqSet = set(motifDict[motifId])
			seqSetList.append(seqSet)
			numAddedSeqs = 0 
			for seqName in motifDict[motifId]:
				if seqName not in seqsAddedList:
					seqsAddedList.append(seqName)
					numAddedSeqs += 1
			numAddedList.append(numAddedSeqs)
		
		depthDict[i] = (motifIdList, seqSetList, numAddedList)
	
	return depthDict	
	
	
	
		
	
	


def main(args):
    pass

if( __name__ == "__main__" ):
    main(sys.argv)
else:
    pass
