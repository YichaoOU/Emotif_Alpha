# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
from __future__ import division
import sys
import os
import re

def processBranchoutput(inFileName):
	""" Read the Branch and cut output file and get the motif IDs found 
	Args:
		inFileName: name of the Branch and cut results file
	"""
	motifIdList = []
	with open(inFileName, 'rb') as handler:
		flag = 0
		for line in handler:
			line = line.strip()
			if re.search(r'Found one of size', line):
				flag = 1
				continue
			if flag == 1:
				if re.search(r'>', line):
					motifId = line[1:]
					motifIdList.append(motifId)
	return motifIdList
	


def callBranchDepth(jid, motifDict, tomtomDict, depth, fileList):
	"""Call the Branch and cut algorithm to select motifs
	To run use: ./bc moitfSeqFile 
	To compile using: g++ branch_cut.cc -lglpk -o bc
	Args:
		motifDict: dict between motif IDs and seqs it occurs in
		tomtomDict: dict between motif Ids and a list of motif Ids that are similar to it
	Returns:
		depthDict: dictionary between the depth and a tuple of (motifIdList, list of sets of seqs per motif, newSeqsAdded)
	"""
	print 'Branch and cut sequence coverage with depth:', depth
	#get the path of the module	
	branch_path = os.path.realpath(__file__)
	#remove the name of the module from the end of it. 
	split = branch_path.split('/')
	branch_path = '/'.join(split[:len(split)-1])
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
		branch_path = branch_path + '/bc ' + motifSeqFileName
		branchOutFile = jid + '_branch_out_depth_' + str(depth)
		fileList.append(branchOutFile)
		cmd = branch_path + ' > ' + branchOutFile
		print 'branch cmd:', cmd 
		try:
			os.system(cmd)
		except:
			print 'branch execution failed. Exiting'
		#read and process the branch file and return a list of motif IDs ordered as output by the algthm
		motifIdList = processBranchoutput(branchOutFile)
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
