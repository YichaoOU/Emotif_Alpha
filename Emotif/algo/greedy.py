from __future__ import division
"""
greedy.py -- greedy algorithm
greedy.py -- greedy algorithm

Changelog:
    greedy.py ; started at 2014-09-24 16:44:57 EDT by Liang Chen
"""
import sys
from copy import deepcopy
##

def findMin(Sdict, R, w):
	"""Find the min cost subset. the subset that covers the most seqiences is the best.
	
	Args:
		Sdict: dict between motif Ids and sequences the motif occurs in
		R: the list of sequences not covered yet
		w: list of weights assigned to each subset/motif. For greedy here all set to 1 (unicost set cover)
	
	Returns:
		subset with lowest cost and cost of it and motif ID of this subset
	
	"""
	minCost = 99999 
	minElement = -1
	for motifId, sSet in Sdict.items():
		if len(sSet.intersection(R)) == 0:#no new 
			continue
		else:
			#print '\tId:', motifId, sSet
			tmpSet = sSet.intersection(R)
			#print '\ttmpSet:', tmpSet
			cost = w[motifId-1]/(len(sSet.intersection(R)))
			if cost < minCost:
				minCost = cost
				minElement = motifId 
	
	#print 'miElem:', minElement,'Sdict[min]', Sdict[minElement]
	return Sdict[minElement], w[minElement-1], minElement
	


def callGreedy(Uset, Sdict):
	"""use the greedy algthm to find the smallest set cover
	
	Follow the pseudo code by the introduction to algorithms book and the wiki page. Very simple
	Just pick the set which covers the most at each iteration.
	Used an online implementation for now. This is greedy with no weights (unicost set cover) so make the weights
	vector all equal to one
	
	Args:
		Uset: set of seq Ids
		Sdict; dict between motif Ids and the set of sequences each motif
		
	Returns:
		motifIdList: list of motif Ids with smallest set coverage
		C: the seq sets found, 
	"""
	print "greedy.py::call_greedy()"
	#see how many motifs (subsets) we have
	numMotifs = len(Sdict)
	#make a list of weigts equal to zero
	wghtList  = [1] * numMotifs
    
   
	R = Uset
	C = []
	costs = []
	counter = 0
	motifIdList = []
	while len(R) != 0:
		S_i, cost, motifId = findMin(Sdict, R, wghtList)
		C.append(S_i)
		motifIdList.append(motifId)
		R = R.difference(S_i)
		costs.append(cost)
		counter += 1
	
	#print 'Cover:', C
	print 'motif ID list:', motifIdList
	#print "Total Cost: ", sum(costs), costs
	return motifIdList, C
    

def findMinDepth(Sdict, R, w, foundIdList):
	"""Find the min cost subset. the subset that covers the most seqiences is the best.
	This has a new list attached to it which is the list of IDs already found
	Args:
		Sdict: dict between motif Ids and set of sequences the motif occurs in
		R: the list of sequences not covered yet
		w: list of weights assigned to each subset/motif. For greedy here all set to 1 (unicost set cover)
	
	Returns:
		subset with lowest cost and cost of it and motif ID of this subset
	
	"""
	minCost = 99999 
	minElement = -1
	newSeqsAdded = 0
	for motifId, sSet in Sdict.items():
		#skip if motif already reported
		if motifId in foundIdList:
			continue
		if len(sSet.intersection(R)) == 0:#no new 
			continue
		else:
			#print '\tId:', motifId, sSet
			#tmpSet = sSet.intersection(R)
			cost = w[motifId-1]/(len(sSet.intersection(R)))
			if cost < minCost:
				minCost = cost
				minElement = motifId
				newSeqsAdded = len(sSet.intersection(R)) 
	
	if minElement == -1:
		return list(), list(), -1, list()
	#print 'miElem:', minElement,'Sdict[min]', Sdict[minElement]
	return Sdict[minElement], w[minElement-1], minElement, newSeqsAdded


def callGreedyDepth(Uset, Sdict, depth):
	"""use the greedy algthm to find the smallest set cover with depth. For each depth discovered remove the motif IDs found
	for that dpeth then do the search again
	
	Args:
		Uset: set of seq Ids
		Sdict; dict between motif Ids and the set of sequences each motif
		
	Returns:
		motifIdList: list of motif Ids with smallest set coverage
		C: the seq sets found, 
	"""
	print 'Greedy sequence coverage with depth:', depth
	#see how many motifs (subsets) we have
	numMotifs = len(Sdict)
	#make a list of weigts equal to zero
	wghtList  = [1] * numMotifs
	#make a list that stores motif IDs found so far
	foundIdList = []
	#make a dpth dictionary between the depth and a tuple of (motifIdList, list of sets of seqs per motif)
	depthDict = {}
    #loop thru the depths 
	for i in range(1,depth+1):
		print 'depth:', i
		R = Uset
		C = []
		costs = []
		counter = 0
		motifIdList = []
		newSeqsAddedList = []
		while len(R) != 0:
			S_i, cost, motifId, newSeqsAdded = findMinDepth(Sdict, R, wghtList,foundIdList)
			C.append(S_i)
			motifIdList.append(motifId)
			R = R.difference(S_i)
			newSeqsAddedList.append(newSeqsAdded)
			costs.append(cost)
			counter += 1
		
		#print 'Cover:', C
		print 'motif ID list:', motifIdList
		#update the found list
		for newId in motifIdList:
			foundIdList.append(newId)
		#make a tuple and store into depth dict
		depthDict[i] = (motifIdList, C, newSeqsAddedList)
		#print "Total Cost: ", sum(costs), costs
		#return motifIdList, C
	
	#return the depth dictionary
	return depthDict




def findMinDepthRefine(Sdict, R, w, foundIdList):
	"""Find the min cost subset. the subset that covers the most seqiences is the best.
	This has a new list attached to it which is the list of IDs already found
	Args:
		Sdict: dict between motif Ids and sequences the motif occurs in
		R: the list of sequences not covered yet
		w: list of weights assigned to each subset/motif. For greedy here all set to 1 (unicost set cover)
	
	Returns:
		subset with lowest cost and cost of it and motif ID of this subset
	
	"""
	minCost = 99999 
	minElement = -1
	newSeqsAdded = 0
	for motifId, sSet in Sdict.items():
		#skip if motif already reported or is simialr to previous motifs
		if motifId in foundIdList:
			continue
		if len(sSet.intersection(R)) == 0:#no new 
			continue
		else:
			#print '\tId:', motifId, sSet
			#tmpSet = sSet.intersection(R)
			cost = w[motifId-1]/(len(sSet.intersection(R)))
			if cost < minCost:
				minCost = cost
				minElement = motifId
				newSeqsAdded = len(sSet.intersection(R))  
	
	if minElement == -1:
		return list(), list(), -1, list()
	#print 'miElem:', minElement,'Sdict[min]', Sdict[minElement]
	return Sdict[minElement], w[minElement-1], minElement, newSeqsAdded


def callGreedyDepthRefine(Uset, Sdict, depth, tomtomDict, motifIdDict, idMotifDict, filterMinNumSeq):
	"""use the greedy algthm to find the smallest set cover with depth. For each depth discovered remove the motif IDs found
	for that dpeth then do the search again
	Here with motif refinement applied such that simialr motifs are removed 
	Args:
		Uset: set of seq Ids
		Sdict: dict between motif Ids and the set of sequences each motif
		depth: coverage depth
		tomtomDict: dict between motif Ids and a list of motif Ids that are similar to it
		motifIdDict: dict between motif name and ID (as found by processMotifHitFile)
		idMotifDict: dict between motif id and motif name (as found by processMotifHitFile)
	Returns:
		depthDict: dictionary between the depth and a tuple of (motifIdList, list of sets of seqs per motif)
	"""
	print 'Greedy sequence coverage with depth:', depth
	#see how many motifs (subsets) we have
	numMotifs = len(Sdict)
	#make a list of weigts equal to 1
	wghtList  = [1] * numMotifs
	#make a depth dictionary between the depth and a tuple of (motifIdList, list of sets of seqs per motif)
	depthDict = {}
	#ID of motifs found per depth
	motifIdList = []
    #loop thru the depths 
	for i in range(1,depth+1):
		print 'depth:', i
		#make a list that stores motif IDs found so far for depth coverage
		foundIdList = []
		#update the found list
		for newId in motifIdList:
			foundIdList.append(newId)
		print 'foundIdList:', foundIdList
		motifIdList = []
		R = Uset
		C = []
		costs = []
		newSeqsAddedList = []
		while len(R) != 0:
			#S_i is a set
			S_i, cost, motifId, newSeqsAdded = findMinDepthRefine(Sdict, R, wghtList, foundIdList)
			#check if nothing found so go to next depth
			if motifId != -1:
				#remove similar motifs to this one 
				for val in tomtomDict:
					if isinstance(val,str):
						tomIdType = 'str'
					if isinstance(val,int):
						tomIdType = 'integer'
					break
					
				motifName = idMotifDict[motifId]
				for matchName in tomtomDict[motifName]:
					if matchName not in motifIdDict:
						continue
					matchId = motifIdDict[matchName]
					foundIdList.append(matchId)
				
				#check how many new seqs added
				print 'newAdded:', newSeqsAdded,'filter:', filterMinNumSeq
				if newSeqsAdded > filterMinNumSeq:
					C.append(S_i)
					newSeqsAddedList.append(newSeqsAdded)
					R = R.difference(S_i)
					costs.append(cost)
					motifIdList.append(motifId)
				else:
					break
					
			else:
				break
		
		#print 'Cover:', C
		#make a tuple and store into depth dict
		print 'after list:', motifIdList
		depthDict[i] = (motifIdList, C, newSeqsAddedList)
		#print "Total Cost: ", sum(costs), costs
		#return motifIdList, C
	
	#return the depth dictionary
	return depthDict

def callGreedy_MultiCover(Uset, Sdict, depth, tomtomDict, motifIdDict, idMotifDict, filterMinNumSeq, seqIdDict, idSeqDict, resultsFileName
				, foreHitDict, backHitDict, totalForeNumSeqs, totalBackNumSeqs):
	"""
	Args:
		Uset: set of seq Ids
		Sdict: dict between motif Ids and the set of sequences each motif occurs in
		depth: coverage depth
		tomtomDict: dict between motif Ids and a list of motif Ids that are similar to it
		motifIdDict: dict between motif name and ID (as found by processMotifHitFile)
		idMotifDict: dict between motif id and motif name (as found by processMotifHitFile)
	"""
	#motifDict, minForeCov, foreHitDict, backHitDict, totalForeNumSeqs, totalBackNumSeqs, outFileName
	print 'in depth function depth:', depth
	
	#go thru all the sequences and make a class to store info for it, a dict between seq Id and Seq object
	seqObjDict = {}
	seqObjDict = makeSeqObjs(Uset, seqIdDict, idSeqDict, depth)
	
	#see how many motifs (subsets) we have
	numMotifs = len(Sdict)
	print 'numMotifs:', numMotifs
	#make a list of weigts equal to zero
	wghtList  = [1] * numMotifs
	print 'filterMinNumSeq:', filterMinNumSeq
	R = Uset
	C = []
	j = 0
	motifIdList = []
	newSeqsAddedList = []
	foundIdList = []
	while len(R) != 0 and j < len(Sdict):
		S_i, motifId, newSeqsAdded, seqCovList = findMax(Sdict, R, foundIdList, idMotifDict)
		j += 1
		#if motifId in motifIdList:
				#continue
		if motifId != -1:
			if newSeqsAdded < filterMinNumSeq:
				foundIdList.append(motifId)
				continue
			#remove similar motifs to this one 
			for val in tomtomDict:
				if isinstance(val,str):
					tomIdType = 'str'
				if isinstance(val,int):
					tomIdType = 'integer'
				break
				
			motifName = idMotifDict[motifId]
			for matchName in tomtomDict[motifName]:
				if matchName not in motifIdDict:
					continue
				matchId = motifIdDict[matchName]
				foundIdList.append(matchId)
			
			
			C.append(S_i)
			motifIdList.append(motifId)
			foundIdList.append(motifId)
			#go thru the sequences and decrement their depth
			print 'len seqCovList:', len(seqCovList)
			for seqId in seqCovList:
				seqObjDict[seqId].depthCount = seqObjDict[seqId].depthCount - 1
				#print 'd:', seqObjDict[seqId].depthCount
				if seqObjDict[seqId].depthCount == 0:
					tmpSet = set()
					tmpSet.add(seqId)
					R = R.difference(tmpSet)
					newSeqsAddedList.append(newSeqsAdded)
		else:
			break
	
	#print 'Cover:', C
	print 'motif ID list:', motifIdList
	
	#check the counts of the sequences:
	#for seqId in seqObjDict:
		#print 'seqId:', seqId, idSeqDict[seqId], 'depth:', seqObjDict[seqId].depthCount
	
	#print "Total Cost: ", sum(costs), costs
	depthDict = {}
	for val in motifIdList:
		print val,idMotifDict[val]
		#for seqName in Sdict[val]:
			#print '\t', seqName
	depthDict[1] = (motifIdList, C, newSeqsAddedList)
	
	#write results
	selectMotifList = writeResults(motifIdList, resultsFileName, foreHitDict, backHitDict, totalForeNumSeqs, totalBackNumSeqs, idMotifDict)
	return selectMotifList
	
def writeResults(motifIdList, resultsFileName, foreHitDict, backHitDict, totalForeNumSeqs, totalBackNumSeqs, idMotifDict):
	
	outFile = open(resultsFileName, 'wb')
	outFile.write('#Name,foreground_motif_seqCount,foreground_motif_seqCoverage(%),foreground_seqs_added,cumulative_foreground_seq_cov(%)' 
		   + ',background_motif_seqCount,background_motif_seqCoverage(%),background_seqs_added,cumulative_background_seq_cov(%)' +'\n')
	
	
	foreAllSeqList = []
	backAllSeqList = []
	backCumCov = 0
	foreCumCovCov = 0
	selectMotifList = []
	for val in motifIdList:
		#print val,idMotifDict[val]
		motifName = idMotifDict[val]
		lineList = []
		foreSeqsAdded = 0
		backSeqsAdded = 0
		motifForeCov = 100*(len(foreHitDict[motifName])/totalForeNumSeqs)
		motifForeCov = round(motifForeCov, 1)
		#fore info
		for seqName in foreHitDict[motifName]:
			if seqName not in foreAllSeqList:
				foreAllSeqList.append(seqName)
				foreSeqsAdded += 1
		foreCumCov = round(100*(len(foreAllSeqList)/totalForeNumSeqs), 1)
		foreCov = round(100*(len(foreHitDict[motifName])/totalForeNumSeqs), 1)
		foreNumSeqs = len(foreHitDict[motifName])
		
		#Background info
		if motifName in backHitDict:
			for seqName in backHitDict[motifName]:
				if seqName not in backAllSeqList:
					backAllSeqList.append(seqName)
					backSeqsAdded += 1
			backCumCov = round(100*(len(backAllSeqList)/totalBackNumSeqs), 1)
			backCov = round(100*(len(backHitDict[motifName])/totalBackNumSeqs), 1)
			backNumSeqs = len(backHitDict[motifName])
		else:
			backCov = 0
			backNumSeqs = 0
		
		selectMotifList.append(motifName)
		lineList = [motifName, str(foreNumSeqs), str(foreCov), str(foreSeqsAdded), str(foreCumCov), str(backNumSeqs), str(backCov), str(backSeqsAdded)
		,str(backCumCov)]
		lineStr = ','.join(lineList)
		outFile.write(lineStr + '\n')
		
	
	outFile.close()	
	print 'number of selected motifs:', len(selectMotifList)
	return selectMotifList
def makeSeqObjs(seqSet, seqIdDict, idSeqDict, depth):
	""" 
	Args:
		seqSet: set of seq IDs
	"""
	seqObjDict = {}
	for seqId in seqSet:
		seqName = idSeqDict[seqId]
		#print seqId, seqName
		#make an object
		seqObjDict[seqId] = Seq(seqId, seqName, depth)
		
		
	print 'len of seq obj dict:', len(seqObjDict)
	return seqObjDict

class Seq:
	"""Class for each sequence to keep count of how many times it is covered """
	def __init__(self, seqId, seqName, depthCount):
		self.seqId = seqId
		self.seqName = seqName
		self.depthCount = depthCount

def findMax(Sdict, R, foundIdList, idMotifDict):
	"""Find the motif(subset) that covers the most uncovered sequences.
	
	Args:
		Sdict: dict between motif Ids and sequences the motif occurs in
		R: the list of sequences not covered yet
		foundIdList: list of motifs simialr to previous motifs reported
		
	Returns:
		subset with which covers the most and motif ID of this subset
	"""
	maxCov = 0 
	newSeqsAdded = 0
	maxCovMotif = -1
	#list which stores the seqs covered by the motif
	seqCovList = []
	for motifId, sSet in Sdict.items():
		#skip if motif already reported or is simialr to previous motifs
		if motifId in foundIdList:
			continue
		if len(sSet.intersection(R)) == 0:#no new 
			continue
		else:
			#print '\tId:', motifId, idMotifDict[motifId], sSet
			#tmpSet = sSet.intersection(R)
			cov = len(sSet.intersection(R))
			if cov > maxCov:
				maxCov = cov
				maxCovMotif = motifId
				newSeqsAdded = len(sSet.intersection(R))
				seqCovList = list(sSet.intersection(R))  
	
	
	if maxCovMotif == -1:
		return list(), -1, list(), list()
	#print 'miElem:', minElement,'Sdict[min]', Sdict[minElement]
	#Sdict[minElement] is the list of seqs covered by motif maxElement 
	print maxCovMotif, idMotifDict[maxCovMotif], 'len seqCovList:', len(seqCovList), 'cov:', maxCov, 'newAdded:', newSeqsAdded
	return Sdict[maxCovMotif], maxCovMotif, newSeqsAdded, seqCovList	
		
def main(args):
    pass

if( __name__ == "__main__" ):
    main(sys.argv)
    pass
else:
    pass

