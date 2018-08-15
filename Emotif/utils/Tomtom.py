# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
from __future__ import division

import os
import re

def callTomtom(queryPwmFile, databasePwmFile, paraDict, outDir):
	"""Call the tomtom executable
	USAGE: tomtom [options] <query file> <target file>+

   Options:
     -o <output dir>		name of directory for output files;
				will not replace existing directory
     -oc <output dir>		name of directory for output files;
				will replace existing directory
     -bfile <background file>	name of background file
     -m <id>	use only query motifs with a specified id; may be repeated
     -mi <index>	use only query motifs with a specifed index; may be repeated
     -thresh <float>		significance threshold; default: 0.5
     -evalue		use E-value threshold; default: q-value
     -dist allr|ed|kullback|pearson|sandelin|blic1|blic5|llr1|llr5
				distance metric for scoring alignments;
				default: pearson
     -internal			only allow internal alignments;
				default: allow overhangs
     -min-overlap <int>		minimum overlap between query and target;
				default: 1
     -query-pseudo <float>	default: 0.0
     -seed <int>
     -incomplete-scores		ignore unaligned columns in computing scores      
				default: use complete set of columns
     -target-pseudo <float>	default: 0.0
     -text			output in text format (default is HTML)
     -png			create PNG logos
				default: don't create PNG logos
     -eps			create EPS logos
				default: don't create EPS logos
     -no-ssc			don't apply small-sample correction to logos
				default: use small-sample correction
     -verbosity [1|2|3|4]	default: 2

	default from website: 
	tomtom -no-ssc -oc . -verbosity 1 -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10 pwm_test_2 db/JASPAR_CORE_2014.meme
	
	Args:
		paraDict: dict of parameters for tomtom
		queryFile: PWM file of queries to serach for
		databaseFile: database of PWMs to serach against
	"""
	#get the path of the module	
	# tomtom_path = os.path.realpath(__file__)
	#remove the name of the module from the end of it. Replace won't work since w ehave DME.py and DME.pyc
	# split = tomtom_path.split('/')
	# tomtom_path = '/'.join(split[:len(split)-1])
	tomtom_path = 'tomtom'
	command = (tomtom_path + ' ' + ' -oc ' + outDir + ' --verbosity 1' + 
		' -no-ssc' + ' -min-overlap 5 -dist pearson -evalue -thresh ' + str(paraDict['evalue']) + ' '
		+ queryPwmFile + ' ' + databasePwmFile  )
	print 'Tomtom commnad:', command
	try:
		os.system(command)
	except:
		print 'Tomtom execution failed. Exiting'
		exit()


def removeMatches(inFile):
	"""remove matches from tomtom. make a dict between Ids and remove the matches
	
	Args:
		The tomtom.txt file
	
	"""
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
			query = int(lineSplit[0])
			target = int(lineSplit[1])
			eVal = float(lineSplit[4])
			overlap = int(lineSplit[6])
			strand = lineSplit[9]
			#if strand != '+':
				#continue
			if query not in matchDict:
				matchDict[query] = []
			if target not in matchDict[query]:
				matchDict[query].append(target)
	
	#check the match Dict and pick the motif with highest rank (score) before other motifs
	removeIdList = []
	selectIdList = []
	for queryId in sorted(matchDict.iterkeys(), key=int):
		#if the motif is in the remove list skip it
		if queryId in removeIdList:
			continue
		for targetId in matchDict[queryId]:
			if targetId == queryId:
				continue
			if targetId not in removeIdList:
				removeIdList.append(targetId)
		selectIdList.append(queryId)
		
			
	return selectIdList
			


def parse(inFile):
	"""Parse the tomtom output file
	Args:
		tomtom.txt output file
	Returns:
		A dict between motif Ids and list of matched motif Ids
	"""	
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
				matchDict[query] = []
			if target not in matchDict[query] and target != query:
				matchDict[query].append(target)
	
	
	
	return matchDict


def makeReport(matchDict, motifDict, targetName, outFileName):
	"""Make a summary report that shows the motif names and the motifs matched to it and wheter the target TF was found and
	and at what depth
	Args:
		matchDict: dict between motif Ids and a list of names/ids that match this motif
		motifDict: the motif dict between motif Ids and MyMotif objects
		targetName: name of a specific target TF we searching for
		outFileName: name of the outptu file to write to
	Returns:
	"""
	outFile = open(outFileName, 'wb')
	#print the header
	outFile.write('#ID,name,Coverage(%),Depth,Matches' + '\n')
	foundCounter = 0
	depthList = []
	for motifId in matchDict:
		lineList = []
		motifObj = motifDict[motifId]
		matches = '///'.join(matchDict[motifId])
		#check if any of the matches match the targetName supplied
		matchTarget = 'No'
		for matchName in matchDict[motifId]:
			if re.search(targetName, matchName):
				foundCounter += 1
				matchTarget = 'Yes'
				if motifObj.depth not in depthList:
					depthList.append(motifObj.depth)
		lineList = [str(motifId), motifObj.regExp,str(motifObj.seqCov), str(motifObj.depth), matches, matchTarget]
		line = ','.join(lineList)
		outFile.write( line + '\n')
		
	outFile.write('#######\n');
	outFile.write('num_matches:' + str(foundCounter) + '\n')
	#sort the list
	depthList.sort(key=int)
	depthLine = [str(i) for i in depthList]
	depthLine = '-'.join(depthLine)
	outFile.write('match_depth(s):' + depthLine + '\n')
	
	
	outFile.close()
	

def makeReportScan(matchDict, motifDict, targetName, outFileName):
	"""Make a summary report that shows the motif names and the motifs matched to it and wheter the target TF was found and
	and at what depth
	Args:
		matchDict: dict between motif Ids and a list of names/ids that match this motif
		motifDict: the motif dict between motif Ids and MyMotif objects
		targetName: name of a specific target TF we searching for
		outFileName: name of the outptu file to write to
	Returns:
	"""
	outFile = open(outFileName, 'wb')
	#print the header
	outFile.write('#ID,Matches,Match_names' + '\n')
	foundCounter = 0
	depthList = []
	for motifId in matchDict:
		lineList = []
		#motifObj = motifDict[motifId]
		matches = '///'.join(matchDict[motifId])
		#check if any of the matches match the targetName supplied
		matchTarget = 'No'
		for matchName in matchDict[motifId]:
			if re.search(targetName, matchName):
				foundCounter += 1
				matchTarget = 'Yes'
				#if motifObj.depth not in depthList:
					#depthList.append(motifObj.depth)
		#lineList = [str(motifId), motifObj.regExp,str(motifObj.seqCov), str(motifObj.depth), matches, matchTarget]
		lineList = [str(motifId),matchTarget,matches]
		line = ','.join(lineList)
		outFile.write( line + '\n')
		
	outFile.write('#######\n');
	outFile.write('num_matches:' + str(foundCounter) + '\n')
	
	
	
	outFile.close()


def main(args):
    pass

##
if( __name__ == '__main__' ):
    main(sys.argv)


