# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
from __future__ import division
"""
File for MyMotif class which is the more general version of the BioPython motif obect
and other motif related functions
"""

class MyMotif:
	"""general motif class
	
	"""
	def __init__(self):
		self.Id = 0 #motif id, an int value from 1
		self.regExp = ''
		self.bioMotifObj = '' #the BioPython motif object
		self.logoFileName = ''
		self.foreSeqList = []#fore=foreground. list of seq names where the motifs discoverd by the tool occur in. This should be as reported by the tool itself
		self.foreWordObjList = []#list of word objects that make up this motif. Of type MyWord
		self.foreSeqCov = 0
		self.foreDepth = 0

		#these below values are filled after running the motifs across the testing file if available
		self.seqList = []#list of seqs a motif occurs in
		self.seqIdList = []#list of seq Ids the motif occurs in 
		self.score = 0		
		self.seqCov = 0 #percent of sequences covered
		self.depth = 0 #the depth coverage of the motif
		
		#PWM lines as a list
		self.pwmLines = []
		

class MyWord:
	"""general word (site) class to store info for a specific word
	
	"""
	def __init__(self):
		self.string = ''#the word(site) string
		self.position = 0 #pos of the word
		self.seqName = '' #name of seq the word occurs in
		self.strand = '' #p for pos or n for neg
		self.motifId = 0#id of the motif that this word belongs to 


#####
#####		
def main(args):
    pass

##
if( __name__ == '__main__' ):
    main(sys.argv)


