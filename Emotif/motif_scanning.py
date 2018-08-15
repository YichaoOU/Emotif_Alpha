'''
input pwm files
output tab-deliminated location table (future work)
output motif hit file (.mhit)
9-29-2015
Yichao combines Rami's code

'''
import shutil
from utils import Fimo
from utils import general_utils
import re
import os
def run_fimo(jid, confDict,fileList):

	MS_folders = []
	MS_dir = jid+"_Motif_Scanning"
	os.makedirs(MS_dir)

	#get the PWM file to scan for 
	pwmFileName = confDict['input']['pwm_file']
	#find the total number of seqs in the testing file
	foreNumSeqs = general_utils.findNumSeqs(confDict['input']['pos_seq']) 
	if confDict['input']['neg_seq'] != 'none':
		backNumSeqs = general_utils.findNumSeqs(confDict['input']['neg_seq'])
	
	#read the PWM file and get the order of the motifs so to be used when reporting the results
	orderMotifList = getPwmOrder(pwmFileName)
	# fileList = []
	#run the FIMO tool for froeground file
	fimoParaDict = {}
	fimoParaDict['thresh'] = confDict['fimo']['pvalue']
	#fimo output directory name
	fimoForeOutDir = jid + '_fore_Fimo'
	Fimo.callFimo(confDict['input']['pos_seq'] , pwmFileName, fimoParaDict, fimoForeOutDir)
	MS_folders.append(fimoForeOutDir)
	#parse the fimo results folder
	strand = confDict['fimo']['strand']
	#fimoDict between motif ID and seq hits in the testing file
	foreFimoDict = Fimo.parseFimo(fimoForeOutDir+'/fimo.txt', strand)
	confDict['input']['pos_fimo'] = MS_dir+"/"+fimoForeOutDir+'/fimo.txt'
	#write the motifs in the format for sequence coverge algorithms 
	foreMotifSeqFileName = jid + '_pos.mhit'
	general_utils.writeMotifSeqFile(foreFimoDict, foreMotifSeqFileName)
	fileList.append(foreMotifSeqFileName)
	
	
	#run the FIMO tool for the background file if it exists
	if confDict['input']['neg_seq'] != 'none':
		#Background fimo output dir name
		fimoBackOutDir = jid + '_back_Fimo'
		Fimo.callFimo(confDict['input']['neg_seq'], pwmFileName, fimoParaDict, fimoBackOutDir)
		MS_folders.append(fimoBackOutDir)
		#parse the fimo results folder
		strand = confDict['fimo']['strand']
		#fimoDict between motif ID and seq hits in background file
		backFimoDict = Fimo.parseFimo(fimoBackOutDir + '/fimo.txt', strand)
		confDict['input']['neg_fimo'] = MS_dir+"/"+fimoBackOutDir + '/fimo.txt'
		#write the motifs in the format for sequence covearge algorithms
		backMotifSeqFileName = jid + '_neg.mhit'
		general_utils.writeMotifSeqFile(backFimoDict, backMotifSeqFileName)
		fileList.append(backMotifSeqFileName) 
	
	
	
	# Move fore and back folders into motif scanning folders
	for outFile in MS_folders:
		print outFile
		shutil.move(outFile, MS_dir)	
	fileList.append(MS_dir)

def run_fimo_testing(jid, confDict,fileList):
	jid =jid + "_testing"
	MS_folders = []
	MS_dir = jid+"_Motif_Scanning"
	os.makedirs(MS_dir)

	#get the PWM file to scan for 
	pwmFileName = confDict['input']['pwm_selected']
	#find the total number of seqs in the testing file
	foreNumSeqs = general_utils.findNumSeqs(confDict['input']['pos_testing_seq']) 
	if confDict['input']['neg_testing_seq'] != 'none':
		backNumSeqs = general_utils.findNumSeqs(confDict['input']['neg_testing_seq'])
	
	#read the PWM file and get the order of the motifs so to be used when reporting the results
	orderMotifList = getPwmOrder(pwmFileName)
	# fileList = []
	#run the FIMO tool for froeground file
	fimoParaDict = {}
	fimoParaDict['thresh'] = confDict['fimo']['pvalue']
	#fimo output directory name
	fimoForeOutDir = jid + '_fore_Fimo'
	Fimo.callFimo(confDict['input']['pos_testing_seq'] , pwmFileName, fimoParaDict, fimoForeOutDir)
	MS_folders.append(fimoForeOutDir)
	#parse the fimo results folder
	strand = confDict['fimo']['strand']
	#fimoDict between motif ID and seq hits in the testing file
	foreFimoDict = Fimo.parseFimo(fimoForeOutDir+'/fimo.txt', strand)
	confDict['input']['pos_fimo_testing'] = MS_dir+"/"+fimoForeOutDir+'/fimo.txt'
	#write the motifs in the format for sequence coverge algorithms 
	foreMotifSeqFileName = jid + '_pos.mhit'
	confDict['input']['testing_pos_mhit'] = foreMotifSeqFileName
	general_utils.writeMotifSeqFile(foreFimoDict, foreMotifSeqFileName)
	fileList.append(foreMotifSeqFileName)
	
	
	#run the FIMO tool for the background file if it exists
	if confDict['input']['neg_testing_seq'] != 'none':
		#Background fimo output dir name
		fimoBackOutDir = jid + '_back_Fimo'
		Fimo.callFimo(confDict['input']['neg_testing_seq'], pwmFileName, fimoParaDict, fimoBackOutDir)
		MS_folders.append(fimoBackOutDir)
		#parse the fimo results folder
		strand = confDict['fimo']['strand']
		#fimoDict between motif ID and seq hits in background file
		backFimoDict = Fimo.parseFimo(fimoBackOutDir + '/fimo.txt', strand)
		confDict['input']['neg_fimo_testing'] = MS_dir+"/"+fimoBackOutDir + '/fimo.txt'
		#write the motifs in the format for sequence covearge algorithms
		backMotifSeqFileName = jid + '_neg.mhit'
		confDict['input']['testing_neg_mhit'] = backMotifSeqFileName
		general_utils.writeMotifSeqFile(backFimoDict, backMotifSeqFileName)
		fileList.append(backMotifSeqFileName) 
	
	
	
	# Move fore and back folders into motif scanning folders
	for outFile in MS_folders:
		print outFile
		shutil.move(outFile, MS_dir)	
	fileList.append(MS_dir)




def getPwmOrder(pwmFileName):
	"""
	Read the PWm file and get the motifs order as in the file
	"""
	motifList = []
	with open(pwmFileName, 'rb') as handler:
		for line in handler:
			if re.search(r'MOTIF', line):
				split = line.split()
				motifId = split[1]
				motifList.append(motifId)
	
	print 'motifList:', motifList	
	return motifList
















