from __future__ import division
import os
import sys
import shutil
import re
# from xml.dom import minidom as mnd
from utils.general_utils import makeMotifDictFromPwmFile as mdp
from utils.general_utils import PWMtoPFM
from Bio import SeqIO as io
import matplotlib as mpl
mpl.use('Agg')
from motif_scanning import *
from motif_scanning import *
from motif_selection import *
from motif_output import *
def run_testing_SE_SP_ACC(jid,confDict):

	

	
	### 2-25-2016 added
	testing_plots_files = []
	testing_plots_folder = jid + "_testing"
	os.makedirs(testing_plots_folder)
	
	output_file = jid + "_testing_result.csv"
	f = open(output_file,"wb")
	testing_plots_files.append(output_file)
	
	
	SolPwmFile = confDict['input']['pwm_selected']
	pos_seq = confDict['input']['pos_testing_seq']
	neg_seq = confDict['input']['neg_testing_seq']
	# folder movement is coded in run_fimo()
	run_fimo_testing(jid,confDict,testing_plots_files)
	# confDict['input']['testing_pos_mhit'] = jid + '_pos.mhit'
	# confDict['input']['testing_neg_mhit'] = jid + '_neg.mhit'
	pos_mhit = confDict['input']['testing_pos_mhit']
	neg_mhit = confDict['input']['testing_neg_mhit']
	# testing_plots_files.append(file)
	posCovered = readHitFile2NumSeqs(pos_mhit)
	negCovered = readHitFile2NumSeqs(neg_mhit)
	posNumSeqs = findNumSeqs(pos_seq)
	negNumSeqs = findNumSeqs(neg_seq)
	tp = len(posCovered)
	fp = len(negCovered)
	# fp = len([mydict[x] for x in mykeys])

	tn = negNumSeqs - fp
	fn = posNumSeqs - tp
	
	[se,sp,acc,fdr,ppv] = cal_metrics(tp,tn,fp,fn)
	print >>f,(str(se)+","+str(sp)+","+str(acc))
	f.close()
	# testing_plots_files.append(pos_mhit)
	# testing_plots_files.append(neg_mhit)
	for outFile in testing_plots_files:
		print outFile
		shutil.move(outFile, testing_plots_folder)
	
	return testing_plots_folder

def run_testing_SE_SP_ACC_genetic(jid,confDict,label):

	

	
	### 2-25-2016 added
	testing_plots_files = []
	testing_plots_folder = jid + "_testing_" + label
	os.makedirs(testing_plots_folder)
	
	output_file = jid + "_testing_result.csv"
	f = open(output_file,"wb")
	testing_plots_files.append(output_file)
	
	
	SolPwmFile = confDict['input']['pwm_selected']
	pos_seq = confDict['input']['pos_testing_seq']
	neg_seq = confDict['input']['neg_testing_seq']
	# folder movement is coded in run_fimo()
	run_fimo_testing(jid,confDict,testing_plots_files)
	# confDict['input']['testing_pos_mhit'] = jid + '_pos.mhit'
	# confDict['input']['testing_neg_mhit'] = jid + '_neg.mhit'
	pos_mhit = confDict['input']['testing_pos_mhit']
	neg_mhit = confDict['input']['testing_neg_mhit']
	# testing_plots_files.append(file)
	posCovered = readHitFile2NumSeqs(pos_mhit)
	negCovered = readHitFile2NumSeqs(neg_mhit)
	posNumSeqs = findNumSeqs(pos_seq)
	negNumSeqs = findNumSeqs(neg_seq)
	tp = len(posCovered)
	fp = len(negCovered)
	# fp = len([mydict[x] for x in mykeys])

	tn = negNumSeqs - fp
	fn = posNumSeqs - tp
	
	[se,sp,acc,fdr,ppv] = cal_metrics(tp,tn,fp,fn)
	print >>f,(str(se)+","+str(sp)+","+str(acc))
	f.close()
	# testing_plots_files.append(pos_mhit)
	# testing_plots_files.append(neg_mhit)
	for outFile in testing_plots_files:
		print outFile
		shutil.move(outFile, testing_plots_folder)
	
	return testing_plots_folder


def read_score_list(pos_score,neg_score,jid):
	pos_neg = []
	dict = {}
	lines = [line.rstrip('\n') for line in open(pos_score)]
	count = 0
	label = "positive_"
	for line in lines:
		count += 1
		key = label + str(count)
		dict[key] = float(line)
	lines = [line.rstrip('\n') for line in open(neg_score)]
	p_count = count
	count = 0
	label = "negative_"
	for line in lines:
		count += 1
		key = label + str(count)
		dict[key] = float(line)
	for w in sorted(dict, key=dict.get, reverse=True):
		if "negative" in w:
			pos_neg.append("n")
		else:
			pos_neg.append("p")
	n_count = count
	
	AUC = roc_auc(pos_neg,jid,p_count,n_count)
	# command = "rm -rf " + output_folder
	# os.system(command)
	
	
	# liyc@liyc-OU:/home/working/program/DiMO_v1.6/R$ python read_convert.py pos_score.txt neg_score.txt 
	# 0.775637504067
	# test my AUC calculation with DiMO AUC, they are the same
	# liyc@liyc-OU:/home/working/program/DiMO_v1.6/R$ Rscript test1.py Drosophila_bcd_1_training_positive.fasta Drosophila_bcd_1_training_negative.fasta bcd_1.pwm 
	# [1] 0.778371
	# The 0.003 difference may due to the order of "P" "N" when they are the same score
	# DiMO AUC is in scorePerceptron.R
	return AUC

def run_roc_auc_array(jid_array,pos_seq,neg_seq,pwm_file_array):
	length = len(jid_array)
	AUC_array = []
	for i in range(length):
		pwm_file = pwm_file_array[i]
		jid = jid_array[i]
	
		AUC = mast_mapping(jid,pwm_file,pos_seq,neg_seq)
		# the jid.png is a file
		AUC_array.append(AUC)
	return AUC_array
	
def run_roc_auc(jid,pos_seq,neg_seq,pwm_file):
	# length = len(jid_array)
	# AUC_array = []
	# for i in range(length):
	# pwm_file = pwm_file_array[i]
	# jid = jid_array[i]
	fileList = []
	[p_size,n_size] = re_output_seq(pos_seq,neg_seq) # pos.fa neg.fa
	R_script = "/usr/local/lib/python2.7/dist-packages/emoti/DiMO/test.py"
	dict = mdp(pwm_file)
	pfm_array_file = "pfm_list.txt"
	f = open(pfm_array_file,"wb")
	fileList.append(pfm_array_file)
	for k in dict.keys():
		
		m_obj = dict[k]
		# print "#############m_obj.Id"
		# print "#############m_obj.Id"
		# print k
		# print m_obj.Id
		# print m_obj.regExp
		# print m_obj.pwmLines
		pfm = PWMtoPFM(m_obj.pwmLines,k)
		fileList.append(pfm)
		print >>f,pfm
	f.close()
	R_command = "Rscript " + R_script + ' ' + \
	str(p_size) + " " + str(n_size) + ''
	print R_command
	os.system(str(R_command))
	AUC = read_score_list("pos_score.txt","neg_score.txt",jid)
	fileList.append("pos_score.txt")
	fileList.append("neg_score.txt")
	fileList.append("pos.fa")
	fileList.append("neg.fa")
	return AUC
	
def run_roc_auc_fimo(jid,p_count,n_count,m_name,dict):
	# fileList = []
	seq_score = {}
	pos_neg=[]
	for m in m_name:
		for s in dict[m].keys():
			seq_score[s] = dict[m][s]
	for w in sorted(seq_score, key=seq_score.get):
		# print w,seq_score[w]
		if "negative" in w:
			pos_neg.append("n")
		else:
			pos_neg.append("p")	
	
	AUC = roc_auc(pos_neg,jid,p_count,n_count)
	
	
	return AUC
		
def mast_mapping(jid,pwm_file,pos_seq,neg_seq):
	# command = "cat "+pos_seq+" "+neg_seq +" > mast.seqs"
	# os.system(command)
	output_folder = jid + "_mast_out"
	p_count, n_count = combine_pos_neg_seqs(pos_seq,neg_seq)
	command = "mast " + pwm_file + " mast.seqs -ev 100000 -notext -nohtml -oc " + output_folder
	os.system(command)
	os.system("rm mast.seqs")
	xml = mnd.parse("./"+output_folder+"/mast.xml")
	seqs = xml.getElementsByTagName('sequence')
	pos_neg=[]
	for i in range(p_count+n_count):
		if "negative" in seqs[i].attributes['name'].value:
			pos_neg.append("n")
		else:
			pos_neg.append("p")	
	AUC = roc_auc(pos_neg,jid,p_count,n_count)
	command = "rm -rf " + output_folder
	os.system(command)
	# delete output_folder later
	
	return AUC
	
def combine_pos_neg_seqs(pos_seq,neg_seq):
	seq_hash={}
	num_pos = 0
	num_neg = 0
	f=open(pos_seq,"rU")
	record = io.parse(f,"fasta")
	for r in record:
		seq_hash[str(r.id)+"_positive"]=str(r.seq)
		num_pos += 1
	f.close()
	f=open(neg_seq,"rU")
	record = io.parse(f,"fasta")
	for r in record:
		seq_hash[str(r.id)+"_negative"]=str(r.seq)
		num_neg += 1
	f.close()
	f=open("mast.seqs","wb")
	for k in seq_hash:
		print >>f,">"+k
		print >>f,seq_hash[k]
	f.close()
	return [num_pos,num_neg]
	
def re_output_seq(pos_seq,neg_seq):
	seq_hash={}
	num_pos = 0
	num_neg = 0
	f=open(pos_seq,"rU")
	record = io.parse(f,"fasta")
	for r in record:
		seq_hash[str(r.id)+"_positive"]=str(r.seq).upper().replace("N","A")
		num_pos += 1
	f.close()
	f=open("pos.fa","wb")
	for k in seq_hash:
		print >>f,">"+k
		print >>f,seq_hash[k]
	f.close()
	seq_hash={}
	f=open(neg_seq,"rU")
	record = io.parse(f,"fasta")
	for r in record:
		# DiMO can't accept N letter, DiMO needs all caps
		seq_hash[str(r.id)+"_negative"]=str(r.seq).upper().replace("N","A")
		num_neg += 1
	f.close()
	f=open("neg.fa","wb")
	for k in seq_hash:
		print >>f,">"+k
		print >>f,seq_hash[k]
	f.close()
	return [num_pos,num_neg]
	
	
def roc_auc(pos_neg,file_name,p_count,n_count):
	SE=[0.0]
	FP=[0.0]
	count_n=0
	auc_value = 0.0
	for i in range(len(pos_neg)):
		if pos_neg[i] == "n":
			count_n+=1
			MF=i+1-count_n
			MB=count_n
			SE.append(MF/float(p_count))
			FP.append(MB/float(n_count))
	#print SE
	#print FP
	SE.append(1.0)
	FP.append(1.0)
	for i in range(1,len(SE)):
			auc_value += (SE[i]+SE[i-1])*(FP[i]-FP[i-1])	
	import matplotlib.pyplot as plt
	plt.plot(FP,SE)
	title = "ROC curve " + file_name
	plt.title(title)
	plt.ylabel('sensitivity')
	plt.xlabel('1 - specificity')
	plt.axis([0, 1, 0, 1])
	out_png = file_name + ".png"
	plt.savefig(out_png)
	plt.close()
	return auc_value/2
	
	

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



def cal_metrics(tp,tn,fp,fn):
	se = tp / (tp + fn)
	sp = tn / (fp + tn)
	acc = (tp + tn) / (tp + fn + fp + tn)
	fdr = fp / (tp + fp)
	ppv = tp / (tp + fp)
	return [se,sp,acc,fdr,ppv]

def readHitFile2NumSeqs(hitFile):
	"""Read the hit file and make a dict between motifnames and their seqs """
	seqList =[]
	# motifName = ''
	with open(hitFile, 'rb') as handler:
		for line in handler:
			line = line.strip()
			if re.search(r'>', line):
				# motifName = line[1:]
				# if motifName not in motifDict:
					# motifDict[motifName] = []
				continue
			seqList.append(line)
	
	# print 'len of motifDict:', len(motifDict)
	return list(set(seqList))


