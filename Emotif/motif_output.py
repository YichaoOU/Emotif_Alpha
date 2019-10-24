from __future__ import division


#script to read the motifs of a specific solution and get its PWMs and make its logos and do scanning if needed
#Also make a latex file to make the tables and logos automatically 
import sys
import os 
import argparse
import shutil
import re
import copy
from copy import deepcopy as dp
from glob import glob
#add the util path to use util files in it
# sys.path.append('/home/rami/Documents/seq_cov/utils')
from utils import general_utils
from utils.Fimo import *
from Bio import SeqIO
import math
from utils.general_utils import makeMotifDictFromPwmFile as readPWM
from motif_evaluation import run_roc_auc_fimo
from motif_filtering import hit2hit
PATH = os.path.dirname(os.path.abspath(__file__))
# print "########################"
# print PATH

def run_motif_output(jid,confDict):
	fileList = []
	
	if confDict['job_type']['motif_output'] == 'true':		
		MOUT_folder = []	
		MOUT_dir = jid+"_Motif_Output"	
		os.makedirs(MOUT_dir)	
		if 	confDict['job_type']['motif_selection'] == 'true':	
			if confDict['motif_selection']['genetic_algo'] == 'true':			
				file,pdf = run_Rami_style(jid,confDict,"genetic_algo")		
				MOUT_folder.append(file)		
				fileList.append(pdf)		
			if confDict['motif_selection']['genetic_algo_multi'] == 'true':			
				file,pdf = run_Rami_style(jid,confDict,"genetic_algo_multiCover")		
				MOUT_folder.append(file)		
				fileList.append(pdf)		
			if confDict['motif_selection']['greedy'] == 'true':			
				greedy_result_file = confDict['input']['greedy_result_file']		
				file,pdf = run_Rami_style_greedy(jid,confDict,greedy_result_file)		
				MOUT_folder.append(file)		
				fileList.append(pdf)		
			if confDict['motif_selection']['greedy_multicover'] == 'true':			
				# latexFile = run_Rami_style_multicover		
				file,pdf = run_Rami_style_multicover(jid,confDict,"greedy_multicover_algo")		
				MOUT_folder.append(file)		
				# makeLatexTable(motifList, motifObjDict, latexFile)		
				fileList.append(pdf)		
		confDict['input']['pwm_selected']	= 	confDict['input']['all_pwm_file']
		file = Motif_discovery_accuracy_AUC(jid,confDict,"all_motif_ranking","")
		MOUT_folder.append(file)
		for outFile in MOUT_folder:	
			print outFile
			shutil.move(outFile, MOUT_dir)
		# file = single_motif_rank(jid,confDict,"test_rank")	
		fileList.append(MOUT_dir)	
	return fileList


def run_Rami_style_greedy(jid,confDict,greedy_result_file):
	confDict['input']['num_sols'] = confDict['motif_selection']['depth']
	confDict['input']['ga_solution'] = greedy_result_file
	label = "greedy_refine"
	file,pdf = run_Rami_style(jid,confDict,label)
	return file,pdf

def run_Rami_style(jid,confDict,label):

	
	jid = jid + "_" + label
	print "run_Rami_style jid : ",jid
	
	### 2-24-2016 added
	PCA_plots_files = []
	PCA_plots_folder = jid + "_final_selection_PCA"
	os.makedirs(PCA_plots_folder)
	
	output_folder = jid + "_rami_style_motif_output"
	os.makedirs(output_folder)
	fileList = [] 
	fileList.append(PCA_plots_folder)
	allPwmFile = confDict['input']['pwm_file']
	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	solFile = confDict['input']['ga_solution']
	num_sols = confDict['input']['num_sols']
	pos_mhit = confDict['input']['pos_mhit']
	neg_mhit = confDict['input']['neg_mhit']
	pos_fimo = confDict['input']['pos_fimo']
	neg_fimo = confDict['input']['neg_fimo']
	dict_msp = {}
	parseFIMO_forAUC(dict_msp,pos_fimo,"positive")
	parseFIMO_forAUC(dict_msp,neg_fimo,"negative")
	all_solution_motifs = {}
	# mdict = readPWM(allPwmFile)
	solNumList = range(1,int(num_sols)+1)
	mast_command = "mast motif_file " + pos_seq + " -o output_folder_mast -ev 9999 -notext -nostatus"
	if "greedy" in jid:
		solObjDict = parse_motif_selection_solution_greedy(solFile,solNumList,pos_mhit)
	else:
		solObjDict = parse_motif_selection_solution2(solFile,solNumList,pos_mhit)
	posNumSeqs = findNumSeqs(pos_seq)
	negNumSeqs = findNumSeqs(neg_seq)
	
	latexFileName = jid + '_latexTables.tex'
	fileList.append(latexFileName)
	latexFile = open(latexFileName, 'wb')
	header=r'''\documentclass{article}
	\usepackage{graphicx} 
	\usepackage{geometry}
	\geometry{margin=0in}
	'''
	latexFile.write(header)
	latexFile.write('\\begin{document}\n')
	#check it
	for solNum, solObj in solObjDict.items():
		print 'S:',solNum,solObj.solNum
		for motifName in solObj.motifObjDict:
			print '\t', motifName
			all_solution_motifs[motifName] = solObj.motifObjDict[motifName]
		#make a pwm file for this solution
		#out file for selected PWMs
		### Yichao Edit here
		mName = solObj.motifObjDict.keys()
		if len(mName) < 1:
			continue
		print mName
		print "###############    motif selection selected moitfs name ##########"
		# makeSolPwmFile(selectPwmFile, mName, mdict)
		### 2-26-2016 added 
		out_pos_hit = jid + '_sol_' + str(solNum) + "_final_selection_pos.mhit"
		out_neg_hit = jid + '_sol_' + str(solNum) + "_final_selection_neg.mhit"
		hit2hit(pos_mhit,neg_mhit,out_pos_hit,out_neg_hit,mName)
		confDict['input']['pos_mhit'] = out_pos_hit
		confDict['input']['neg_mhit'] = out_neg_hit
		file = single_motif_rank(jid,confDict,'sol_' + str(solNum)+"final_selection")
		PCA_plots_files.append(out_pos_hit)
		PCA_plots_files.append(out_neg_hit)
		PCA_plots_files.append(file)
		file = mhit2csv(jid,confDict,'sol_' + str(solNum)+"_final_selection_")
		PCA_plots_files.append(file)
		output_file = jid + '_sol_' + str(solNum) +"_"+"final_selection"+"_PCAplot.png"
		# print __file__
		R = "Rscript %s/pca_plot.R "%(PATH)+ file +" "+ output_file
		os.system(R)
		PCA_plots_files.append(output_file)
		confDict['input']['pos_mhit'] = pos_mhit
		confDict['input']['neg_mhit'] = neg_mhit
		
		
		### finish ###
		selectPwmFile = jid + '_sol_' + str(solNum) + '_select_motifs.pwm'
		
		makeSolPwmFile2(selectPwmFile, solObj.motifObjDict, allPwmFile)
		 
		fileList.append(selectPwmFile)
		
		#make the logos
		#make motif logos from the selected PWMs
		logoDirName = jid + '_sol_' + str(solNum) + '_logos'
		makeMotifLogoFromPwm(selectPwmFile, logoDirName)
		fileList.append(logoDirName)
		
		#make the latex file
		latexFile.write('\graphicspath{ {./'+ logoDirName +'/} }\n')
		latexFile.write('\clearpage \n')
		# print 'posCum:', solObjDict[solNum].posCumCov
		posCumNumSeqs = str(round(math.ceil(solObjDict[solNum].posCumCov*posNumSeqs), 1))
		# negCumNumSeqs = str(round(math.ceil(solObjDict[solNum].negCumCov*posNumSeqs), 1))
		# print 'posCumNumSeqs:', posCumNumSeqs
		
		negCumNumSeqs = str(round(math.ceil(solObjDict[solNum].negCumCov*negNumSeqs), 1))
		tp = int(float(posCumNumSeqs))
		fp = int(float(negCumNumSeqs))
		# tn = int(negCumNumSeqs)
		tn = negNumSeqs - fp
		fn = posNumSeqs - tp
		[se,sp,acc,fdr,ppv] = cal_metrics(tp,tn,fp,fn)
		new_jid = label + "_" + str(solNum)
		posCummCov = round(solObjDict[solNum].posCumCov*100, 1)
		negCummCov = round(solObjDict[solNum].negCumCov*100, 1)
		numMotifs = len(solObj.motifObjDict)
		latexFile.write('Motif solution number ' + str(solNum) + ' Num motifs:'+ str(numMotifs) + '\\newline PosCov:' + str(posCummCov) + ' NumSeqs:' + str(posCumNumSeqs) + '\n') 
		latexFile.write('\\newline NegCov:' + str(negCummCov) + ' NumSeqs:' + str(negCumNumSeqs) + '\n')
		makeLatexTable2(solObj.motifObjDict, logoDirName, latexFile, solNum, posNumSeqs, negNumSeqs)
		#fileList.append(latexFileName)
		# makeLatexTable_machine_learning_metrics
		# TP & TN & FP & FN & SE & SP & ACC & FDR & PPV & AUC
		# auc = run_roc_auc(new_jid,pos_seq,neg_seq,selectPwmFile)
		auc = run_roc_auc_fimo(new_jid,posNumSeqs,negNumSeqs,mName,dict_msp)
		makeLatexTable_machine_learning_metrics(latexFile,[tp,tn,fp,fn,se,sp,acc,fdr,ppv,auc])
		latexFile.write('\\begin{center} \n')
		latexFile.write('\includegraphics[scale=0.7]{' + new_jid + '.png}\n')
		
		latexFile.write('\end{center} \n')
		mast_out = jid +'_sol_' + str(solNum)+ "_mastVIS"
		run_mast_command = mast_command.replace("motif_file",selectPwmFile)
		run_mast_command = run_mast_command.replace("output_folder_mast",mast_out)
		os.system(run_mast_command)
		try:
			shutil.move(mast_out+"/mast.html", mast_out+".html")
			shutil.rmtree(mast_out)
			fileList.append(mast_out+".html")
			fileList.append(new_jid+".png")
		except:
			print "no motifs!"
			break
		
	
	latexFile.write('\end{document}')
	#close files
	latexFile.close()
	# visualization for all motifs in motif selection method
	selectAllPwmFile = jid + '_sol_all'  + '_select_motifs.pwm'
	makeSolPwmFile2(selectAllPwmFile, all_solution_motifs, allPwmFile)
	mast_out = jid +'_sol_all' + "_mastVIS"
	run_mast_command = mast_command.replace("motif_file",selectAllPwmFile)
	run_mast_command = run_mast_command.replace("output_folder_mast",mast_out)
	os.system(run_mast_command)
	try:
		shutil.move(mast_out+"/mast.html", mast_out+".html")
		shutil.rmtree(mast_out)
		fileList.append(mast_out+".html")
		fileList.append(selectAllPwmFile)
	except:
		print "no motifs!"
		
	
	result_pdf = jid + '_latexTables.pdf'
	pdflatex_command = 'pdflatex "'+latexFileName+'"'
	os.system(pdflatex_command)
	
	for outFile in PCA_plots_files:
		try:
			shutil.move(outFile, PCA_plots_folder)
		except:
			print outFile,"is not exist"
	
	#move files to results folder
	for outFile in fileList:
		shutil.move(outFile, output_folder)
	
	
	
	return output_folder,result_pdf
	
def Motif_discovery_accuracy_AUC(jid,confDict,label,csv_file):

	
	jid = jid + "_" + label
	print "Motif_discovery_accuracy_AUC jid : ",jid
	
	### 3-25-2016 added
	
	output_folder = jid + "_Motif_discovery_html_report"
	os.makedirs(output_folder)
	fileList = []
	
	# input files
	allPwmFile = confDict['input']['pwm_selected'] # for mast
	pos_seq = confDict['input']['pos_seq'] # for accuracy
	neg_seq = confDict['input']['neg_seq']# for accuracy
	
	posNumSeqs = findNumSeqs(pos_seq)
	negNumSeqs = findNumSeqs(neg_seq)
	
	
	# csv file
	csv_file = single_motif_rank(jid,confDict,"final_selection")
	fileList.append(csv_file)
	# images file
	motif_logo_folder = "motif_discovery_logos"
	meme_command = "meme2images -png " + allPwmFile + " " + motif_logo_folder
	os.system(meme_command)
	fileList.append(motif_logo_folder)
	# output motif report html file
	template = "%s/_templates/template.html"%(PATH)
	template_html = open(template).readlines()
	motif_logo_dict = copy_motif_logo(motif_logo_folder)
	rows = output_table(csv_file,motif_logo_dict)
	
	template_html[46] = "\n".join(rows)
	template_html[13] +=  str(posNumSeqs)
	template_html[16] +=  str(negNumSeqs)
	# line = template_html[19]
	# line = line.replace("SE_motif",str(se))
	# line = line.replace("SP_motif",str(sp))
	# line = line.replace("ACC_motif",str(acc))
	# line = line.replace("AUC_motif",str(auc))
	# template_html[19] = line
	
	html_file = open(jid+".html","wb")
	html_file.write("".join(template_html))
	html_file.close()
	fileList.append(jid+".html")
	
	for outFile in fileList:	
		print outFile
		shutil.move(outFile, output_folder)
		
	return output_folder

def output_row(list):
	out = "<tr>\n"
	ss = map(lambda x:"\t<td>"+x+"</td>",list)
	out += "\n".join(ss)
	out += "</tr>\n"
	return out


def output_table(csv_file,motif_logo_dict):
	rows = []
	for line in open(csv_file).readlines()[1:]:
		line = line.strip()
		line = line.split(",")
		motif_name = line[0]
		motif_name = motif_name.replace("_",".")
		path = motif_logo_dict[motif_name]
		line = ['<img src="'+path+'" height=40 width=100></img>'] + line
		rows.append(output_row(line))
	return rows	
	
	
def copy_motif_logo(folder_name):
	motif_logo = glob(folder_name+"/*.png")
	new_motif_logo_path = {}
	# copy_command = "cp f1i1l1e " + out_folder
	for m in motif_logo:
		motif_file_name = m.split("/")[-1]
		motif_name = motif_file_name[4:-4]
		motif_name = motif_name.replace('_',".")
		print "copy_motif_logo",motif_name
		# command = copy_command.replace("f1i1l1e",m)
		# os.system(command)
		new_motif_logo_path[motif_name] =  folder_name + "/"+motif_file_name
	return new_motif_logo_path	
	
	
def single_motif_rank(jid,confDict,label):
	#Motif obj dict between motif name and its object
	# motifObjDict = {}
	jid = jid + "_" + label
	outFileName = jid + "_SE_SP_ACC.csv"
	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	pos_mhit = confDict['input']['pos_mhit']
	neg_mhit = confDict['input']['neg_mhit']
	foreMotifDict = readHitFile(pos_mhit)
	backMotifDict = readHitFile(neg_mhit)
	# mdict = readPWM(allPwmFile)
	totalNumPosSeqs = findNumSeqs(pos_seq)
	totalNumNegSeqs = findNumSeqs(neg_seq)
	
	outFile = open(outFileName, 'wb')
	outFile.write('Motif,ForeCov,ForeNumSeqs,BackCov,BackNumSeqs,TP,FN,FP,TN,SN,SP,ACC\n')
	for motifName in foreMotifDict:
		posNumSeqs = len(foreMotifDict[motifName])
		posCov = round( 100*(posNumSeqs/totalNumPosSeqs), 1)
		if motifName in backMotifDict:
			negNumSeqs = len(backMotifDict[motifName])
			negCov = round( 100*(negNumSeqs/totalNumNegSeqs), 1)
		else:
			negNumSeqs = 0
			negCov= 0
		TP = int(posNumSeqs)
		FN = int(totalNumPosSeqs - TP)
		FP = int(negNumSeqs)
		TN = int(totalNumNegSeqs - FP)
		SN, SP, ACC = findAccMeas(TP, FN, FP, TN, totalNumPosSeqs, totalNumNegSeqs)
		
		lineList = [motifName, str(posCov), str(posNumSeqs), str(negCov), str(negNumSeqs), str(TP), str(FN), str(FP), str(TN), str(SN), str(SP), str(ACC)]
		lineStr = ','.join(lineList)
		outFile.write(lineStr + '\n')
		# motifObjDict[motifName] = MyMotif(TP, TN, FN, FP, SN, SP, ACC, float(posCov), float(negCov), float(posNumSeqs), float(negNumSeqs))
		
	outFile.close()
	
	# print 'number of motifs in dictionary:', len(motifObjDict)
	return outFileName

	
def single_motif_rank_unused(jid,confDict,label):
	jid = jid + "_" + label
	print "single_motif_rank jid : ",jid
	# output_folder = jid + "_single_motif_rank"
	# os.makedirs(output_folder)
	# fileList = [] 
	allPwmFile = confDict['input']['pwm_file']
	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	pos_mhit = confDict['input']['pos_mhit']
	neg_mhit = confDict['input']['neg_mhit']
	posMotifDict = readHitFile(pos_mhit)
	negMotifDict = readHitFile(neg_mhit)
	# mdict = readPWM(allPwmFile)
	posNumSeqs = findNumSeqs(pos_seq)
	negNumSeqs = findNumSeqs(neg_seq)
	##############
	# logoDirName = jid + '_sol_' + str(solNum) + '_logos'
	# makeMotifLogoFromPwm(selectPwmFile, logoDirName)
	# fileList.append(logoDirName)
	################
	table_file = jid + "_SE_SP_ACC.csv"
	fTable = open(table_file,"wb")
	# fileList.append(fTable)
	for mName in posMotifDict.keys():
		tp = len(posMotifDict[mName])
		try:
			fp = len(negMotifDict[mName])
		except:
			fp = 0
		tn = negNumSeqs - fp
		fn = posNumSeqs - tp
		[se,sp,acc,fdr,ppv] = 	cal_metrics(tp,tn,fp,fn)
		array = [mName,tp,fp,se,sp,acc]
		print >> fTable,",".join(map(lambda x:str(x),array))
	fTable.close()

	
	return table_file
	
def mhit2csv(jid,confDict,label):
	print "mhit2csv jid : ",jid
	# output_folder = jid + "_" + label + "_csv_file"
	# os.makedirs(output_folder)
	# fileList = [] 
	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	pos_mhit = confDict['input']['pos_mhit']
	neg_mhit = confDict['input']['neg_mhit']

	posMotifDict = readHitFile(pos_mhit)
	negMotifDict = readHitFile(neg_mhit)
	# assume posmhit file contains all the motifs, which is supposed to be true
	mName = list(posMotifDict.keys())

	pos_name = findName(pos_seq) # fasta file
	neg_name = findName(neg_seq)
	motifs_csv_file = jid + "_" + label + "_hit.csv"
	motifs_csv = open(motifs_csv_file,"wb")
	title = dp(mName)
	title.append("Class")
	print >>motifs_csv,",".join(title)
	for item in pos_name:				
		output_line = []			
		for motif in mName:			
			if item in posMotifDict[motif]:		
				output_line.append(1)	
			else:		
				output_line.append(0)	
		output_line.append("pos")			
		output_line = map(lambda x: str(x),output_line)			
		print >>motifs_csv,",".join(output_line)			
	for item in neg_name:				
		output_line = []			
		for motif in mName:			
			try:		
				if item in negMotifDict[motif]:	
					output_line.append(1)
				else:	
					output_line.append(0)
			except:		
				output_line.append(0)	
		output_line.append("neg")			
		output_line = map(lambda x: str(x),output_line)			
		print >>motifs_csv,",".join(output_line)			
	# fileList.append(motifs_csv)
	motifs_csv.close()
	
	return motifs_csv_file

def findName(fastaFile):
	handle = open(fastaFile, "rU")
	numSeqs = []
	for record in SeqIO.parse(handle, "fasta"):
		seqId = str(record.id)
		numSeqs.append(seqId)
	handle.close()
	
	return numSeqs
	
def cal_metrics(tp,tn,fp,fn):
	se = tp / (tp + fn)
	sp = tn / (fp + tn)
	acc = (tp + tn) / (tp + fn + fp + tn)
	fdr = fp / (tp + fp)
	ppv = tp / (tp + fp)
	return [se,sp,acc,fdr,ppv]
	
	
def parse_motif_selection_solution(solFile,solNumList):
	solObjDict = {}
	flag = 0
	with open(solFile, 'rb' ) as handler:
		for line in handler:
			line = line.strip()
			if re.search(r'#', line) or not line.strip():
				continue
			if re.search(r'Solution', line):
				split = line.split(':')
				num = int(split[1])
				if num in solNumList:
					flag = 1
					#make a list of motif objects between motif name and motif object; it is a motifObjDict = {}
					solObjDict[num] = Solution(num)
					print 'sol:',num
				else:
					flag = 0
				continue
			if flag == 1:
				split = line.split(',')
				motifName = split[0]
				posCov = float(split[1])*100
				posCov = round(posCov,1)
				negCov = float(split[3])*100
				negCov = round(negCov, 1)
				posCumCov = float(split[2])
				negCumCov = float(split[4])
				solObjDict[num].setPosCov(posCumCov)
				solObjDict[num].setNegCov(negCumCov)
				solObjDict[num].motifObjDict[motifName] = MotifObj(motifName, posCov, negCov, num)
	return solObjDict
# def parse_motif_selection_solution3(solFile,solNumList,pos_mhit):
def parse_motif_selection_solution_multicover(solFile,solNumList,pos_mhit):
	solObjDict = {}
	flag = 0
	num=1
	posMotifDict = readHitFile(pos_mhit)
	solObjDict[num] = Solution(num)
	with open(solFile, 'rb' ) as handler:
		for line in handler:
			line = line.strip()
			if re.search(r'#', line) or not line.strip():
				continue
			print "multicover################################"
			print line
			split = line.split(',')
			motifName = split[0]
			posCov = float(split[2])
			posCov = round(posCov,1)
			negCov = float(split[6])
			negCov = round(negCov, 1)
			posCumCov = float(split[4])
			negCumCov = float(split[8])
			
			solObjDict[num].setPosCov(posCumCov)
			solObjDict[num].setNegCov(negCumCov)
			solObjDict[num].motifObjDict[motifName] = MotifObj(motifName, posCov, negCov, num,posMotifDict[motifName])
	return solObjDict
def parse_motif_selection_solution_greedy(solFile,solNumList,pos_mhit):
	solObjDict = {}
	flag = 0
	posMotifDict = readHitFile(pos_mhit)
	with open(solFile, 'rb' ) as handler:
		for line in handler:
			line = line.strip()
			if re.search(r'#', line) or not line.strip():
				continue
			if re.search(r'Depth', line):
				split = line.split(':')
				num = int(split[1])
				if num in solNumList:
					flag = 1
					#make a list of motif objects between motif name and motif object; it is a motifObjDict = {}
					solObjDict[num] = Solution(num)
					print 'sol:',num
				else:
					flag = 0
				continue
			if flag == 1:
				split = line.split(',')
				motifName = split[0]
				posCov = float(split[2])
				posCov = round(posCov,1)
				negCov = float(split[6])
				negCov = round(negCov, 1)
				posCumCov = float(split[4]) / 100
				negCumCov = float(split[8]) / 100# doesn't make sense
				solObjDict[num].setPosCov(posCumCov)
				solObjDict[num].setNegCov(negCumCov)
				solObjDict[num].motifObjDict[motifName] = MotifObj(motifName, posCov, negCov, num,posMotifDict[motifName])
	return solObjDict
	
class MotifObj:
	"""Class to make an motif object
	
	"""
	def __init__(self):
		self.motifName = '' #name of the motif
		self.posCov = 0 
		self.negCov = 0 
		self.solNum = 0 #which solution does it belong to
		self.startPosList = []
		self.endPosList = []
		self.seqName = []
		
	def __init__(self,motifName, posCov, negCov, solNum,seq_array):
		self.motifName = motifName
		self.posCov = posCov
		self.negCov = negCov
		self.solNum = solNum
		self.seqName = seq_array
	
	def setVals_1(self, motifName, posCov, negCov, solNum):
		self.motifName = motifName
		self.posCov = posCov
		self.negCov = negCov
		self.solNum = solNum

	def setVals_2(self, motifName, posCov, negCov, solNum, startPos, endPos):
		self.motifName = motifName
		self.posCov = posCov
		self.negCov = negCov
		self.solNum = solNum
		self.startPosList.append(startPos)
		self.endPosList.append(endPos)
		
		
class Solution:
	def __init__(self, solNum):
		self.solNum = solNum
		self.motifObjDict = {} #dict between motif names and theri objects
		self.posCumCov = 0
		self.negCumCov = 0
		
	def setPosCov(self, inCov):
		self.posCumCov = inCov
		
	def setNegCov(self, inCov):
		self.negCumCov = inCov
	
def makeMotifLogoFromPwm(pwmFileName, outDirName):
	"""
	Given the pwm file in MEME format make the logos
	"""	
	command =  '/home/liyc/working/bin/meme2images -png ' + pwmFileName + ' ' + outDirName
	print 'meme2images command:', command
	try:
		os.system(command)
	except:
		print 'meme2images execution failed. Exiting'
		exit()


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


def readConfFile(confFile):
	""" Get the fasta file names from the config file""" 
	posFastaFile = ''
	negFastaFile = ''
	with open(confFile, 'rb') as handler:
		for line in handler:
			line = line.strip()
			if re.search(r'posFastaFile', line):
				split = line.split('=')
				posFastaFile = split[1]
				continue
			if re.search(r'negFastaFile', line):
				split = line.split('=')
				negFastaFile = split[1]
	
	return posFastaFile, negFastaFile

def makeLatexTable(motifObjDict, logoDirName, latexOut, solNum, totalPosNumSeqs, totalNegNumSeqs):
	"""Make a latex file for the table """
	#logoCAKGCAWT.png
	#header
	#latexOut = open(latexFileName, 'wb')
	lineList = [] 
	latexOut.write('\\begin{center} \n \\begin{tabular}{|c c c c c c|} \n \hline \n')
	latexOut.write('Logo & Motif & PosCov(\%) & PosNumSeq & NegCov(\%) & NegNumSeqs \\\ \n \hline \n')
	for key,value in sorted(motifObjDict.iteritems(), key = lambda (k,v):v.posCov, reverse=True):
		motifName = key
		motifObj = value
		lineList = []
		print motifName, motifObj.posCov
		motifLogo = '\includegraphics[scale=0.2]{' + 'logo'+motifName.replace('.','_')+'.png' +'}'
		logoPath = motifLogo
		if re.search(r'_', motifName):
			split = motifName.split('_')  
			motifName = '-'.join(split)
		posNumSeqs = round((motifObj.posCov/100)*totalPosNumSeqs)
		negNumSeqs = round((motifObj.negCov/100)*totalNegNumSeqs)
		lineList = [logoPath, motifName, str(motifObj.posCov), str(posNumSeqs), str(motifObj.negCov), str(negNumSeqs)]
		lineStr = ' & '.join(lineList)
		print lineStr
		latexOut.write(lineStr + '\\\ \n')
		
	
	latexOut.write('\hline \n \end{tabular} \n \end{center} \n')
	
	
	#close file
	#latexOut.close()
	
def makeLatexTable2(motifObjDict, logoDirName, latexOut, solNum, totalPosNumSeqs, totalNegNumSeqs):
	"""Make a latex file for the table """
	#logoCAKGCAWT.png
	#header
	#latexOut = open(latexFileName, 'wb')
	lineList = [] 
	latexOut.write('\\begin{center} \n \\begin{tabular}{|c c c c c c c c|} \n \hline \n')
	latexOut.write('Logo & Motif & PosCov(\%) & \#PosSeq & PosCumCov(\%) & \#PosCumSeq & NegCov(\%) & NegNumSeqs \\\ \n \hline \n')
	first_motif_flag = True
	Cum_seq_array = []
	Cum_seq_num = 0
	# Cum_seq_percent = 0.0
	
	for key,value in sorted(motifObjDict.iteritems(), key = lambda (k,v):v.posCov, reverse=True):
		motifName = key
		motifName = motifName.replace('.','_')
		motifName = motifName.replace('-','_')
		motifObj = value
		lineList = []
		print motifName, motifObj.posCov
		motifLogo = '\includegraphics[scale=0.2]{' + 'logo'+motifName+'.png' +'}'
		logoPath = motifLogo
		if re.search(r'_', motifName):
			split = motifName.split('_')  
			motifName = '-'.join(split)
		### 2-24-2016 added
		if len(motifName) >15:
			motifName = motifName[:4]+"-"+motifName[-10:]
		posNumSeqs = round((motifObj.posCov/100)*totalPosNumSeqs)
		negNumSeqs = round((motifObj.negCov/100)*totalNegNumSeqs)
		if first_motif_flag:
			length_new_add_seq = posNumSeqs
			percent_new_add_seq = motifObj.posCov
			first_motif_flag = False
			Cum_seq_array = motifObj.seqName
			Cum_seq_num = len(Cum_seq_array)
		else:
			pre_set =  set(Cum_seq_array)
			cur_set =  set(motifObj.seqName)
			new_add_seq = cur_set.difference(pre_set)
			length_new_add_seq = len(new_add_seq) + Cum_seq_num
			percent_new_add_seq = (length_new_add_seq/totalPosNumSeqs)*100
			Cum_seq_num = length_new_add_seq
			Cum_seq_array += motifObj.seqName
		lineList = [logoPath, motifName, str(motifObj.posCov), str(posNumSeqs),str(round(percent_new_add_seq,1)),
		str(length_new_add_seq),
		str(motifObj.negCov), str(negNumSeqs)]
		lineStr = ' & '.join(lineList)
		print lineStr
		latexOut.write(lineStr + '\\\ \n')
		# pre_motif_name = motifName
		# pre_motif_obj = motifObj
		# Cum_seq_num = len(motifObj.seqName)
		# for s in motifObj.seqName:
			# Cum_seq_array.append(s)
	
	latexOut.write('\hline \n \end{tabular} \n \end{center} \n')
	

def makeLatexTable_machine_learning_metrics(latexOut, lineList):
	"""Make a latex file for the table """
	#logoCAKGCAWT.png
	#header
	#latexOut = open(latexFileName, 'wb')
	latexOut.write('\\begin{center} \n \\begin{tabular}{|c c c c c c c c c c|} \n \hline \n')
	latexOut.write('TP & TN & FP & FN & SE & SP & ACC & FDR & PPV & AUC \\\ \n \hline \n')
	# lineStr = ' & '.join(lineList)
	lineStr = ' & '.join(map(str, lineList))
	print lineStr
	latexOut.write(lineStr + '\\\ \n')
	# http://www.cs.rpi.edu/~leen/misc-publications/SomeStatDefs.html
	latexOut.write('\hline \n \end{tabular} \n \end{center} \n')
	

def makeSolPwmFile2(pwmFileName, motifObjDict, allPwmFile):
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
				if motifName in motifObjDict:
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
	
def makeSolPwmFile(selectPwmFile, mName, mdict):
	pwmFile = open(pwmFileName, 'wb')
	pwmFile.write('MEME version 4.4\nALPHABET= ACGT\nstrands: + -\nBackground letter frequencies (from web form):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n')
	for name in mName:
		print >>pwmFile,"MOTIF",name
		print >>pwmFile,"letter-probability matrix:"
		print >>pwmFile,"\n".join(mdict[name].pwmLines)
		print >>pwmFile,"\n"
	pwmFile.close()
	
	
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

from decimal import Decimal	
def parseFIMO_forAUC(dict_msp,file,label):
	
	for line in open(file).readlines()[1:]:
		line = line.split()
		# print line
		m = line[0] # motif name
		s = line[1] + "_" + label # seq name
		p = Decimal(line[6]) # p value
		# print "m s p",m,s,p
		try:
			if dict_msp[m][s] > p:
				dict_msp[m][s] = p
		except:
			try:
				dict_msp[m][s] = p
			except:
				dict_msp[m] = {}
				dict_msp[m][s] = p
	# print dict.keys()
	# for k in dict.keys():
		# print dict[k].keys()
	
	
	
	
def findAccMeas(TP, FN, FP, TN, P, N):
	"""
	Find different measures
	Sn = TP / (TP + FN)
	Sp = TN / (FP + TN)
	Acc =  True positive + True negative / Total population OR TP + TN / P + N
	precision = Tp/Tp+FP
	"""
	sn = TP/(TP + FN)
	sp = TN / (FP + TN)
	acc = (TP+TN) / (P+N)
	
	return (sn*100), (100*sp), (100*acc)
	
class MyMotif:
	def __init__(self, tp, tn, fn, fp, sn, sp, acc, posCov, negCov, posNumSeqs, negNumSeqs):
		self.tp = tp
		self.tn = tn
		self.fn = fn
		self.fp = fp
		self.sn = sn
		self.sp = sp
		self.acc = acc
		self.posCov = posCov
		self.negCov = negCov
		self.posNumSeqs = posNumSeqs
		self.negNumSeqs = negNumSeqs
		
	
def run_Rami_style_multicover(jid,confDict,label):
	jid = jid + "_" + label
	print "run greedy_multicover_output jid : ",jid
	
	### 2-24-2016 added
	PCA_plots_files = []
	PCA_plots_folder = jid + "_final_selection_PCA"
	os.makedirs(PCA_plots_folder)
	
	output_folder = jid + "_rami_style_output"
	os.makedirs(output_folder)
	fileList = [] 
	fileList.append(PCA_plots_folder)
	allPwmFile = confDict['input']['pwm_file']
	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	solFile = confDict['input']['greedy_multicover_sol']
	# num_sols = confDict['input']['num_sols']
	pos_mhit = confDict['input']['pos_mhit']
	neg_mhit = confDict['input']['neg_mhit']
	pos_fimo = confDict['input']['pos_fimo']
	neg_fimo = confDict['input']['neg_fimo']
	dict_msp = {}
	parseFIMO_forAUC(dict_msp,pos_fimo,"positive")
	parseFIMO_forAUC(dict_msp,neg_fimo,"negative")
	all_solution_motifs = {}
	# mdict = readPWM(allPwmFile)
	solNumList = range(1,int(1)+1)
	mast_command = "mast motif_file " + pos_seq + " -o output_folder_mast -ev 9999 -notext -nostatus"
	solObjDict = parse_motif_selection_solution_multicover(solFile,solNumList,pos_mhit)
	posNumSeqs = findNumSeqs(pos_seq)
	negNumSeqs = findNumSeqs(neg_seq)
	
	latexFileName = jid + '_latexTables.tex'
	fileList.append(latexFileName)
	latexFile = open(latexFileName, 'wb')
	header=r'''\documentclass{article}
	\usepackage{graphicx} 
	\usepackage{geometry}
	\geometry{margin=0in}
	'''
	latexFile.write(header)
	latexFile.write('\\begin{document}\n')
	#check it
	for solNum, solObj in solObjDict.items():
		print 'S:',solNum,solObj.solNum
		for motifName in solObj.motifObjDict:
			print '\t', motifName
			all_solution_motifs[motifName] = solObj.motifObjDict[motifName]
		#make a pwm file for this solution
		#out file for selected PWMs
		### Yichao Edit here
		mName = solObj.motifObjDict.keys()
		if len(mName) < 1:
			continue
		# print mName
		# makeSolPwmFile(selectPwmFile, mName, mdict)
		### 2-26-2016 added 
		out_pos_hit = jid + "_final_selection_pos.mhit"
		out_neg_hit = jid + "_final_selection_neg.mhit"
		hit2hit(pos_mhit,neg_mhit,out_pos_hit,out_neg_hit,mName)
		confDict['input']['pos_mhit'] = out_pos_hit
		confDict['input']['neg_mhit'] = out_neg_hit
		file = single_motif_rank(jid,confDict,"final_selection")
		PCA_plots_files.append(out_pos_hit)
		PCA_plots_files.append(out_neg_hit)
		PCA_plots_files.append(file)
		file = mhit2csv(jid,confDict,"_final_selection_")
		PCA_plots_files.append(file)
		output_file = jid +"_"+"final_selection"+"_PCAplot.png"
		R = "Rscript %s/pca_plot.R "%(PATH)+ file +" "+ output_file
		os.system(R)
		PCA_plots_files.append(output_file)
		
		### finish ###
		selectPwmFile = jid + '_select_motifs.pwm'
		
		makeSolPwmFile2(selectPwmFile, solObj.motifObjDict, allPwmFile)
		 
		fileList.append(selectPwmFile)
		
		#make the logos
		#make motif logos from the selected PWMs
		
		# ----------------------------------------------
		
		logoDirName = jid + '_logos'
		makeMotifLogoFromPwm(selectPwmFile, logoDirName)
		fileList.append(logoDirName)
		
		#make the latex file
		latexFile.write('\graphicspath{ {./'+ logoDirName +'/} }\n')
		latexFile.write('\clearpage \n')
		# print 'posCum:', solObjDict[solNum].posCumCov
		posCumNumSeqs = str(round(math.ceil(solObjDict[solNum].posCumCov*posNumSeqs/100), 1))
		# negCumNumSeqs = str(round(math.ceil(solObjDict[solNum].negCumCov*posNumSeqs), 1))
		# print 'posCumNumSeqs:', posCumNumSeqs
		
		negCumNumSeqs = str(round(math.ceil(solObjDict[solNum].negCumCov*negNumSeqs/100), 1))
		tp = int(float(posCumNumSeqs))
		fp = int(float(negCumNumSeqs))
		# tn = int(negCumNumSeqs)
		tn = negNumSeqs - fp
		fn = posNumSeqs - tp
		[se,sp,acc,fdr,ppv] = cal_metrics(tp,tn,fp,fn)
		new_jid = label + "_" + str(solNum)
		posCummCov = round(solObjDict[solNum].posCumCov, 1)
		negCummCov = round(solObjDict[solNum].negCumCov, 1)
		numMotifs = len(solObj.motifObjDict)
		latexFile.write('Motif solution number ' + str(solNum) + ' Num motifs:'+ str(numMotifs) + '\\newline PosCov:' + str(posCummCov) + ' NumSeqs:' + str(posCumNumSeqs) + '\n') 
		latexFile.write('\\newline NegCov:' + str(negCummCov) + ' NumSeqs:' + str(negCumNumSeqs) + '\n')
		makeLatexTable2(solObj.motifObjDict, logoDirName, latexFile, solNum, posNumSeqs, negNumSeqs)
		#fileList.append(latexFileName)
		# makeLatexTable_machine_learning_metrics
		# TP & TN & FP & FN & SE & SP & ACC & FDR & PPV & AUC
		# auc = run_roc_auc(new_jid,pos_seq,neg_seq,selectPwmFile)
		auc = run_roc_auc_fimo(new_jid,posNumSeqs,negNumSeqs,mName,dict_msp)
		makeLatexTable_machine_learning_metrics(latexFile,[tp,tn,fp,fn,se,sp,acc,fdr,ppv,auc])
		latexFile.write('\\begin{center} \n')
		latexFile.write('\includegraphics[scale=0.7]{' + new_jid + '.png}\n')
		
		latexFile.write('\end{center} \n')
		mast_out = jid +'_sol_' + str(solNum)+ "_mastVIS"
		run_mast_command = mast_command.replace("motif_file",selectPwmFile)
		run_mast_command = run_mast_command.replace("output_folder_mast",mast_out)
		os.system(run_mast_command)
		try:
			shutil.move(mast_out+"/mast.html", mast_out+".html")
			shutil.rmtree(mast_out)
			fileList.append(mast_out+".html")
			fileList.append(new_jid+".png")
		except:
			print "no motifs"
		
	
	latexFile.write('\end{document}')
	#close files
	latexFile.close()
	# visualization for all motifs in motif selection method
	selectAllPwmFile = jid + '_sol_all'  + '_select_motifs.pwm'
	makeSolPwmFile2(selectAllPwmFile, all_solution_motifs, allPwmFile)
	mast_out = jid +'_sol_all' + "_mastVIS"
	run_mast_command = mast_command.replace("motif_file",selectAllPwmFile)
	run_mast_command = run_mast_command.replace("output_folder_mast",mast_out)
	os.system(run_mast_command)
	try:
		shutil.move(mast_out+"/mast.html", mast_out+".html")
		shutil.rmtree(mast_out)
		fileList.append(mast_out+".html")
		fileList.append(selectAllPwmFile)
	except:
		print "no motifs"
	
	result_pdf = jid + '_latexTables.pdf'
	pdflatex_command = 'pdflatex "'+latexFileName+'"'
	os.system(pdflatex_command)
	
	
	for outFile in PCA_plots_files:
		try:
			shutil.move(outFile, PCA_plots_folder)
		except:
			print outFile,"is not exist"
	
	#move files to results folder
	for outFile in fileList:
		shutil.move(outFile, output_folder)
	
	
	
	return output_folder,result_pdf
	
def parse_motif_selection_solution2(solFile,solNumList,pos_mhit):
	solObjDict = {}
	flag = 0
	posMotifDict = readHitFile(pos_mhit)
	with open(solFile, 'rb' ) as handler:
		for line in handler:
			line = line.strip()
			if re.search(r'#', line) or not line.strip():
				continue
			if re.search(r'Solution', line):
				split = line.split(':')
				num = int(split[1])
				if num in solNumList:
					flag = 1
					#make a list of motif objects between motif name and motif object; it is a motifObjDict = {}
					solObjDict[num] = Solution(num)
					print 'sol:',num
				else:
					flag = 0
				continue
			if flag == 1:
				split = line.split(',')
				motifName = split[0]
				posCov = float(split[1])*100
				posCov = round(posCov,1)
				negCov = float(split[3])*100
				negCov = round(negCov, 1)
				posCumCov = float(split[2])
				negCumCov = float(split[4])
				solObjDict[num].setPosCov(posCumCov)
				solObjDict[num].setNegCov(negCumCov)
				solObjDict[num].motifObjDict[motifName] = MotifObj(motifName, posCov, negCov, num,posMotifDict[motifName])
	return solObjDict

	
def mhit2csv_with_label(jid,confDict,label):
	print "mhit2csv jid : ",jid
	# output_folder = jid + "_" + label + "_csv_file"
	# os.makedirs(output_folder)
	# fileList = [] 
	pos_seq = confDict['input']['pos_seq']
	pos_mhit = confDict['input']['pos_mhit']

	posMotifDict = readHitFile(pos_mhit)
	# assume posmhit file contains all the motifs, which is supposed to be true
	mName = list(posMotifDict.keys())

	pos_name = findName(pos_seq) # fasta file
	motifs_csv_file = jid + "_" + label + "_for_clustering.csv"
	motifs_csv = open(motifs_csv_file,"wb")
	title = dp(mName)
	title.append("Class")
	print >>motifs_csv,",".join(title)
	for item in pos_name:				
		output_line = []			
		for motif in mName:			
			if item in posMotifDict[motif]:		
				output_line.append(1)	
			else:		
				output_line.append(0)	
		output_line.append(item)			
		output_line = map(lambda x: str(x),output_line)			
		print >>motifs_csv,",".join(output_line)			
	motifs_csv.close()
	
	return motifs_csv_file
