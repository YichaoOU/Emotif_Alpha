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
from motif_filtering import *
from motif_output import *


def run_motif_selection(jid, confDict):
	fileList = []
	# filtering step

	if confDict['job_type']['motif_filtering'] == 'true':			
				
				
		all_motifs_file_list = []		
		label = "all_motifs"		
		label_folder = jid + "_"+label		
		os.makedirs(label_folder)		
		# output all motifs SE SP ACC measurement table		
		file = single_motif_rank(jid,confDict,label) # right now this returns a folder?		
		all_motifs_file_list.append(file)		
				
		# output motif-seq csv hit table		
		file = mhit2csv(jid,confDict,label)		
		all_motifs_file_list.append(file)		
				
		# output PCA plot		
		output_file = jid +"_"+label+"_PCAplot.png"		
		R = "Rscript /usr/local/lib/python2.7/dist-packages/Emotif/pca_plot.py "+ file +" "+ output_file		
		os.system(R)		
		all_motifs_file_list.append(output_file)		
				
				
				
		arff_file = jid +"_"+label+".arff"		
		python_script = "python /usr/local/lib/python2.7/dist-packages/Emotif/csv2arff.py " + file + " " + arff_file		
		os.system(python_script)		
		fileList.append(arff_file)		
		for outFile in all_motifs_file_list:		
			print outFile
			try:
				shutil.move(outFile, label_folder)
			except:
				print outFile,"is not exist!"
		fileList.append(label_folder)		
		final_motifs = []		
		MF_folder = []		
		MF_dir = jid+"_Motif_Filtering"		
		os.makedirs(MF_dir)		
		# if confDict['motif_filtering']['coverage_filter'] == 'true':		
			# file,poshit,neghit = prefilter(jid,confDict)	
			# fileList.append(file)	
			# confDict['input']['pos_mhit'] = file + "/" + poshit	
			# confDict['input']['neg_mhit'] = file + "/" + neghit	
		if confDict['motif_filtering']['rf_gini_filter'] == 'true':		
			top_gini_motifs = RF_gini_filter(jid,confDict,arff_file,MF_folder)	
			# fileList.append(file)	
			final_motifs += top_gini_motifs	
		if confDict['motif_filtering']['rf_entropy_filter'] == 'true':		
			top_entropy_motifs = RF_entropy_filter(jid,confDict,arff_file,MF_folder)	
			# fileList.append(file)	
			final_motifs += top_entropy_motifs	
		# not impletemented yet		
		if confDict['motif_filtering']['svm_filter'] == 'true':		
			top_motifs = SVM_filter(jid,confDict,arff_file)	
			# fileList.append(file)	
				
		for outFile in MF_folder:		
			print outFile	
			shutil.move(outFile, MF_dir)			
		fileList.append(MF_dir)		
		# get final motif PWM files		
		final_motifs = list(set(final_motifs))
		if len(final_motifs) == 0:
			return fileList
		allPwmFile = confDict['input']['pwm_file']		
		final_PWM_file = jid + "_final_motifs.pwm"		
		subPwmFile(final_PWM_file,final_motifs,allPwmFile)		
		confDict['input']['pwm_file'] = final_PWM_file		
		fileList.append(final_PWM_file)		
		# get final mhit files		
		in_pos_hit = confDict['input']['pos_mhit']		
		in_neg_hit = confDict['input']['neg_mhit']		
		out_pos_hit = jid + "_final_pos.mhit"		
		out_neg_hit = jid + "_final_neg.mhit"		
		hit2hit(in_pos_hit,in_neg_hit,out_pos_hit,out_neg_hit,final_motifs)		
		fileList.append(out_pos_hit)		
		fileList.append(out_neg_hit)		
		confDict['input']['pos_mhit'] = out_pos_hit		
		confDict['input']['neg_mhit'] = out_neg_hit		
		confDict['input']['before_pos_mhit'] = in_pos_hit		
		confDict['input']['before_neg_mhit'] = in_neg_hit		
		# output final motifs SE SP ACC measurement table		
		final_motifs_file_list = []		
		label = "filtered_motifs"		
		label_folder = jid + "_"+label		
		os.makedirs(label_folder)		
		csv_file = single_motif_rank(jid,confDict,label)		
		final_motifs_file_list.append(csv_file)		
		# output final motif-seq csv hit table		
		file = mhit2csv(jid,confDict,label)		
		final_motifs_file_list.append(file)		
				
		# output final motifs PCA plot		
		output_file = jid +"_"+label+"_PCAplot.png"		
		R = "Rscript /usr/local/lib/python2.7/dist-packages/Emotif/pca_plot.py "+ file +" "+ output_file		
		os.system(R)		
		final_motifs_file_list.append(output_file)		
		
		# html output
		confDict['input']['pwm_selected'] = final_PWM_file
		file = Motif_discovery_accuracy_AUC(jid,confDict,"all_motif_ranking",csv_file)
		final_motifs_file_list.append(file)
		for outFile in final_motifs_file_list:		
			print outFile	
			try:
				shutil.move(outFile, label_folder)
			except:
				print outFile,"is not exist!"
		fileList.append(label_folder)	
		
	if confDict['job_type']['motif_selection'] == 'true':
		if confDict['motif_selection']['genetic_algo'] == 'true':
			file = run_GA(jid,confDict)
			fileList.append(file)
			confDict['input']['num_sols'] = get_GA_num_sol(file)
			confDict['input']['ga_solution'] = get_GA_sol_file(file)
		if confDict['motif_selection']['genetic_algo_multi'] == 'true':
			file = run_GA_multiCover(jid,confDict)
			fileList.append(file)
			confDict['input']['num_sols'] = get_GA_num_sol(file)
			confDict['input']['ga_solution'] = get_GA_sol_file(file)
		file,temp = run_others(jid,confDict)
		greedy_result_file = get_Greedy_sol_file(file,jid)
		confDict['input']['greedy_multicover_sol'] = file+"/"+temp
		confDict['input']['greedy_result_file'] = greedy_result_file
		
		
		fileList.append(file)
		
		
	return fileList


def run_GA(jid, confDict):

	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	output_folder = jid + '_motif_selection_GA'
	pos_mhit = confDict['input']['pos_mhit']
	neg_mhit = confDict['input']['neg_mhit']
	GA_path = confDict['genetic_algo']['path']
	
	niter = confDict['genetic_algo']['numevals']
	maxNeg = confDict['genetic_algo']['maxnegcovpermotif']
	filter = confDict['genetic_algo']['filtercutoff']
	penalty = confDict['genetic_algo']['penaltyvalue']
	
	GA_confFile = jid + '_GA.conf'
	f=open(GA_confFile,'w')
	print >>f,"posFastaFile="+pos_seq
	print >>f,"negFastaFile="+neg_seq
	print >>f,"posMotifHitFile="+pos_mhit
	print >>f,"negMotifHitFile="+neg_mhit
	print >>f,"numEvals="+niter
	print >>f,"maxNegCovPerMotif="+maxNeg
	print >>f,"filterCutoff="+filter
	print >>f,"penaltyValue="+penalty
	
	f.close()
	
	
	GA_command = GA_path + ' ' + output_folder + ' ' + GA_confFile
	print GA_command
	os.system(GA_command)
	os.system('mv '+ GA_confFile + ' ' + output_folder)
	return output_folder
	
def run_GA_multiCover(jid, confDict):

	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	output_folder = jid + '_motif_selection_GA_multiCover'
	pos_mhit = confDict['input']['pos_mhit']
	neg_mhit = confDict['input']['neg_mhit']
	GA_path = confDict['genetic_algo_multi']['path']
	
	niter = confDict['genetic_algo_multi']['numevals']
	maxNeg = confDict['genetic_algo_multi']['maxnegcovpermotif']
	minPos = confDict['genetic_algo_multi']['minposcov_threshold']
	filter = confDict['genetic_algo_multi']['filtercutoff']
	penalty = confDict['genetic_algo_multi']['penaltyvalue']
	
	GA_confFile = jid + '_GA_multiCover.conf'
	f=open(GA_confFile,'w')
	print >>f,"posFastaFile="+pos_seq
	print >>f,"negFastaFile="+neg_seq
	print >>f,"posMotifHitFile="+pos_mhit
	print >>f,"negMotifHitFile="+neg_mhit
	print >>f,"numEvals="+niter
	print >>f,"maxNegCovPerMotif="+maxNeg
	print >>f,"filterCutoff="+filter
	print >>f,"penaltyValue="+penalty
	print >>f,"minPosCov_threshold="+minPos
	
	f.close()
	
	
	GA_command = GA_path + ' ' + output_folder + ' ' + GA_confFile
	print GA_command
	os.system(GA_command)
	os.system('mv '+ GA_confFile + ' ' + output_folder)
	return output_folder
	
def run_SVM_RFE(jid, confDict):

	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	output_folder = jid + '_motif_selection_GA'
	pos_mhit = confDict['input']['pos_mhit']
	neg_mhit = confDict['input']['neg_mhit']
	GA_path = confDict['genetic_algo']['path']
	
	niter = confDict['genetic_algo']['numevals']
	maxNeg = confDict['genetic_algo']['maxnegcovpermotif']
	filter = confDict['genetic_algo']['filtercutoff']
	penalty = confDict['genetic_algo']['penaltyvalue']
	
	GA_confFile = jid + '_GA.conf'
	f=open(GA_confFile,'w')
	print >>f,"posFastaFile="+pos_seq
	print >>f,"negFastaFile="+neg_seq
	print >>f,"posMotifHitFile="+pos_mhit
	print >>f,"negMotifHitFile="+neg_mhit
	print >>f,"numEvals="+niter
	print >>f,"maxNegCovPerMotif="+maxNeg
	print >>f,"filterCutoff="+filter
	print >>f,"penaltyValue="+penalty
	
	f.close()
	
	
	GA_command = GA_path + ' ' + output_folder + ' ' + GA_confFile
	print GA_command
	os.system(GA_command)
	os.system('mv '+ GA_confFile + ' ' + output_folder)
	return output_folder



def run_others(jid, confDict):

	pos_seq = confDict['input']['pos_seq']
	# may need to change later for other algorithms
	output_folder = jid + '_motif_selection_greedy'
	pos_mhit = confDict['input']['pos_mhit']
	depth = int(confDict['motif_selection']['depth'])
	filterThresh = float(confDict['motif_selection']['filter_percentage'])
	epsilon = confDict['tomtom']['evalue']
	
	neg_seq = confDict['input']['neg_seq']
	neg_mhit = confDict['input']['neg_mhit']
	
	fileList = []

	fimoDict = read_mhit(pos_mhit)
	backFimoDict = read_mhit(neg_mhit)
	motifIdDict, idMotifDict, seqIdDict, idSeqDict, Uset, Sdict = general_utils.processMotifHitFile_1(pos_mhit)
	foreNumSeqs = general_utils.findNumSeqs(pos_seq)
	backNumSeqs = general_utils.findNumSeqs(neg_seq)
	filterMinNumSeq = filterThresh * foreNumSeqs
	
	### 2-24-2016 added
	# greedy_multicover_selected_motifs = []	
	depthOutFileName=""
	
	result_dict = {} 
	# result_dict[ilp] = dict_obj
	
	# ------------------------------------------
	#check if motif refinement requested to remove redundant motifs when applying the seq_cov algthms
	#the tomtomDict will be used with the seq coverage algorithm. The dict is between all the motifs to check their similarity 
	tomtomDict = {}
	if confDict['motif_selection']['tomtom'] == 'true':
		print 'Peform motif refine using tomtom'
		#dict that stores parameters for the tomtom tool
		tomParaDict = {}
		tomParaDict['evalue'] = epsilon
		#tomtom output directory name
		tomtomOutDir = jid + '_Tomtom_refine'
		#we are comparing the motifs against each other, so use the same PWM file
		Tomtom.callTomtom(confDict['input']['pwm_file'], confDict['input']['pwm_file'], tomParaDict, tomtomOutDir)
		fileList.append(tomtomOutDir)
		#parse the tomtom output file and return a dict between motif Name and list[] of motifs that match it
		tomtomDict = Tomtom.parse(tomtomOutDir+'/tomtom.txt')
	# ------------------------------------------
	if confDict['motif_selection']['kmer'] == 'true':
		motif_names = list(motifIdDict.keys())
		for m in motif_names:
			tomtomDict[m]=[]
		
	
	#####
	#ILP 
	####
	if confDict['motif_selection']['ilp'] == 'true':
		print 'ILP algorithm'
		depthDict = ILP.callILPDepth(jid, fimoDict, tomtomDict, depth, fileList)
		result_dict["ilp"] = depthDict
		
	####	
	#branch and cut	
	####
	if confDict['motif_selection']['branch_cut'] == 'true':
		print 'Branch and cut algorithm'
		depthDict = branch_cut.callBranchDepth(jid, fimoDict, tomtomDict, depth, fileList)
		result_dict["branch_cut"] = depthDict
		
	
	####	
	#required cover	
	####
	if confDict['motif_selection']['required_cover'] == 'true':
		print 'Required cover algorithm'
		depthDict = required_cover.callRequiredDepth(jid, fimoDict, tomtomDict, depth, fileList)
		result_dict["required_cover"] = depthDict		
	
	####
	#greedy and/or refinement	
	####
	if confDict['motif_selection']['greedy'] == 'true':
		print 'Perform greedy sequence coverage with refinement and depth:', depth
		depthDict = greedy.callGreedyDepthRefine(Uset, Sdict, depth, tomtomDict, motifIdDict, idMotifDict, filterMinNumSeq)
		result_dict["greedy_refine"] = depthDict
		# Greedy without refinement is not working properly
		# print 'Perform greedy sequence coverage with no refinement and depth:', depth
		# depthDict = greedy.callGreedyDepth(Uset, Sdict, depth)
		# result_dict["greedy_NO_refine"] = depthDict
		# final_motif_dict['greedy_multicover'] = depthDict
	for key in result_dict.keys():
		depthOutFileName_main = jid + '_depth_results_' + key
		depthDict = result_dict[key]
		for depthVal in sorted(depthDict.iterkeys(), key=int):
			depthOutFileName = depthOutFileName_main + "_" + str(depthVal) + '.csv'
			finalSelectMotifList = outputResults(depthOutFileName, depthDict, backFimoDict, filterMinNumSeq, idMotifDict,foreNumSeqs, backNumSeqs)
			print 'final:', finalSelectMotifList
			fileList.append(depthOutFileName)
			
	####
	#greedy MultiCover Algorithm Rami Dec-2015	
	####
	if confDict['motif_selection']['greedy_multicover'] == 'true':
		print 'Perform greedy multi-sequence coverage with refinement and depth:', depth
		# depthDict = greedy.callGreedy_MultiCover(Uset, Sdict, depth, tomtomDict, motifIdDict, idMotifDict, filterMinNumSeq)
		depthOutFileName = jid + '_greedy_multicovery_' + str(depth) + "_result.csv"
		selectMotifList = greedy.callGreedy_MultiCover(Uset, Sdict, depth, tomtomDict, motifIdDict, idMotifDict, filterMinNumSeq, seqIdDict, idSeqDict, depthOutFileName, fimoDict, backFimoDict, foreNumSeqs, backNumSeqs)
		# greedy_multicover_selected_motifs = selectMotifList
		# result_dict["greedy_refine"] = depthDict
		# Greedy without refinement is not working properly
		# print 'Perform greedy sequence coverage with no refinement and depth:', depth
		# depthDict = greedy.callGreedyDepth(Uset, Sdict, depth)
		# result_dict["greedy_NO_refine"] = depthDict
		fileList.append(depthOutFileName)	
			
	for outFile in fileList:
		print outFile,output_folder
		shutil.move(outFile, output_folder)
	return output_folder,depthOutFileName
	
def get_GA_num_sol(folder):
	file = folder + "/" + folder + "_motif_num_sols.csv"
	line = open(file).readlines()[0]
	num = line.split(":")[-1]
	return num
	
	
def get_GA_sol_file(folder):
	file = folder + "/" + folder + "_results.csv"
	return file
	
def get_GA_multiCover_num_sol(folder):
	file = folder + "/" + folder + "_motif_num_sols.csv"
	line = open(file).readlines()[0]
	num = line.split(":")[-1]
	return num
	
	
def get_GA_multiCover_sol_file(folder):
	file = folder + "/" + folder + "_results.csv"
	return file
	
def get_Greedy_sol_file(folder,jid):
	file = folder + "/" + jid + "_depth_results_greedy_refine_1.csv"
	return file

	
def prefilter(jid, confDict):

	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	output_folder = jid + '_prefilter'
	pos_mhit = confDict['input']['pos_mhit']
	neg_mhit = confDict['input']['neg_mhit']
	maxNegCov = float(confDict['prefilter']['maxnegcov'])
	minPosCov = float(confDict['prefilter']['minposcov'])
	
	
	os.makedirs(output_folder)
	fileList = []
	
	totalNumPosSeqs = general_utils.findNumSeqs(pos_seq)
	totalNumNegSeqs = general_utils.findNumSeqs(neg_seq)
	print 'numPosSeqs:', totalNumPosSeqs,'numNegSeqs:', totalNumNegSeqs
	#dict between motif names and their seq lists
	posMotifDict = readHitFile(pos_mhit)
	negMotifDict = readHitFile(neg_mhit)
	# for k in negMotifDict.keys():
		# print "Neg Check, ",k,len(negMotifDict[k])
	
	print 'maxNegCov:', maxNegCov,'minPosCov:', minPosCov
	filterPosHitFile = jid + '_pos_filter_hitFile_maxNegCov_' + str(int(maxNegCov*100)) + '.mhit'
	filterNegHitFile = jid + '_neg_filter_hitFile_maxNegCov_' + str(int(maxNegCov*100)) + '.mhit'
	fileList.append(filterPosHitFile)
	fileList.append(filterNegHitFile)
	filtMotifList = filterHits(posMotifDict, negMotifDict, maxNegCov, minPosCov, totalNumPosSeqs, totalNumNegSeqs)
	
	#write hit files
	writeFiltFiles(filtMotifList, filterPosHitFile, filterNegHitFile, posMotifDict, negMotifDict)
	
	#write the pos and neg cov of the filtered motifs
	covInfoFileName = jid + '_cov_info.csv'
	fileList.append(covInfoFileName)
	writeCovFiltFiles(filtMotifList, posMotifDict, negMotifDict, totalNumPosSeqs, totalNumNegSeqs, covInfoFileName)
	
	#posFiltFileName = jid + '_pos_filtered_hits_' + str(minNegCov)
	#fileList.append(posFiltFileName)
	#negFiltFileName = jid + '_neg_filtered_hits_' + str(minNegCov)
	#fileList.append(negFiltFileName)
	
	
	#move files to results folder
	for outFile in fileList:
		shutil.move(outFile, output_folder)
	
	return output_folder,filterPosHitFile,filterNegHitFile
	
	
	
	
	
	
	
	
	
	
	
	
	

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

def output_filtering_diagram(jid,twoD_array,t):
	import numpy as np
	import matplotlib.pylab as plt
	labels = ["0%","5%","10%","15%","20%","25%","30%","35%"]
	X_name = "Max Neg Cov"
	Y_name = "Min Pos Cov"
	matrix = np.matrix(twoD_array)
	output_csv = jid+"_filtering_metrics.csv"
	output_svg = jid+"_filtering_metrics.svg"
	np.savetxt(output_csv, matrix, delimiter=",")
	
	
def filterHits(posMotifDict, negMotifDict, maxNegCov, minPosCov, totalNumPosSeqs, totalNumNegSeqs):
	"""
	maxNegCov: the maximum coverage of the motif
	minPosCov: the minimum positive coverage of the motif
	"""
	##dict between filtered motifs and their sequences
	metrics_array = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35]
	filtMotifList = []
	for motifName, seqList in posMotifDict.items():
		 
		posNumSeqs = len(seqList)
		# print motifName,posNumSeqs
		posCov = posNumSeqs/totalNumPosSeqs
		# print "posCov",posCov
		if posCov < minPosCov:
			# print "if posCov < minPosCov:"
			continue 
		if motifName in negMotifDict:
			negNumSeqs = len(negMotifDict[motifName])
			# print "Negative: ",motifName,negNumSeqs
			negCov = negNumSeqs/totalNumNegSeqs
		else:
			negCov = 0
		if negCov > maxNegCov:
			continue
		filtMotifList.append(motifName)
	
	print 'filtMotfList:', len(filtMotifList)
	##dict between filtered motifs and their sequences
	#filtMotifList = []
	#c = 0
	#if maxNegCov == 0:
		##for the zero neg case
		##posList = []
		##negList = []
		##d = 0
		#for motifName, seqList in posMotifDict.items():
			#if motifName not in negMotifDict:
				#c+=1
				#filtMotifList.append(motifName)
				##d+=1
				##for seqName in seqList:
					##if seqName not in posList:
						##posList.append(seqName)
		
		##print 'len Pos:', len(posList), 'd:',d
		##loop thru the negative motifs and do the filtering
		#print 'c:', c,'filtMotfList:', len(filtMotifList)
		#return filtMotifList
	
	
	
	#for motifName, seqList in negMotifDict.items():
		#negNumSeqs = len(seqList)
		#negPer = negNumSeqs/totalNumNegSeqs
		#if negPer < maxNegCov:
			#c+=1
			##print 'motif:', motifName,'negNumSeqs:', negNumSeqs,'negPer:', negPer
			#filtMotifList.append(motifName)
	#print 'c:', c,'filtMotfList:', len(filtMotifList)
	
	negSeqList = []
	posSeqList = []
	for motifName in filtMotifList:
		#print 'motif:',motifName
		for seqName in posMotifDict[motifName]:
			if seqName not in posSeqList:
				posSeqList.append(seqName)
		if motifName in negMotifDict:
			for seqName in negMotifDict[motifName]:
				if seqName not in negSeqList:
					negSeqList.append(seqName)
					
		#break
		
	print 'negPer:', maxNegCov, 'posSeqList len:',len(posSeqList),'posPer:', round(100*len(posSeqList)/totalNumPosSeqs,1),'negSeqList len:', len(negSeqList),'Per:',round(100*len(negSeqList)/totalNumNegSeqs,1)

	return filtMotifList

	#
	
	
def writeFiltFiles(filtMotifList, filterPosHitFile, filterNegHitFile, posMotifDict, negMotifDict):
	posFile = open(filterPosHitFile, 'wb')
	negFile = open(filterNegHitFile, 'wb')
	
	d = 0
	for motifName in filtMotifList:
		# print "in motif_selection, motif name:",motifName
		d+=1
		posFile.write('>' + motifName + '\n')
		for seqName in posMotifDict[motifName]:
			posFile.write(seqName + '\n')
		if motifName in negMotifDict:
			negFile.write('>' + motifName + '\n')
			for seqName in negMotifDict[motifName]:
				negFile.write(seqName + '\n')
	print 'd:',d
	#close files
	posFile.close()
	negFile.close()


def writeCovFiltFiles(filtMotifList, posMotifDict, negMotifDict, totalNumPosSeqs, totalNumNegSeqs, covInfoFileName):
	outFile = open(covInfoFileName, 'wb')
	outFile.write('Motif,PosCov,PosNumSeqs,NegCov,NegNumSeqs\n')
	for motifName in filtMotifList:
		posNumSeqs = len(posMotifDict[motifName])
		posCov = round( 100*(posNumSeqs/totalNumPosSeqs), 1)
		if motifName in negMotifDict:
			negNumSeqs = len(negMotifDict[motifName])
			negCov = round( 100*(negNumSeqs/totalNumNegSeqs), 1)
		else:
			negNumSeqs = 0
			negCov= 0
		lineList = [motifName, str(posCov), str(posNumSeqs), str(negCov), str(negNumSeqs)]
		lineStr = ','.join(lineList)
		outFile.write(lineStr + '\n')
		
		 
	outFile.close()



def svg_row_text(x,y,text):
	output = """<text xml:space="preserve" x="x_value" y="y_value" font-family="Calibri" font-size="15" text-anchor="end">text_value</text>"""
	output = output.replace("x_value",str(x))
	output = output.replace("y_value",str(y))
	output = output.replace("text_value",text)
	return output
def svg_col_text(x,y,text):
	output = """<text xml:space="preserve" x="x_value" y="y_value" font-family="Calibri" font-size="15" transform="rotate(270 x_value y_value)">text_value</text>"""
	output = output.replace("x_value",str(x))
	output = output.replace("y_value",str(y))
	output = output.replace("text_value",text)
	return output
def svg_numeric_polygon(x1,x2,y1,y2,numeric):
	# color_value = 4*(100-int(auc_value*100))
	# color_value = int(1/(auc_value*auc_value)*10)
	text = str(numeric)
	if text <= 8:
		color_value = int(255-text*70) 
	elif text <20:
		color_value = int(255-text*100)
	else:
		color_value = int(255-text*150)
	output = """<polygon points="x1,y1 x2,y1 x2,y2 x1,y2 " fill="rgb(255,color_value,color_value)"/><text xml:space="preserve" x="pos_x" y="pos_y" font-family="Calibri" font-size="7.5" fill="rgb(0,0,0)"  font-weight="bold" text-anchor="middle">text_value</text>"""
	output = output.replace("x1",str(x1))
	output = output.replace("x2",str(x2))
	output = output.replace("y1",str(y1))
	output = output.replace("y2",str(y2))
	output = output.replace("pos_x",str((x1+x2)/2))
	output = output.replace("pos_y",str((y1+y2)/2+3))
	output = output.replace("color_value",str(color_value))
	output = output.replace("text_value",str(text))
	return output

def read_mhit(file):
	dict={}
	lines = open(file).readlines()
	length = len(lines)
	temp_name = ""
	for i in range(length):
		line = lines[i].strip()
		if ">" in line:
			array = line.split(">")
			temp_name = array[-1]
			dict[temp_name] = []
		else:
			dict[temp_name].append(line)
				
	return dict

	
	
def outputResults(outFileName, depthDict, backFimoDict, filterMinNumSeq, idMotifDict, foreNumSeqs, backNumSeqs):
	"""Output the results from the depth dictionary
	Args:
		outFileName: name of file to output results to
	"""
	# a list that stores the fianl list of IDs that passed the threshold filter
	finalSelectMotifList = []
	depthOut = open(outFileName, 'wb')
	for depthVal in sorted(depthDict.iterkeys(), key=int):
		print 'depth is:', depthVal
		depthOut.write('\nDepth:' + str(depthVal) + '\n')
		depthOut.write('#Name,foreground_motif_seqCount,foreground_motif_seqCoverage(%),foreground_seqs_added,cumulative_foreground_seq_cov(%)' 
		   + ',background_motif_seqCount,background_motif_seqCoverage(%),background_seqs_added,cumulative_background_seq_cov(%)' +'\n')
		cumCov = 0
		backCumCov = 0
		backSeqsAdded = []
		motifIdList = depthDict[depthVal][0]
		motifSeqSet = depthDict[depthVal][1]
		newSeqsAddedList = depthDict[depthVal][2]
		#add new information to the motif objects as well
		for i in range(len(motifIdList)):
			lineList = []
			motifId = motifIdList[i]
			motifSeqset = motifSeqSet[i]
			seqsAdded = newSeqsAddedList[i]
			#apply the filter threshold
			if seqsAdded < filterMinNumSeq:
				continue
			if motifId not in idMotifDict:
				motifName = motifId
			else:
				motifName = idMotifDict[motifId]
			finalSelectMotifList.append(motifName)
			
			lineList.append(motifName)
			lineList.append(str(len(motifSeqset)))
			cov = 100 * (len(motifSeqset)/foreNumSeqs)
			lineList.append(str(cov))
			lineList.append(str(seqsAdded))
			cumCov += (seqsAdded/foreNumSeqs)*100
			lineList.append(str(cumCov))
			#add the background information for this motif
			#in the backFimoDict the IDs are string type
			
			if motifName in backFimoDict:
				backSeqList = backFimoDict[motifName]
				backSeqCount = len(backSeqList)
				backSeqCov = (backSeqCount/backNumSeqs)*100
				backSeqsAddedCount = 0
				for backSeq in backSeqList:
					if backSeq not in backSeqsAdded:
						backSeqsAddedCount += 1
						backSeqsAdded.append(backSeq)
				backCumCov += (backSeqsAddedCount/backNumSeqs)*100
			else:
				backSeqCount = 0
				backSeqCov = 0
				backSeqsAddedCount=0
				backCumCov += 0
			lineList.append(str(backSeqCount))
			lineList.append(str(backSeqCov))
			lineList.append(str(backSeqsAddedCount))
			lineList.append(str(backCumCov))
			line = ','.join(lineList)
			depthOut.write(line + '\n')
	depthOut.close()
	return finalSelectMotifList
	
	
	
	
	
	
	