
'''
output pwm files


9-29-2015
Yichao combines Rami's code


3-30-2016
added gkm-SVM
motif_discovery is coded all within this script.

'''
from __future__ import division

import re
from utils import motif_utils
import os
import shutil
#from Bio import motifs
# from Bio.Motif import Motif as motifs
from Bio import motifs
from Bio.Seq import Seq
from utils.general_utils import *
from utils import general_utils
import joblib
from multiprocessing import Pool, cpu_count
# import subprocess
# ./convert-matrix -from meme -to transfac -i /home/working/program/Emotif3/test/test_all_motifs.pwm -o test.pssm -v
# ./convert-matrix -from transfac -to cb -i test_clustering3_cluster_root_motifs.tf -pseudo 1 -multiply 1  -decimals 1 -perm 0 -bg_pseudo 0.01 -return counts -to cb -o test3.cb
# ./matrix-clustering -v 1 -max_matrices 900 -quick -matrix_format transfac -i test.pssm -hclust_method average -title 'sdf' -motif_collection_name 'asf' -metric_build_tree 'Ncor' -lth w 5 -lth cor 0.8 -lth Ncor 0.8 -label_in_tree name -return root_matrices  -o test_clustering3
PATH = os.path.dirname(os.path.abspath(__file__))
def convolution(jid,confDict):
	convo_PWM_file = jid + "_convo_motifs.pwm"
	
	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	pos_dict = read_fasta(pos_seq)
	pos_size = len(pos_dict)
	pos_gene_list = pos_dict.keys()
	neg_dict = read_fasta(neg_seq)
	neg_size = len(neg_dict)
	neg_gene_list = neg_dict.keys()
	
	stride = int(confDict['convolution']['stride'])
	receptive_field_size = int(confDict['convolution']['receptive_field_size'])
	stop_flag = False
	
	MD_convo_folders = []						
	MD_convo_dir = jid+"_Convolution_Motif_Discovery"						
	os.makedirs(MD_convo_dir)						
	PWM_full_list = []
	for i in range(0,pos_size,stride):
		if stop_flag:
			break
		index_end = receptive_field_size+i
		current_pos_list = pos_gene_list[i:index_end]
		if index_end > pos_size:
			current_pos_list = pos_gene_list[i:pos_size] + pos_gene_list[0:(index_end-pos_size)]
			stop_flag = True
		neg_index_end = neg_size*(index_end/pos_size)
		neg_index_start = neg_size*(i/pos_size)
		current_neg_list = neg_gene_list[int(neg_index_start):int(neg_index_end)]
		
	
		pos_seq_convo = jid + "_convo_"+str(i)+"_pos.fa"
		neg_seq_convo = jid + "_convo_"+str(i)+"_neg.fa"
		
		confDict['input']['pos_seq'] = pos_seq_convo	
		confDict['input']['neg_seq'] = neg_seq_convo
		f_pos = open(pos_seq_convo,"wb")
		for s in current_pos_list:
			print >>f_pos,(">"+s)
			print >>f_pos , pos_dict[s]
		f_pos.close()
		f_neg = open(neg_seq_convo,"wb")
		for s in current_neg_list:
			print >>f_neg,(">"+s)
			print >>f_neg , neg_dict[s]
		f_neg.close()
		PWM_file,MD_dir = start_motif_discovery(jid+"_convo_"+str(i),confDict)
		# PWM_full_list.append(PWM_file)
		MD_convo_folders.append(pos_seq_convo)
		MD_convo_folders.append(neg_seq_convo)
		MD_convo_folders.append(MD_dir)
		MD_convo_folders.append(PWM_file)
		### added for fast pooling
		pwm_file,MD_dir = sub_pooling(jid+"_pooling_"+str(i),confDict)
		PWM_full_list.append(pwm_file)
		MD_convo_folders.append(MD_dir)
		MD_convo_folders.append(pwm_file)
	
	count = False						
	with open(convo_PWM_file, 'w') as outfile:						
		for fname in PWM_full_list:					
			with open(fname) as infile:				
				flag = True			
				for line in infile:			
					if flag & count:		
						if 'MOTIF' in line:	
							flag = False
							outfile.write("\n")
							outfile.write(line)
					else:		
								
						outfile.write(line)	
				count = True			
	confDict['input']['pwm_file'] = convo_PWM_file						
	confDict['input']['all_pwm_file'] = convo_PWM_file
	for outFile in MD_convo_folders:						
		print outFile					
		shutil.move(outFile, MD_convo_dir)					
	confDict['input']['pos_seq'] = pos_seq
	confDict['input']['neg_seq'] = neg_seq
	return convo_PWM_file,MD_convo_dir
	
		
def pooling(jid,confDict):
	MD_pool_dir = jid+"_Pooling_Motif_Discovery"
	try:
		os.makedirs(MD_pool_dir)
	except:
		print MD_pool_dir,"exist!"
	fileList = []
	# meme to transfac
	pwm_file = confDict['input']['pwm_file']
	output_tf = jid+"_meme2tf.tf"
	fileList.append(output_tf)
	meme_2_tf_command = "convert-matrix -from meme -to transfac -i input -o output"
	meme_2_tf_command = meme_2_tf_command.replace("input",pwm_file)
	meme_2_tf_command = meme_2_tf_command.replace("output",output_tf)
	os.system(meme_2_tf_command)
	print meme_2_tf_command
	print "======================= meme_2_tf_command ================"
	# clustering
	output_clustering = jid
	
	matrix_clustering_command = "matrix-clustering -max_matrices 2000 -quick -matrix_format transfac -i input -hclust_method average -title 'clustering_title' -motif_collection_name 'motif_collection' -metric_build_tree 'Ncor' -lth w 5 -lth cor 0.7 -lth Ncor 0.6 -label_in_tree name -return root_matrices  -o output"
	matrix_clustering_command = matrix_clustering_command.replace("input",output_tf)
	matrix_clustering_command = matrix_clustering_command.replace("output",output_clustering)
	os.system(matrix_clustering_command)
	output_clustering += "_cluster_root_motifs.tf"
	fileList.append(output_clustering)
	print matrix_clustering_command
	print "=============== matrix_clustering_command =============="
	# transfac to cb, cb is similar to pfm
	output_cb = jid+"_tf2cb.cb"
	fileList.append(output_cb)
	tf_2_cb_command = "convert-matrix -from transfac -to cb -i input -pseudo 1 -multiply 1 -decimals 1 -perm 0 -bg_pseudo 0.01 -return counts -to cb -o output"
	tf_2_cb_command = tf_2_cb_command.replace("input",output_clustering)
	tf_2_cb_command = tf_2_cb_command.replace("output",output_cb)
	os.system(tf_2_cb_command)	
	print tf_2_cb_command
	print "========= tf_2_cb_command    ================="
	# pfm to pwm
	pwmFileName = jid + '_pooled_motifs.pwm'
	convertCBToPwm(output_cb, pwmFileName)
	confDict['input']['pwm_file'] = pwmFileName						
	confDict['input']['all_pwm_file'] = pwmFileName
	for outFile in fileList:						
		print outFile					
		shutil.move(outFile, MD_pool_dir)					
	return pwmFileName,MD_pool_dir
	
def sub_pooling(jid,confDict):
	MD_pool_dir = jid+"_Motif_Discovery"
	try:
		os.makedirs(MD_pool_dir)
	except:
		print MD_pool_dir,"exist!"
	fileList = []
	# meme to transfac
	pwm_file = confDict['input']['pwm_file']
	output_tf = jid+"_meme2tf.tf"
	fileList.append(output_tf)
	meme_2_tf_command = "convert-matrix -from meme -to transfac -i input -o output"
	meme_2_tf_command = meme_2_tf_command.replace("input",pwm_file)
	meme_2_tf_command = meme_2_tf_command.replace("output",output_tf)
	os.system(meme_2_tf_command)
	print meme_2_tf_command
	print "======================= meme_2_tf_command ================"
	# clustering
	output_clustering = jid
	
	matrix_clustering_command = "matrix-clustering -max_matrices 900 -quick -matrix_format transfac -i input -hclust_method average -title 'clustering_title' -motif_collection_name 'motif_collection' -metric_build_tree 'Ncor' -lth w 5 -lth cor 0.6 -lth Ncor 0.4 -label_in_tree name -return root_matrices  -o output"
	matrix_clustering_command = matrix_clustering_command.replace("input",output_tf)
	matrix_clustering_command = matrix_clustering_command.replace("output",output_clustering)
	os.system(matrix_clustering_command)
	output_clustering += "_cluster_root_motifs.tf"
	fileList.append(output_clustering)
	print matrix_clustering_command
	print "=============== matrix_clustering_command =============="
	# transfac to cb, cb is similar to pfm
	output_cb = jid+"_tf2cb.cb"
	fileList.append(output_cb)
	tf_2_cb_command = "convert-matrix -from transfac -to cb -i input -pseudo 1 -multiply 1 -decimals 1 -perm 0 -bg_pseudo 0.01 -return counts -to cb -o output"
	tf_2_cb_command = tf_2_cb_command.replace("input",output_clustering)
	tf_2_cb_command = tf_2_cb_command.replace("output",output_cb)
	os.system(tf_2_cb_command)	
	print tf_2_cb_command
	print "========= tf_2_cb_command    ================="
	# pfm to pwm
	pwmFileName = jid + '_pooled_motifs.pwm'
	convertCBToPwm(output_cb, pwmFileName)
	# confDict['input']['pwm_file'] = pwmFileName						
	# confDict['input']['all_pwm_file'] = pwmFileName
	for outFile in fileList:						
		print outFile					
		shutil.move(outFile, MD_pool_dir)					
	return pwmFileName,MD_pool_dir
	
	
	

def start_motif_discovery(jid,confDict):
	### motif discovery ###							
	if confDict['job_type']['motif_discovery'] == 'true':							
		MD_folders = []						
		MD_dir = jid+"_Motif_Discovery"						
		os.makedirs(MD_dir)						
								
		PWMfileList = [] 						
		if confDict['motif_discovery']['gimme'] == 'true':						
			output_folder,file = run_gimme(jid,confDict)					
			PWMfileList.append(file)					
			MD_folders.append(output_folder)					
								
		if confDict['motif_discovery']['decod'] == 'true':						
			output_folder,file = run_DECOD(jid,confDict)					
			PWMfileList.append(file)					
			MD_folders.append(output_folder)					
								
		if confDict['motif_discovery']['dme'] == 'true':						
			output_folder,file = run_DME(jid,confDict)					
			PWMfileList.append(file)					
			MD_folders.append(output_folder)					
		if confDict['motif_discovery']['gkm_svm'] == 'true':						
			output_folder,file_list = run_gkm_SVM(jid,confDict)					
			PWMfileList += file_list					
			MD_folders.append(output_folder)					
		if confDict['motif_discovery']['info_gibbs'] == 'true':						
			output_folder,file_list = run_info_gibbs(jid,confDict)					
			PWMfileList += file_list					
			MD_folders.append(output_folder)					
		# print PWMfileList						
		# combine all three PWM files into one PWM file						
		count = False
		output_PWM_file = jid + "_all_motifs.pwm"
		with open(output_PWM_file, 'w') as outfile:						
			for fname in PWMfileList:					
				with open(fname) as infile:				
					flag = True			
					for line in infile:			
						if flag & count:		
							if 'MOTIF' in line:	
								flag = False
								outfile.write("\n")
								outfile.write(line)
						else:		
								
							outfile.write(line)	
					count = True			
		confDict['input']['pwm_file'] = output_PWM_file						
		confDict['input']['all_pwm_file'] = output_PWM_file						
								
								
								
		# Move Gimme, DECOD, DME folders into motif discovery folders						
		for outFile in MD_folders:						
			print outFile					
			shutil.move(outFile, MD_dir)					
								
		# fileList.append(output_PWM_file)						
		# fileList.append(MD_dir)						
	return output_PWM_file,MD_dir
	
def run_info_gibbs(jid, confDict):

	pos_seq = confDict['input']['pos_seq']
	strand = confDict['input']['strand']
	
	gibbs_folder_list = []
	gibbs_output_dir = jid + '_info_gibbs_motifs'
	os.makedirs(gibbs_output_dir)	
	
	motif_size = confDict['info_gibbs']['motif_size']
	motif_size = motif_size.split()
	num_inter = confDict['info_gibbs']['num_inter']
	num_run = confDict['info_gibbs']['num_run']
	num_motifs = confDict['info_gibbs']['num_motifs']
	zoops = confDict['info_gibbs']['zoops']
	temperature = confDict['info_gibbs']['temperature']
	
	gibbs_pwm_file_list = []
	gibbs_command = "/home/working/program/rsat/bin/info-gibbs --strand=+- -i input -w motif_size --iter=num_inter --nrun=num_run --motifs=num_motifs --zoops -t temperature "
	gibbs_command = gibbs_command.replace("input",pos_seq)
	gibbs_command = gibbs_command.replace("num_inter",num_inter)
	gibbs_command = gibbs_command.replace("num_run",num_run)
	gibbs_command = gibbs_command.replace("num_motifs",num_motifs)
	gibbs_command = gibbs_command.replace("temperature",temperature)
	gibbs_command_list = []
	
	if zoops != "true":
		gibbs_command = gibbs_command.replace("--zoops","")
	for l in motif_size:
		outputFile = jid + '_info_gibbs_' + str(l) + '.gibbs'
		gibbs_folder_list.append(outputFile)
		try:
			current_command = gibbs_command.replace("motif_size",l)
			current_command += " > "+outputFile
			# os.system(current_command)
			
			
		except:
			print 'gibbs_command execution failed. Exiting'
			exit()
		# convert matrix		
		convert_command = "/home/working/program/rsat/perl-scripts/convert-matrix -i " + outputFile + " -from infogibbs -to cb  -pseudo 1 -multiply 1 -decimals 1 -perm 0 -bg_pseudo 0.01 -return counts -to cb -o output"
		
		output_converted_File = jid + '_info_gibbs_' + str(l) + '.cb'
		gibbs_folder_list.append(output_converted_File)
		convert_command = convert_command.replace("output",output_converted_File)
		# os.system(convert_command)
		gibbs_command_list.append(current_command+";"+convert_command)
		pwmFileName = gibbs_output_dir+"/" + 'Info_gibb_motifs_'+str(l)+'.pwm'
		# convertCBToPwm(output_converted_File, pwmFileName)
		gibbs_pwm_file_list.append(pwmFileName)
	
	# run parallel
	pool = Pool(8)
	pool.map(os.system,gibbs_command_list)
	pool.close()
	pool.join()
	
	for l in motif_size:
		output_converted_File = jid + '_info_gibbs_' + str(l) + '.cb'
		pwmFileName = gibbs_output_dir+"/" + 'Info_gibb_motifs_'+str(l)+'.pwm'
		convertCBToPwm(output_converted_File, pwmFileName)
	for outFile in gibbs_folder_list:
		shutil.move(outFile, gibbs_output_dir)
	
	return gibbs_output_dir,gibbs_pwm_file_list

def run_gimme(jid, confDict):

	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	output_folder = jid + '_gimme_motifs'
	os.makedirs(output_folder)
	if confDict['motif_discovery']['convolution'] != 'true':
		
		neg_dict = read_fasta(neg_seq)
		neg_seq_gimme = jid + "_gimme_neg.fa"
		f_neg = open(neg_seq_gimme,"wb")
		for s in neg_dict.keys():
			print >>f_neg,(">"+s)
			print >>f_neg , neg_dict[s]
		f_neg.close()
		neg_seq = neg_seq_gimme
		shutil.move(neg_seq, output_folder)
		neg_seq = output_folder + "/" + neg_seq_gimme
	strand = confDict['input']['strand']
	
	
	m_size = confDict['gimme']['motif_size']
	fraction = confDict['gimme']['fraction']
	tools = confDict['gimme']['tools']
	
	
	'''
	-n name of the folder to be created
	-s search for single strand only
	-f fraction of sequences to be used for motif discovery
	-b user -u file.fasta use neg_seq
	-k output intermediate results
	
	'''
	
	## latest gimme, no -b user option
	gimme_command = 'gimme motifs '	+ pos_seq + \
		' -k -a '+ m_size + ' -n ' + output_folder + \
		' -b user -u ' + neg_seq + \
		' -f ' +  fraction
	if strand == 'single':
		gimme_command += ' -s '
	if tools != 'false':
		gimme_command += ' -t ' + tools
	
	print 'gimme motifs command is:'
	print gimme_command
	
	try:
		os.system(gimme_command)
	except:
		print 'gimmemotifs execution failed. Exiting'
		exit()
	
	gimmeAllMotifPfm = output_folder + '/intermediate_results/all_motifs.pfm'
	pwmFileName = output_folder + '/' + jid + '_all_gimme_motifs.pwm'
	convertFileToPwm(gimmeAllMotifPfm, pwmFileName)
	return output_folder,pwmFileName


def run_DECOD(jid, confDict):
	
	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	strand = confDict['input']['strand']
	output_folder = jid + '_DECOD_motifs'
	os.makedirs(output_folder)
	motif_size = confDict['DECOD']['motif_size']
	motif_size = motif_size.split()
	motif_number = confDict['DECOD']['motif_number']
	niter = confDict['DECOD']['niteration']
	
	'''
	
	-pos pos_seq
	-neg neg_seq
	-w motif width
	-stand single/both, both is default
	-o output file
	-nmotif <#motifs to find, default 10>
	-niter <#iterations, default 50>
	
	http://www.sb.cs.cmu.edu/DECOD/README.txt
	
	'''
	# java -jar /home/working/program/decod/DECOD-20111024.jar -nogui 
	fileList = []
	decod = 'java -jar %s/DECOD-20111024.jar -nogui '%(PATH)
	pwmFileName = output_folder + '/' + jid + '_decod_motifs.pwm'
	# fileList.append(pwmFileName)
	pwmFile = open(pwmFileName, 'wb')
	pwmFile.write('MEME version 4.4\nALPHABET= ACGT\nstrands: + -\nBackground letter frequencies:\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n')	
	decod_command = decod+' -pos ' + pos_seq + ' -neg ' + neg_seq + \
		' -niter ' + niter + ' -nmotif ' + motif_number
	if strand == 'single':
		decod_command += ' -strand forward '
	decod_command_list = []
	output_decod_file_list = []
	for l in motif_size:
		outputFile = jid + '_DECOD_' + str(l) + '.decod'
		output_decod_file_list.append(outputFile)
		try:
			decod_command_temp = decod_command + ' -w ' + str(l) + \
				' -o ' + outputFile   
			print "decod_command:"
			print decod_command_temp
			# os.system(decod_command_temp)
			decod_command_list.append(decod_command_temp)
			# sp = subprocess.Popen(["/bin/bash", "-i", "-c", "decod_command"])
			# sp.communicate()
			fileList.append(outputFile)
			
			# writeDecodPwm(outputFile, pwmFile, l)
		except:
			print 'decod_command execution failed. Exiting'
			exit()
	
	# run parallel
	pool = Pool(8)
	pool.map(os.system,decod_command_list)
	pool.close()
	pool.join()
	# write file
	map(lambda x:writeDecodPwm(jid + '_DECOD_' + str(l) + '.decod', pwmFile, str(x)),motif_size)
	pwmFile.close()
	# print fileList
	for outFile in fileList:
		shutil.move(outFile, output_folder)
		
		
		
	return output_folder,pwmFileName

	
def run_DME(jid, confDict):

	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	strand = confDict['input']['strand']
	output_folder = jid + '_DME_motifs'
	os.makedirs(output_folder)
	motif_size = confDict['DME']['motif_size']
	motif_number = confDict['DME']['motif_number']
	motif_size = motif_size.split()
	
	'''
	Usage: dme2 [OPTIONS] filename
	-z, --zoops                 use the ZOOPS model (default: hybrid)
	-t, --tcm                   use the TCM model (default: hybrid)
	-b, --background=STRING     background sequence file (FASTA format)
	-o, --output=STRING         output file name (default: stdout)
	-n, --number=INT            number of motifs to produce. (default: 1)
	-p, --prefix=STRING         motif accession prefix (default: "DME")
	-w, --width=INT             minimum desired motif width (default: 8)
	-i, --bits=FLOAT            min bits per column (default depends on width)
	-r, --refine=FLOAT          refinement granularity (default depends on width)
	-a, --adjust=FLOAT          adjust contribution of fg and bg
	-C, --changes=INT           changes per refinement (default: 1)
	-v, --verbose               print more run information

	Help options:
	-?, --help                  Show this help message
		--usage                 Display brief usage message

	'''	
	fileList = []
	#init a list to store name of dme output files
	dmeOutFileList = []
	dme_command1 = 'dme2 -b ' + neg_seq + ' -n ' + motif_number
	# if strand == 'single':
		# dme_command1 = dme_command1.replace('dme2','dme2_both')
	dme_command_list = []
	for l in motif_size:
		dme_command = dme_command1
		outputFile = jid + '_DME_' + str(l)
		try:
			dme_command += ' -w ' + str(l) + \
				' -o ' + outputFile + ' ' + pos_seq
			print "dme_command:"
			print dme_command
			# os.system(dme_command)
			dme_command_list.append(dme_command)
			fileList.append(outputFile)
			dmeOutFileList.append(outputFile)
		except:
			print 'dme_command execution failed. Exiting'
	# run parallel
	pool = Pool(8)
	pool.map(os.system,dme_command_list)
	pool.close()
	pool.join()
	
	motifDict = processDmeFiles(dmeOutFileList)
	pwmFileName = output_folder + '/' + jid + '_dme_motifs.pwm'
	general_utils.writePWMFromMotifs(motifDict, pwmFileName,"DME")
	# shutil.move(pwmFileName, output_folder)
	
	for outFile in fileList:
		shutil.move(outFile, output_folder)
	
	return output_folder,pwmFileName


def run_gkm_SVM(jid, confDict):



	pos_seq = confDict['input']['pos_seq']
	neg_seq = confDict['input']['neg_seq']
	strand = confDict['input']['strand']
	motif_size = confDict['gkm_svm']['motif_size']
	motif_size = motif_size.split()
	
	
	output_folder = jid + '_gkm_SVM_motifs'
	os.makedirs(output_folder)
	fileList=[]
	
	if confDict['gkm_svm']['default'] == "true":
		kmer_size = 10
		maxMismatch = 3
		filterType = 1
		numTreads = 2
		informative_columns = 6
		alpha = 3.0
		nMaxPWM = 5
		nMinKmers = 10
		top_frac = 1
		cutoff = 5
	else:
		kmer_size = confDict['gkm_svm']['kmer_size']
		maxMismatch = confDict['gkm_svm']['maxmismatch']
		filterType = confDict['gkm_svm']['filtertype']
		numTreads = confDict['gkm_svm']['numtreads']
		informative_columns = confDict['gkm_svm']['informative_columns']
		alpha = float(confDict['gkm_svm']['alpha'])
		nMaxPWM = confDict['gkm_svm']['nmaxpwm']
		nMinKmers = confDict['gkm_svm']['nminkmers']
		cutoff = int(confDict['gkm_svm']['cutoff'])
		top_frac = int(confDict['gkm_svm']['top_frac'])
	
		
	
	
	
	
	# compute kernel
	kernel_out = jid+".kernel"
	fileList.append(kernel_out)
	kernel_command = "gkmsvm_kernel -l kmer_size -k informative_columns -d maxMismatch -t filterType -R -T  numTreads pos_fasta neg_fasta kernel_out"
	if strand != "single":
		kernel_command = kernel_command.replace("-R","")
	
	kernel_command = kernel_command.replace("kmer_size",str(kmer_size))
	kernel_command = kernel_command.replace("informative_columns",str(informative_columns))
	kernel_command = kernel_command.replace("maxMismatch",str(maxMismatch))
	kernel_command = kernel_command.replace("filterType",str(filterType))
	kernel_command = kernel_command.replace("numTreads",str(numTreads))
	kernel_command = kernel_command.replace("pos_fasta",str(pos_seq))
	kernel_command = kernel_command.replace("neg_fasta",str(neg_seq))
	kernel_command = kernel_command.replace("kernel_out",str(kernel_out))
	os.system(kernel_command)
	
	# train SVM
	SVM_train_out = jid+"_svmtrain"
	SVM_train_command = "gkmsvm_train kernel_out pos_fasta neg_fasta SVM_train_out"
	SVM_train_command = SVM_train_command.replace("kernel_out",str(kernel_out))
	SVM_train_command = SVM_train_command.replace("pos_fasta",str(pos_seq))
	SVM_train_command = SVM_train_command.replace("neg_fasta",str(neg_seq))
	SVM_train_command = SVM_train_command.replace("SVM_train_out",str(SVM_train_out))
	os.system(SVM_train_command)
	alpha_out = SVM_train_out+"_svalpha.out"
	seq_out = SVM_train_out+"_svseq.fa"
	fileList.append(alpha_out)
	fileList.append(seq_out)
	# make k mer
	# kmer_command = "nrkmers.py 10 10mers.fa"
	kmer_command = "nrkmers.py "+str(kmer_size)+" kmers.fa"
	os.system(kmer_command)
	fileList.append("kmers.fa")
	# classify SVM
	kmers_scores = jid + "_kmer_scores.txt"
	SVM_classify_command = "gkmsvm_classify -d maxMismatch -t filterType -R kmers.fa svseq_fa alpha_out kmers_scores"
	SVM_classify_command = SVM_classify_command.replace("maxMismatch",str(maxMismatch))
	SVM_classify_command = SVM_classify_command.replace("filterType",str(filterType))
	SVM_classify_command = SVM_classify_command.replace("svseq_fa",str(seq_out))
	SVM_classify_command = SVM_classify_command.replace("alpha_out",str(alpha_out))
	SVM_classify_command = SVM_classify_command.replace("kmers_scores",str(kmers_scores))
	fileList.append(kmers_scores)
	os.system(SVM_classify_command)
	# kerm to PWM
	emalign_command = "svmw_emalign.py -i 40 -a alpha -n nMaxPWM -c cutoff -m nMinKmers -f top_frac -p prefix kmers_scores motif_size kmer_pwm"
	PWM_file_list = []
	gkm_svm_command_list = []
	for l in motif_size:
		kmer_pwm = output_folder+"/"+jid + "_size_" + str(l)
		temp = emalign_command.replace("alpha",str(alpha))
		temp = temp.replace("nMaxPWM",str(nMaxPWM))
		temp = temp.replace("cutoff",str(cutoff))
		temp = temp.replace("nMinKmers",str(nMinKmers))
		temp = temp.replace("top_frac",str(top_frac))
		temp = temp.replace("kmers_scores",kmers_scores)
		temp = temp.replace("kmer_pwm",kmer_pwm)
		temp = temp.replace("motif_size",str(l))
		temp = temp.replace("prefix","length."+str(l))
		# os.system(temp)
		gkm_svm_command_list.append(temp)
		out_file = kmer_pwm+"_models.meme"
		PWM_file_list.append(out_file)
		alpha += 2
		#cutoff += 2

		
	# run parallel
	pool = Pool(8)
	pool.map(os.system,gkm_svm_command_list)
	pool.close()
	pool.join()		
	for outFile in fileList:
		shutil.move(outFile, output_folder)
	return output_folder,PWM_file_list
		

def convertFileToPwm(pfmFile, pwmFileName):
	# Rami's code
	""" read the pfm file produced by gimmemotifs and convert it to a pwm file"""
	#a dict between motif name and its list of pfm lines
	motifPfmDict = {}
	#a dict between motif name and its pwm
	motifPwmDict = {}
	flag = 0
	with open(pfmFile, 'rb') as handler:
		for line in handler:
			line = line.strip()
			if re.search(r'>', line):
				motifName = line[1:]
				motifPfmDict[motifName] = []
				continue
			motifPfmDict[motifName].append(line)
			
	#check the motif dict
	for motifName, pfmList in motifPfmDict.items():
		# print motifName,len(pfmList)
		#print pfmList
		nucs = {"A":0,"C":1,"G":2,"T":3}
		pfm = [[0 for x in range(4)] for x in range(len(pfmList))]
		#print pfm
		for rowIndex in range(len(pfmList)):
			row = pfmList[rowIndex]
			split = row.split('\t')
			#print 'spl:',split
			for j in range(4):
				pfm[rowIndex][j] = float(split[j])
				#print rowIndex,j,split[j]
			
		#pwm = pfm_to_pwm(pfm)
		motifPwmDict[motifName] = pfm_to_pwm(pfm)
		#break
	
	pwmFile = open(pwmFileName, 'wb')
	pwmFile.write('MEME version 4.4\nALPHABET= ACGT\nstrands: + -\nBackground letter frequencies:\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n')
	for motifName, pwm in motifPwmDict.items():
		pwmFile.write('\n\nMOTIF ' + motifName + '\n')
		pwmFile.write('letter-probability matrix: alength= 4 w= '+ str(len(pwm)) + ' nsites= 100 E= 0.001\n') 
		pwmLines = "\n".join(["\t".join(["%s" % x for x in row]) for row in pwm])
		pwmFile.write(pwmLines)
	
	pwmFile.close()
	return 1
	
def convertCBToPwm(pfmFile, pwmFileName):
	# modify Rami's code
	""" read the pfm file produced by gimmemotifs and convert it to a pwm file"""
	#a dict between motif name and its list of pfm lines
	motifPfmDict = {}
	#a dict between motif name and its pwm
	motifPwmDict = {}
	flag = 0
	pattern = re.compile(r'^>.*consensus=(?P<name>[A-Z]+).*size=(?P<size>[0-9]+)')
	with open(pfmFile, 'rb') as handler:
		for line in handler:
			line = line.strip()
			match = pattern.search(line)
			if match:
				name = match.group('name')
				size = match.group('size')
				
				motifName = name+"_"+size 
				motifPfmDict[motifName] = []
				continue
			motifPfmDict[motifName].append(line)
			
	#check the motif dict
	for motifName, pfmList in motifPfmDict.items():
		# print motifName,len(pfmList)
		#print pfmList
		nucs = {"A":0,"C":1,"G":2,"T":3}
		pfm = [[0 for x in range(4)] for x in range(len(pfmList))]
		#print pfm
		for rowIndex in range(len(pfmList)):
			row = pfmList[rowIndex]
			split = row.split('\t')
			#print 'spl:',split
			for j in range(4):
				pfm[rowIndex][j] = float(split[j])
				#print rowIndex,j,split[j]
			
		#pwm = pfm_to_pwm(pfm)
		motifPwmDict[motifName] = pfm_to_pwm(pfm)
		#break
	
	pwmFile = open(pwmFileName, 'wb')
	pwmFile.write('MEME version 4.4\nALPHABET= ACGT\nstrands: + -\nBackground letter frequencies:\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n')
	for motifName, pwm in motifPwmDict.items():
		pwmFile.write('\n\nMOTIF ' + motifName + '\n')
		pwmFile.write('letter-probability matrix: alength= 4 w= '+ str(len(pwm)) + ' nsites= 100 E= 0.001\n')  
		pwmLines = "\n".join(["\t".join(["%s" % x for x in row]) for row in pwm])
		pwmFile.write(pwmLines)
	
	pwmFile.close()
	return 1
	
def pfm_to_pwm(pfm, pseudo=0.001):
	"""Used form the Gimmemotifs tool """
	return [[(x + pseudo)/(float(sum(row)) + pseudo * 4) for x in row] for row in pfm]

def writeDecodPwm(inFile, outFile, wordLen):
	"""Read the motifs and put them in meme PWM format """
	flag = 0
	#dict between motif names and list of pfm lines
	motifDict = {}
	motifName = ''
	count = 0
	with open(inFile, 'rb') as handler:
		for line in handler:
			line = line.strip()
			if re.search(r'>Motif', line):
				flag = 1
				motifName = line[1:]
				motifDict[motifName] = []
				#motifDict[motifName].append([])
				count = 0
				continue
			if flag == 1:
				newLine = line[3:len(line)-1]
				split = newLine.split()
				#for val in split:
					#motifDict[motifName][count].append(val)
				motifDict[motifName].append(split)
				count += 1
			if count == 4:
				flag = 0
	
		
	#check the motif dict
	for motifName, matrix in motifDict.items():
		#print motifName
		outFile.write('\n\nMOTIF DECOD_' + motifName + '_' + str(wordLen) + '\n')
		
		#print '\t', matrix
		trans = zip(*matrix)
		# print '\t', trans
		# outFile.write('letter-probability matrix: alength= 4 w= '+ str(len(trans)) + '\n') 
		outFile.write('letter-probability matrix: alength= 4 w= '+ str(len(trans)) + ' nsites= 100 E= 0.001\n') 
		for row in trans:
			print '\t\t', row
			for val in row:
				print '\t\t\t', val
				outFile.write(str(val) + ' ')
			outFile.write('\n')
		
	return 1

def parseDME(outFile):
	"""
	
	Args:
		outFile: name of DME output file
	
	Returns:
		 a motif dict between motif Ids (motif names as specified in the DME output file) and a BioPython motif object
	
	"""
	motifId = ''
	motifScore = 0
	#dict between a motif rank (ID) and a BioPython motif object
	motifDict = {}
	motifRank = 0
	#list that stores the binding sites of the motif
	siteList = []
	#list that stores the seqs that the motif occurs in
	motifSeqList = []
	#list that stores word objs per motif
	wordObjList = []
	maxNumMotifs = 5
	with open(outFile, 'rd') as handler:
		for line in handler:
			line = line.strip()
			if re.search(r'ID', line):
				#make the motif instance/object
				if len(siteList) != 0:
					motifRank += 1
					#initialize the MyMotif object
					motifDict[motifId] = motif_utils.MyMotif()
					#make and store the BioPython motif object
					# print [x.upper() for x in siteList]
					motifDict[motifId].bioMotifObj = motifs.create([x.upper() for x in siteList])
					motifDict[motifId].regExp = motifId
					motifDict[motifId].score = motifScore
					motifDict[motifId].foreSeqList = list(motifSeqList)#make a copy of the list
					#add the motif ID to each word obj
					for wordObj in wordObjList:
						wordObj.motifId = motifId
					motifDict[motifId].foreWordObjList = list(wordObjList)
					#re initialize for next motif
					siteList = []
					motifSeqList = []
					wordObjList = []
				#if motifRank == maxNumMotifs:
					#break
				split = line.split()
				motifId = split[1]
				motifId = motifId.strip()
				continue
			if re.search(r'SCORE', line):
				split = line.split()
				score = split[1]
				split = score.split('=')
				motifScore = split[1]
				continue
			if re.search(r'BS', line):
				split_1 = line.split(';')
				site = split_1[0]
				seq = split_1[1].strip()
				if seq not in motifSeqList:
					motifSeqList.append(seq)
				position = split_1[2].strip()
				strand = split_1[5].strip()
				#siteSplit
				split = site.split()
				site = split[1]
				site = site.strip()
				mySite = Seq(site)
				siteList.append(mySite)
				#make a word obj
				wordObj = motif_utils.MyWord()
				wordObj.string = site
				wordObj.position = position
				wordObj.strand = strand
				wordObj.seq = seq
				wordObjList.append(wordObj)
				continue
			
				
	
	#store the last motif info
	if len(siteList) != 0:
		motifRank += 1
		motifDict[motifId] = motif_utils.MyMotif()
		#make and store the BioPython motif object
		# print (siteList)
		motifDict[motifId].bioMotifObj = motifs.create([x.upper() for x in siteList])
		motifDict[motifId].regExp = motifId
		motifDict[motifId].score = motifScore
		motifDict[motifId].foreSeqList = list(motifSeqList)#make a copy of the list
		#add the motif ID to each word obj
		for wordObj in wordObjList:
			wordObj.motifId = motifId
		motifDict[motifId].foreWordObjList = list(wordObjList)
	
	##check the motif dict
	#counter = 3
	#for motifRank, motifObj in motifDict.items():
		#print 'rank:', motifRank,'score:', motifObj.score, 'regExp:',motifObj.regExp
		#print motifObj.bioMotifObj
		##check the seq names
		#for seqName in motifObj.foreSeqList:
			#print '\t', seqName
		##check the word information
		#for wordObj in motifObj.foreWordObjList:
			#print '\t', wordObj.motifId, wordObj.string, wordObj.strand, wordObj.position, wordObj.seqName
		#counter += 1
		#if counter == 1:
			#break
				
	return motifDict		


def processDmeFiles(dmeFileList):
	"""Given a list of DME output files, read them and store them as motif PWMs
	Args:
		dmeFileList: list of the DME output file names per motif length
	Returns:
		a motif dictionary of motifs of all lengths 
	"""
	motifDict = {}
	#motifId = 1
	for fileName in dmeFileList:
		# print "parsing:", fileName
		tmpMotifDict = parseDME(fileName)
		for mId in tmpMotifDict:
			motifDict[mId] = tmpMotifDict[mId]
			#motifId += 1
	
	
	return motifDict









