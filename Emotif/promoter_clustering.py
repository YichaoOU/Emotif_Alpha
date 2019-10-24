from __future__ import division
import os
import sys
import argparse
import shutil
import re
from utils import *
from utils.general_utils import *
import warnings
warnings.filterwarnings("ignore")
import sys
import numpy as np
from motif_filtering import *
from scipy.io.arff import loadarff as la
from sklearn.cluster import AgglomerativeClustering as HC
from sklearn.cluster import SpectralClustering as SC
from sklearn.metrics import silhouette_samples, silhouette_score
from copy import deepcopy as dp
# incomplete
# just trying to see an initial result
def jacards_HC(jid,confDict,csv):
	fileList = [] 
	# run R script to do clustering
	dir = os.path.dirname(os.path.realpath(__file__))
	R_script_dir = dir + "/jacard_HC_ward.R "
	HC_figure = jid + "_hc_ward.png"
	HC_clusters = jid + "_clusters.txt"
	R_command = R_script_dir + csv + " " + jid + " > " + HC_clusters
	print R_command
	os.system("Rscript " + R_command)
	fileList.append(HC_figure)
	fileList.append(HC_clusters)
	# parse the cluster assignment file
	clusters = parse_cluster_assignment(HC_clusters)
	cluster_size = map(lambda x:len(x),clusters)
	total_sum = sum(cluster_size)
	cluster_percentage = map(lambda x:x/total_sum,cluster_size)
	cluster_percentage_cum = list(np.cumsum(cluster_percentage))
	print "########## cluster_percentage_cum ##########"
	print cluster_percentage_cum
	# partition pos mhit & fasta files
	pos_hit = confDict['input']['pos_mhit']
	pos_seq = confDict['input']['pos_seq']
	posMotifDict = readHitFile(pos_hit)
	posSeqDict = read_fasta(pos_seq)
	for i in range(len(clusters)):
		seq_list = clusters[i]
		# mhit = dp(posMotifDict)
		out_f = jid + "_cluster_"+str(i)+"_pos.mhit"
		out_fasta_f = jid + "_cluster_"+str(i)+"_pos.fa"
		out = open(out_f,"wb")
		out_fasta = open(out_fasta_f,"wb")
		fileList.append(out_f)
		fileList.append(out_fasta_f)
		for s in seq_list:
			print >>out_fasta,(">"+s)
			print >>out_fasta , posSeqDict[s]
		for m in posMotifDict:
			mhit = list(set(posMotifDict[m]).intersection(set(seq_list)))
			if len(mhit) >= 1:
				print >>out,(">"+m)
				print >>out,"\n".join(mhit)
		out.close()
		out_fasta.close()	
	
	# partition neg mhit & fasta files
	neg_hit = confDict['input']['neg_mhit']
	neg_seq = confDict['input']['neg_seq']
	
	negMotifDict = readHitFile(neg_hit)
	negSeqDict = read_fasta(neg_seq)
	neg_seq_list = negSeqDict.keys()
	neg_clusters = do_split(neg_seq_list,cluster_percentage_cum)
	
	for i in range(len(neg_clusters)):
		seq_list = neg_clusters[i]
		# mhit = dp(negMotifDict)
		
		out_f = jid + "_cluster_"+str(i)+"_neg.mhit"
		out_fasta_f = jid + "_cluster_"+str(i)+"_neg.fa"
		out = open(out_f,"wb")
		out_fasta = open(out_fasta_f,"wb")
		fileList.append(out_f)
		fileList.append(out_fasta_f)
		for s in seq_list:
			print >>out_fasta,(">"+s)
			print >>out_fasta , negSeqDict[s]
		for m in negMotifDict:
			mhit = list(set(negMotifDict[m]).intersection(set(seq_list)))
			if len(mhit) >= 1:
				print >>out,(">"+m)
				print >>out,"\n".join(mhit)
		out.close()
		out_fasta.close()	


	return fileList
	
	
def parse_cluster_assignment(result):
	lines = open(result).readlines()
	# print lines
	num = int((lines[0].strip().split(" ")[-1]).replace('"',''))
	# print num
	cluster = [0]*num
	for i in range(num):
		cluster[i] = []
	for l in lines[1:]:
		line = l.strip().split(" ")
		# print line
		seq_name = (line[1]).replace('"','')
		k = int((line[-1]).replace('"',''))
		cluster[k-1].append(seq_name)
	return cluster
	
def jacards_Spectral_Clustering(jid,confDict,arff):
	data = la(arff)
	# pos_seq = confDict['input']['pos_seq']
	# neg_seq = confDict['input']['neg_seq']
	# pos_name = findName(pos_seq) # fasta file
	# neg_name = findName(neg_seq)
	### get all features
	features=list(data[1])[:-1]
	# top = int(confDict['RF_gini_filter']['top'])
	# total_ranking_file = jid + "_gini_total_ranking.tsv"
	# total_ranking = open(total_ranking_file,"wb")
	# Y=np.array(data[0]["Class"])
	# print Y
	# sum_pos = np.sum(Y)
	sum_pos = 334
	X=np.array(map(lambda x:list(x),data[0][features].tolist()))
	X = X[1:sum_pos,:]
	n,m = X.shape
	print n,m
	Sim = np.ones((n,n))
	print Sim
	for i in range(n):
		# print i
		for j in range(i+1,n):
			Sim[i,j] = jacards_score(X[i],X[j],m)
			Sim[j,i] = Sim[i,j]
	print Sim
	Cluster = []
	a=SC(n_clusters=2,affinity="precomputed",assign_labels = "discretize")
	Clusters2 = a.fit_predict(Sim)
	a=SC(n_clusters=3,affinity="precomputed",assign_labels = "discretize")
	Clusters3 = a.fit_predict(Sim)
	a=SC(n_clusters=4,affinity="precomputed",assign_labels = "discretize")
	Clusters4 = a.fit_predict(Sim)
	a=SC(n_clusters=5,affinity="precomputed",assign_labels = "discretize")
	Clusters5 = a.fit_predict(Sim)
	Cluster.append(Clusters2)
	Cluster.append(Clusters3)
	Cluster.append(Clusters4)
	Cluster.append(Clusters5)
	sil_score = []
	silhouette_avg2 = silhouette_score(X, Clusters2)
	silhouette_avg3 = silhouette_score(X, Clusters3)
	silhouette_avg4 = silhouette_score(X, Clusters4)
	silhouette_avg5 = silhouette_score(X, Clusters5)
	sil_score.append(silhouette_avg2)
	sil_score.append(silhouette_avg3)
	sil_score.append(silhouette_avg4)
	sil_score.append(silhouette_avg5)
	ind = np.argmax(sil_score)
	print sil_score
	return Cluster[ind]

# vectorization
def jacards_score(x,y,m):
	# x is 1*m 
	# y is 1*m
	# jacards score is defined in Srinivas Veerla, 2006, BMC Bioinformatics
	A=0
	B=0
	C=0
	A = np.dot(x,y)
	sum_x = np.sum(x)
	sum_y = np.sum(y)
	B = (sum_x+sum_y-A)/2
	C = B
	return A/(A+B+C)



def calc_pairwise_stability( clusterings, metric ):
	''' 
	Calculate mean pairwise stability between a list of disjoint clusterings,
	using the specified metric for measuring the similarity of two disjoint clusterings.
	'''
	sim_values = []
	for i in range(len(clusterings)):
		for j in range(i+1,len(clusterings)):
			sim_values.append( metric( clusterings[i], clusterings[j] ) )
	return mean(sim_values)

def HC_clustering( X, k, sampling_ratio ):
	''' 
	Apply HC to a subset of samples from the specified dataset, 
	and return a predicted clustering for the complete dataset based on the original centroids.
	'''
	# create a matrix with subset of samples
	n_samples = X.shape[0]
	indices = np.arange(n_samples)
	np.random.shuffle( indices )
	n_subset = int(n_samples * sampling_ratio) 
	X_subset = X[indices[0:n_subset]] 
	# cluster the subset
	clusterer = KMeans(n_clusters=k, n_init=1, init='random', max_iter = 100)
	clusterer.fit(X_subset)
	# produce an assignment for all samples
	return clusterer.predict(X)
	
def do_split(x, percent):
    L = len(x)
    idx1 = [0] + list(int(L * p) for p in percent[:-1])
    idx2 = idx1[1:] + [L]
    return list(x[i1:i2] for i1,i2 in zip(idx1, idx2))
	   
# import sys
# print parse_cluster_assignment(sys.argv[1])