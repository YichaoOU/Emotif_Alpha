# R code to perform hierarchical clustering using jacard similarity and ward algorithm 

# need to test if the calculations are correct
library("cluster")
library(ggplot2)
library(ggdendro)
jacards_score = function(x,y){
	A=0
	B=0
	C=0
	A = x %*% y
	sum_x = sum(x)
	sum_y = sum(y)
	B = (sum_x+sum_y-A)/2
	C = B
	return (A/(A+B+C))
}
# library(stringi)
hc_ward = function(f,jid){
	# jid = str_replace_all(string=jid, pattern=" ", repl="")
	# this file f only contains the positive sequences, where the last column is the sequence name
	df=read.table(f,sep=",",header=T)
	size_df = dim(df)
	data = data.matrix(df[,1:(size_df[2]-1)])
	size_data = dim(data)
	# print (size_data)
	sim = matrix(1, size_data[1], size_data[1])
	for (i in 1:size_data[1]){			
		for (j in (i+1):size_data[1]){		
			if (j>size_data[1]){	
				break
			}	
			sim[i,j] = jacards_score(data[i,],data[j,])	
			sim[j,i] = sim[i,j]	
			
		}		
				
	}
	d = 1-sim
	label = df[,size_df[2]]
	rownames(d) <- label 
	d = as.dist(d)
	hc <- hclust(d, "ward.D")
	output = paste(jid,"_hc_ward.png",sep="")
	# png(filename=,res=600)
	# par(mar=c(1,1,1,1))
	# plot(hc)
	# dev.off()
	ggdendrogram(hc, rotate = FALSE, size = 1)
	ggsave(output)	
	clusters = sapply(2:15,best_k_clusters,hc=hc,d=d)
	# print ("clusters#########################")
	# print (clusters)
	index = which( clusters == max(clusters) , arr.ind = T )
	# print ("#_----------------------index")
	# print (index) 
	k = index+1
	cluster_assignment = cutree(hc,k=k)
	print (paste("Number of clusters:",k))
	for (i in 1:size_df[1]){
		print (paste(label[i],cluster_assignment[i]))
	}
	
}
library(clusterSim)

best_k_clusters = function(k,hc,d){
	hc_cut_k = cutree(hc,k=k)
	# print ("hc_cut_k")
	# print (hc_cut_k)
	# score = summary(silhouette(hc_cut_k,d))
	# print ("score$avg.width")
	# print (k)
	# print (score$avg.width)
	# output = score$avg.width
	# print ("output")
	# print (output)
	output = index.G2(d,hc_cut_k)
	# print (k)
	# print ("output")
	# print (output)
	
	return (output)
		
}

file <- commandArgs(trailingOnly = TRUE);
hc_ward(file[1],file[2])

















