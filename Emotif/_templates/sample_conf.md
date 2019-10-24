
; -------  Global settings --------
[input]
; # when entering the file names put the full path name of the files
; # the foreground and background files for motif discovery

; Positive genes and negative genes are required all the time.
pos_seq=Input_pos_fasta
neg_seq=Input_neg_fasta

; Motif discovery on both strand or single strand
strand= both

; When you already have the motifs, and you want to skip the motif discovery step and do the remaining analysis, please provide the pwm file.
pwm_file= /usr/local/lib/python2.7/dist-packages/Emotif/_dataset/jid_all_motifs.pwm

; When you already have the motifs, and you want to skip the motif discovery & scanning steps and do the remaining analysis, please provide the motif hit file.
pos_mhit=/usr/local/lib/python2.7/dist-packages/Emotif/_dataset/jid_pos.mhit
neg_mhit=/usr/local/lib/python2.7/dist-packages/Emotif/_dataset/jid_neg.mhit

; This functionality is only for in-house use.
pos_testing_seq=/usr/local/lib/python2.7/dist-packages/Emotif/_dataset/pos.fa
neg_testing_seq=/usr/local/lib/python2.7/dist-packages/Emotif/_dataset/neg.fa


pwm_selected= /usr/local/lib/python2.7/dist-packages/Emotif/_dataset/jid_all_motifs.pwm



[job_type]
; # you can control pipeline here
; # it is able to do the following things with approperiate input settings

; - motif discovery only [need two fasta files]
; - motif discovery ---> motif scanning ---> motif selection [need two fasta files]
; - motif scanning only [need two fasta files and one pwm files]
; - motif scanning ---> motif selection [need two fasta files and one pwm files]
; - motif selection only [need two mhit files and a pwm file]
; - motif selection is now a little complex, greedy algorithm requires pwm file while genetic algorithm does not. 
; - motif filtering requires [two fasta files,two mhit files]

motif_discovery = true
motif_scanning = true
sequence_clustering = false
motif_filtering = true
motif_selection = false
motif_output = false
motif_testing = false

; -------  motif_discovery settings --------

[motif_discovery]
gimme = true
DECOD = true
DME = true
gkm_SVM = false
pooling = false
convolution = false
info_gibbs = false


[convolution]
stride = 10
receptive_field_size = 20

[gimme]
motif_size = small
fraction = 0.5
tools = MEME,Weeder,BioProspector,AMD,Homer,GADEM,MDmodule,Improbizer

[DECOD]
;motif_size = 8 10 12
motif_size = 8
motif_number = 1
niteration = 1

[DME]
;motif_size = 8 10 12 14
;motif_number = 40
motif_size = 8
motif_number = 50

[gkm_svm]
motif_size = 8 10 12 14
; motif_numer is specified below, as nMaxPWM

; if you set defaul to be true, anything below won't be used
default = false

;kmer size should be less or equal to 10.
;kmer size can't be an array, e.g 8 10
kmer_size = 8
; kernel specification
maxMismatch = 2
;set filter type: 0(use full filter), 1(use truncated filter:this gaurantees non-negative counts for all L-mers), 2(use h[m],gkm count vector), 3(wildcard), 4(mismatch), default=1
filterType = 1
numTreads = 4
;set number of informative columns, default= 6 
informative_columns = 6

; kmer to motifs
;set the alpha (multiplier) for svm weights (default=3.0)
;control the degenerative of PWM
; when it is large, the PWM will be like a kmer without any degenerativity
alpha = 30.0
;set the maximum number of pwms to build (default=5)
nMaxPWM = 5
;set the minium number of kmers for each pwm (default=10)
nMinKmers = 10 
;set the percentage of the top k-mers to be evaluated as seed
;the larger, the slower
top_frac = 1
;set the cut-off of minimum log oddratio sum from the model      when model rebuilding (default=5)
cutoff = 5


[info_gibbs]
motif_size = 8 10 12 14 16 18 20
num_inter = 200
num_run = 5
num_motifs = 10
zoops = true
temperature = 1.0
;the te mperature parameter T defined in Equation (7) influences the explored space and the convergence speed of the algorithm. A low temperature (T < 1.0) allows a fast convergence but reduces the number of explored solutions. A high temperature (T > 1.0) on the other hand allows a wider exploration, but unfavorably affects the convergence speed. [The temperature should be in range [0.6 1.4]


; -------  motif_scanning settings --------

[motif_scanning]
fimo = true

[fimo]
;default p-value
pvalue = 0.0001
;starnd single or double
strand=double





; -------  sequence clustering settings --------



[coverage_filter]
maxnegcov = 0.5
minposcov = 0.1



; -------  motif_filtering settings --------
[motif_filtering]
; # not impletemented
coverage_filter = false
RF_gini_filter = true
RF_entropy_filter = true
; # not impletemented
SVM_filter = false

[coverage_filter]
maxnegcov = 0.99
minposcov = 0.0005

[RF_gini_filter]
; get top number of motifs
top = 5 

[RF_entropy_filter]
top = 5 

[SVM_filter]
top = 20 

; -------  motif_selection settings --------

[motif_selection]
; Algorithm for motif selection based on sequence coverage. Available options: greedy, ILP, branch_cut, required_cover


greedy = false
greedy_multicover = true
; Compile them on your machine before using them
iLP = false
branch_cut = false
required_cover = false

genetic_algo = false
genetic_algo_multi = false
prefilter = true
tomtom = true
filter_percentage = 0.5
;depth of search
depth=3
kmer = true
[prefilter]
maxnegcov = 0.99
minposcov = 0.0005

[tomtom]
evalue = 0.05
[genetic_algo]
; currently have no parameters
path = java -classpath "/home/working/program/motif_selection/PNPSC/:/home/liyc/working/lib/moea/MOEAFramework-2.5/lib/*" Motif_select_example

; # maximum number of iterations
numEvals=20

; # space is not allowed in the following assignments
penaltyValue=0.6
; this is the alpha value used in the cost function of Positive Negative Partial set cover (PNPSC). This is per solution. the higher alpha the more
;penalized the solution will be if it has haigh negative coverge. 0.6 worked well.

; # maximum motif coverage to be 30%
maxNegCovPerMotif=0.3
; maximum negative coverage per motif. If no motif comply weith this then it will run for a long time;


filterCutoff=0.02
; the filtering threshold where if a motif does not add this minimum number of newly added seqeunces to the set cover, it will be penalized
;5% threhold = 0.005. If there no motifs that comply with this option , then the code will run for a long time. Might be useful
;to check the motifs before assigning this value 


[genetic_algo_multi]
; currently have no parameters
path = java -classpath "/usr/local/lib/python2.7/dist-packages/emoti3/algo/genetic_multicover/:/home/liyc/working/lib/moea/MOEAFramework-2.5/lib/*" Motif_select_example

; # maximum number of iterations
numEvals=20

; # space is not allowed in the following assignments
penaltyValue=0.6
; this is the alpha value used in the cost function of Positive Negative Partial set cover (PNPSC). This is per solution. the higher alpha the more
;penalized the solution will be if it has haigh negative coverge. 0.6 worked well.

; # this is an important parameter
minPosCov_threshold=0.5

; # maximum motif coverage to be 30%
maxNegCovPerMotif=0.3
; maximum negative coverage per motif. If no motif comply weith this then it will run for a long time;


filterCutoff=0.02
; the filtering threshold where if a motif does not add this minimum number of newly added seqeunces to the set cover, it will be penalized
;5% threhold = 0.005. If there no motifs that comply with this option , then the code will run for a long time. Might be useful
;to check the motifs before assigning this value 





; -------  motif_output settings --------
[motif_output]
mast = true
solution_table = true







;below are not considered



[motif.logo]
; make motif logos
opt=true
