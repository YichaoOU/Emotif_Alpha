
'''
input: fa 
output: k_mer_mhit file


12-1-2015

'''
import itertools

# assume sequences are stored in a dict
# dict[seq_name] = seq
bases=['A','T','G','C']

def get_mhit(seq_dict,k,label):
	motif_dict = {}
	kmer = [''.join(p) for p in itertools.product(bases, repeat=k)]
	output = open(label+".mhit","wb")
	for mer in kmer:
		if not motif_dict.has_key(mer):
			motif_dict[mer]=[]
		for name in seq_dict:
			seq = seq_dict[name]
			if mer in seq:
				motif_dict[mer].append(name)
	for mer in motif_dict:
		if len(motif_dict[mer]) == 0:
			continue
		print >>output,(">"+mer)
		for name in motif_dict[mer]:
			print >>output,name
	output.close()
	return 1
	
	
	
	
	
	
	
	
	
	
