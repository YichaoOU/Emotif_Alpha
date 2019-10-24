# 2/6/2016 
# csv to arff format
# arff format is the input for using scikit learning
# python csv2arff.py in.csv out.arff
import sys

infile = sys.argv[1]
f = open(infile).readlines()
features = f[0].strip()
features = features.split(",")
out_file = sys.argv[2]
out = open(out_file,"wb")
print >>out,"@relation",out_file
print >>out," "
for feature in features[0:-1]:
	print >>out,"@attribute",feature,"numeric"
print >>out,"@attribute Class {pos,neg}"
print >>out," "
print >>out,"@data"

print >>out,"".join(f[1:])	








