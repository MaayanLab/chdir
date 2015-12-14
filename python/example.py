from chdir import chdir

# The example.txt is a gene expression file. In the file there are a control 
# group of 20 replicates and an experiment group of 6 replicates. This script 
# shows how to calculate the characteristic direction vector from the data
# using the chdir module.

# Author: Qiaonan Duan
# Ma'ayan Lab, Icahn School of Medicine at Mount Sinai
# Oct. 15, 2013
#

filename = 'example.txt'

with open(filename) as nf:
	header = next(nf).rstrip('\r\n').split('\t')
	header = header[1:]
	ctrlIdx = [i for i in range(len(header)) if header[i]=='0']
	expIdx = [i for i in range(len(header)) if header[i]=='1']
	assert((len(ctrlIdx)+len(expIdx))==len(header))
	#next(nf) #skip 2nd line

	identifiers = []
	ctrlMat = []
	expMat = []
	for line in nf:
		words = line.rstrip('\r\n').split('\t')
		identifiers.append(words[0])
		values = words[1:]
		ctrlMat.append([float(values[i]) for i in ctrlIdx])
		expMat.append([float(values[i]) for i in expIdx])


chdirVector = chdir(ctrlMat,expMat,identifiers,1)
print 'chdirVector:', chdirVector

