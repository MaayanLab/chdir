# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 11:15:29 2015

@author: Luke
"""
import numpy as np
from chdir import chdir

with open('input.txt','r') as nf:
    next(nf)
    genes = []
    ctrls = []
    expms = []
    for line in nf:
        words = line.strip('\r\n').split(' ')
        genes.append(words[0])
        ctrls.append([float(x) for x in words[1:4]])
        expms.append([float(x) for x in words[4:7]])

b = chdir(ctrls,expms)
b = np.squeeze(b)
annot_b = list(zip(genes,b.tolist()))
res = sorted(annot_b,key=lambda x:x[1]*x[1],reverse=True)
assert(res[0][0]=='MCL1')
assert(res[1][0]=='LIMD2')
assert((res[1][1]-0.379125)<0.0001)