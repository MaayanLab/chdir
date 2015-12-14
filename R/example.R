
# The example.txt has gene expression data. In the file there is a control 
# group of 20 replicates and an experiment group of 6 replicates. This script 
# shows how to calculate the characteristic direction vector from the data
# using the chdir module.

# Author: Qiaonan Duan
# Ma'ayan Lab, Icahn School of Medicine at Mount Sinai
# Jan. 13, 2014

source('chdir.R')
source('nipals.R')

data <- read.table('example.txt',sep="\t")
header <- as.list(data[1,2:dim(data)[2]])
genes <- as.vector(data[,1])[2:dim(data)[1]]
mat <- as.matrix(data[2:dim(data)[1],2:dim(data)[2]])

ctrlMat <- mat[,header==0]
expmMat <- mat[,header==1]
# unitV is the characteristic direction.
unitV <- chdir(ctrlMat,expmMat,genes)
