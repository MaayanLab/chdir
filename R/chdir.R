"chdir" <-
function(ctrl,expm,genes,r=1)
# This function caclulates the characteristic direction for a gene expression dataset.
#  	ctrl: control gene expressoion data, a matrix object
#  	expm: experiment gene expression data, a matrix object
#  	b: return value, a vector of n-components, representing the characteristic
#          direction of the gene expression dataset. n equals to the number of genes in the 
#          expression dataset. b is also a matrix object. b is sorted by its components' 
#          absolute values in descending order.
#  	r: regularized term. A parameter that smooths the covariance matrix and reduces
#          potential noise in the dataset. The default value for r is 1, no regularization.
#
#       For the input matrix rows are genes and columns are gene expression profiles.
#       r is the regulization term ranging [0,1]. b is the characteristic direction.
#       ctrl(control) and expm(experiment) matrices should have the same number
#       of genes(rows). 
#
#       Author: Qiaonan Duan
#       Ma'ayan Lab, Icahn School of Medicine at Mount Sinai
#       Jan.13, 2014
#
#		Add gene symbols to results. Apr. 4, 2014

{

     if(dim(ctrl)[1]!=dim(expm)[1]){
	stop('Control expression data must have equal number of genes as experiment expression data!')
}

     if(any(is.na(ctrl))||any(is.na(expm))){
	stop('Control expression data and experiment expression data have to be real numbers. NA was found!')
}


# There should be variance in expression values of each gene. If  
# gene expression values of a gene are constant, it would dramatically
# affect the LDA caculation and results in a wrong answer.
constantThreshold <- 1e-5;
ctrlConstantGenes <- diag(var(t(ctrl))) < constantThreshold
expmConstantGenes <- diag(var(t(expm))) < constantThreshold

if (any(ctrlConstantGenes)){
	errMes <- sprintf('%s row(s) in control expression data are constant. Consider Removing the row(s).',paste(as.character(which(ctrlConstantGenes)),collapse=','))
	stop(errMes)
}else if(any(expmConstantGenes)){
	errMes <- sprintf('%s row(s) in experiment expression data are constant. Consider Removing the row(s).',paste(as.character(which(expmConstantGenes)),collapse=','))
	stop(errMes)
}

# place control gene expression data and experiment gene expression data into
# one matrix
combinedData <- cbind(ctrl,expm)

# get the number of samples, namely, the total number of replicates in  control 
# and experiment. 
dims <- dim(combinedData)
samplesCount <- dims[2]

# the number of output components desired from PCA. We only want to calculate
# the chdir in a subspace that capture most variance in order to save computation 
# workload. The number is set 20 because considering the number of genes usually 
# present in an expression matrix 20 components would capture most of the variance.
componentsCount <- min(c(samplesCount-1,20))


# use the nipals PCA algorithm to calculate R, V, and pcvars. nipals algorithm
# has better performance than the algorithm used by R's builtin PCA function.
# R are scores and V are coefficients or loadings. pcvars are the variances 
# captured by each component 
pcaRes <- nipals(t(combinedData),componentsCount,1e5,1e-4)
R <- pcaRes$T
V <- pcaRes$P
pcvars <- pcaRes$pcvar


# we only want components that cpature 95% of the total variance or a little above.
# cutIdx is the index of the compoenent, within which the variance is just equal
# to or a little greater than 95% of the total.
cutIdx <- which(cumsum(pcvars)>0.95)
if(length(cutIdx)==0){
	cutIdx <- componentsCount
}else{
	cutIdx <- cutIdx[1]
}

# slice R and V to only that number of components.
R <- R[,1:cutIdx]
V <- V[,1:cutIdx]

# the difference between experiment mean and control mean.
meanvec <- rowMeans(expm) - rowMeans(ctrl)


# all the following steps calculate shrunkMats. Refer to the ChrDir paper for detail.
# ShrunkenMats are the covariance matrix that is placed as denominator 
# in LDA formula. Notice the shrunkMats here is in the subspace of those components
# that capture about 95% of total variance.
Dd <- t(R)%*%R/samplesCount
Dd <- diag(diag(Dd))
sigma <- mean(diag(Dd))
shrunkMats <- r*Dd + sigma*(1-r)*diag(dim(R)[2])

# The LDA formula.
#  V%*%solve(shrunkMats)%*%t(V) transforms the covariance matrix from the subspace to full space.
b <- V%*%solve(shrunkMats)%*%t(V)%*%meanvec

# normlize b to unit vector
b <- b*as.vector(sqrt(1/t(b)%*%b))

# sort b to by its components' absolute value in decreasing order and get the 
# sort index
sortRes <- sort(abs(b),decreasing=TRUE,index.return=TRUE)

# sort b by the sort index
bSorted <- as.matrix(b[sortRes$ix])
# sort genes by the sort index
genesSorted <- genes[sortRes$ix]
# assign genesSorted as the row names of bSorted
rownames(bSorted) <- genesSorted

# return bSorted
bSorted <- bSorted
}

