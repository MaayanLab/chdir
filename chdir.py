# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 17:18:15 2015

chdir script customized for the LINCS L1000 pipeline.

Updated on Dec. 2 2015

@author: qn
"""

import numpy as np

def chdir(A, B, r=1):
#  calculate a characteristic direction for a gene expression dataset 
#  A: control numeric gene expressoion data
#  B: experiment numeric gene expression data
#  b: return value, a vector of n-components, representing the characteristic
#     direction of that gene expression dataset. n equals to the number of genes 
#     in the expression dataset.
#  r: regulaized term. A parameter that smooths the covariance matrix and reduces
#     potential noise in the dataset.

# * CHANGED *
# 1. All validation happens in a common function that validates data for all
#    differential expression methods.
# 2. Use numpy to concatenate A and B.


    ctrlCount = np.shape(A)[1]
    expmCount = np.shape(B)[1]
# place control gene expression data and experiment gene expression data into
# one matrix X.
    X = np.concatenate((A,B), axis = 1).T

# get the number of samples (rowCount), namely, the total number of replicates in   
# control and experiment. Also get the number of genes (colCount)
    (rowCount,colCount) = np.shape(X)

#  the number of output components desired from PCA. We only want to calculate
#  the chdir in a subspace that capture most variance in order to save computation 
#  workload. The number is set 20 because considering the number of genes usually 
#  present in an expression matrix 20 components would  capture most of the variance.
    if 30 > rowCount-1:
        maxComponentsNum = rowCount - 1
    else:
        maxComponentsNum = 30

# use the nipals PCA algorithm to calculate scores, loadings, and explained_var. 
# explained_var are the variances captured by each component 
    scores, loadings, explained_var = nipals(X,maxComponentsNum,1e5,1e-4)
    scores = scores.T
    loadings = loadings.T

# We only want components that cpature 95% of the total variance or a little above.
    captured_variance = 0
    for i in range(len(explained_var)):
        captured_variance += explained_var[i]
        if captured_variance > 0.999:
            break

# slice scores and loadings to only that number of components.
    scores = scores[0:i+1] # R in Neil's algorithm
    loadings = loadings[0:i+1] # V in Neil's algorithm

    scores = scores.T
    loadings = loadings.T

# the difference between experiment mean vector and control mean vector.
    meanvec = np.mean(B,axis=1,keepdims=True) - np.mean(A,axis=1,keepdims=True)

# All the following steps calculate shrunkMats. Refer to Neil's paper for detail.
# ShrunkenMats are the covariance matrix that is placed as denominator 
# in LDA formula. Notice the shrunkMats here is in the subspace of those components
# that capture about 95% of total variance.
    ctrlScores = scores[0:ctrlCount,:]
    expmScores = scores[(ctrlCount):(ctrlCount+expmCount),:]
    Dd = (np.dot(ctrlScores.T,ctrlScores)+np.dot(expmScores.T,expmScores))/(ctrlCount+expmCount-2)
#    Dd = np.diag(np.diag(Dd))
    sigma = np.mean(np.diag(Dd))
    shrunkMats = np.dot(r,Dd)+ sigma*(1-r)*np.eye(np.shape(Dd)[0])
    invMat = np.linalg.inv(shrunkMats)

# The LDA formula.
# np.dot(np.dot(loadings,shrunkMats),loadings.T) transforms the covariance 
# matrix from the subspace to full space.
    b = np.dot(loadings,np.dot(invMat,np.dot(loadings.T,meanvec)))

# normlize b to unit vector
    b /= np.linalg.norm(b)

# ! CHANGED !
# Do not sort the genes before returning them. This is because Neil's unit
# test data is not sorted.
# 
# ! CHANGED!
# Return values as Python list for consistent user interface.
    print('Done chdir')
    return b



def nipals(X,a,it=100,tol=1e-4):
    # Nipals algorithm for Principal Component Analysis
    # This function is written largely based on nipals function from R chemometrics package.

    X = np.array(X)
    (obsCount,varCount) = np.shape(X)
    Xh = X - np.tile(np.mean(X,axis=0),(obsCount,1))

    T = np.zeros((obsCount,a))
    P = np.zeros((varCount,a))
    pcvar = np.zeros(varCount)
    varTotal = np.sum(np.var(Xh,axis=0))
    currVar = varTotal
    nr = 0

    for h in range(a):
        th = np.reshape(Xh[:,0],(obsCount,-1))
        ende = False

        while not ende:
            nr = nr + 1
            ph = np.dot(Xh.T,th)/np.dot(th.T,th)
            ph = ph/np.sqrt(np.dot(ph.T,ph))
            thnew = np.dot(Xh,ph)/np.dot(ph.T,ph)
            prec = np.dot((thnew-th).T,(thnew-th))
            th = thnew

            if prec <= np.power(tol,2):
                ende = True
            if it <= nr:
                ende = True
                print('Iteration stops without convergence')

        Xh = Xh - np.dot(th,ph.T)
        T[:,h] = th[:,0]
        P[:,h] = ph[:,0]
        oldVar = currVar
        currVar = np.sum(np.var(Xh,axis=0))
        pcvar[h] = ( oldVar - currVar )/varTotal
        nr = 0

    return T, P, pcvar