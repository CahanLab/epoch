## example Workflow with support functions
#exprMat-expression matrix with rows as  sample columns as genes

# first we normalize and use default values for sFactor, rowSamps, etc
normExpr=normScaleData(exprMat)

# Then, we cluster data to get clusters for TI inference and use
# mclust mixed gaussian modeling with no parameters and all default params of method
cl= clusterData(normExpr)

# Then, we perform TI and get pseudotime score of each cell for ecah branch along with alignemnt 
# or % involvement in each branch\
psuedo=getPsuedoScore(normExpr, rowSamp = TRUE, lowDim = NULL, pcaComp = 3, cl)

# Next, we smooth data, usually for all genes (using apply fucntion) 
# but for this example two genes are used.We also choose to analyze expression dynamics along
# one branch at a time so we filter out cells by branch at beginning. Default bw is used for kernel.  
isBranch1= psuedo$scoreBranch1>0.95; # threshold can be different for assigning cells to branches
smoothedGene1=kernelSmoothData(pseudoScore =psuedo$pseudoBranch1[isBranch1], geneExpr = normExpr[isBranch1,"Gene1"])
smoothedGene2=kernelSmoothData(pseudoScore =psuedo$pseudoBranch1[isBranch1], geneExpr = normExpr[isBranch1,"Gene2"])

# Then, we do dynamic time warping on the two genes. Note that the output of kernel is an [x, y] matrix
# where x is psuedoscore and y is gene expression level at that time so we only input y column into DTW
# output of the function is DTW score before and after flipping the second gene.
aScore=getDtwScore(pTf = smoothedGene1$y, pTrgt = smoothedGene2)

# Lastly, we would implement thresholding function to classify Tfs and assemble network, etc
