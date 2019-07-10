# normalizes and scales data using Seurat normalization function
# rowSamp indicates whether rows are samples-TRUE or rows are genes-FALSE (Default TRUE)
# sFact is caling factor (Default 10000)
# normMeth is method of normalziation- logNormalize, CLR or RC (Default "LogNormalize")
# scale is for whether or not the user wants the data scaled (Default FALSE)
# ... used for scale method for model.use, vars.to.regress, etc 

normScaleData= function(exprMat, rowSamp=TRUE, sFact=10000, normMeth="LogNormalize",scale=TRUE, ...){
  library(Seurat)
  expr=exprMat
  if (rowSamp==TRUE){
    expr=t(expr)
  }
  norm=NormalizeData(expr, normalization.method =normMeth, scale.factor = sFact)
  if (scale){
    norm=ScaleData(norm, ...)
  }
  if (rowSamp==TRUE){
    t(as.matrix(norm))
  } else {
  norm
  }
}
# gets pseudotime score for each branch of trajectory along with cell 
# score/alignment for each branch
# exprMat is expression matrix (normalized and scaled)
# rowSamp indicates whether rows are samples-TRUE or rows are genes-FALSE (Default TRUE)
# lowDim allows user to input lower dimensional embedding of data and bypass PCA if needed (Default NULL)
# pcaComp is number of PCA components to use (Default 3)
# ... allows for passing arguments such as cluster labels and start.clus and end.clus
# if known to slingshot and aid in TI
getPsuedoScore= function (exprMat, rowSamp=TRUE, lowDim=NULL, pcaComp=3,...){
  library(slingshot)
  expr=exprMat
  if (rowSamp==FALSE){
    expr=t(expr);
  }
  if (is.null(lowDim)){
      lowDim=prcomp(expr)$x[,1:pcaComp];
  }
  plot(lowDim)
  lin=getCurves(getLineages(lowDim,...));
  lines(lin)
  res=data.frame(row.names=row.names(expr))
  for (i in 1:length(lin@curves)){
    res=cbind(res, data.frame(lin@curves[[i]]$lambda, lin@curves[[i]]$w))
  }
  colIds=rep(1:(length(colnames(res))/2), times=1, each=2)
  colLabels=rep(c("pseudoBranch", "scoreBranch"),times=length(colIds)/2, each=1)
  heading=paste(as.character(colLabels), as.character(colIds), sep = "")
  colnames(res)=heading
  res
}


# smoothesData using gaussian kernel smoothing 
# pseudoScore is the inferred pseudotime score for the branch at hand
# geneExpr is the gene Expression values for a particular gene corresponding to he pseudotime
# bw is kernel banwidth use. Should be adjusted in relation to pseudotime range (Defualt 10)
# output is a two column vector with pseudocodore and value of gene at that point after smoothing
kernelSmoothData= function(pseudoScore, geneExpr, bw=10){
  library(stats)
  res=ksmooth(pseudoScore, geneExpr, kernel="normal", bandwidth = bw);
}

# clusters the data using gaussian mixed model clustering or kmeans
# exprMat is the expression matrix
# rowSamp indicates whether rows are samples-TRUE or rows are genes-FALSE (Default TRUE)
# clusMethod is method for clustering either "mclust" or "kmeans (Default "mclust")
# pcaComp is number of PCA components to use for dim. reduction (Default 3)
# lowDim allows user to input custom lower dim. embedding of data to use for clustering
# numClus is number of centers to use if kmeans method is chosen (Default 5)
# output is a vector with n cluster ids 1..n for each cell
clusterData=function(exprMat, rowSamp=TRUE, clusMeth="mclust", pcaComp=3, lowDim=NULL, numClus=5){
  expr=exprMat
  if (rowSamp==FALSE){
    expr=t(expr);
  }
  if (is.null(lowDim)){
    lowDim=prcomp(expr)$x[,1:pcaComp];
  }
  if (clusMeth=="mclust"){
    library(mclust)
    res= Mclust(lowDim)$classification
  } else if (clusMeth=="kmeans"){
    library(stats)
    res= kmeans(lowDim, centers = numClus)$cluster
  }
  plot(lowDim, col=res);
  res;
}

# gets the DTW score and reverse score between a potential Tf/target pair
# pTf is potential Tf expression level (ordered by pseudotime score)
# pTrgt is potenital target expression level (ordered by pseudotime score)
# ... used for additional parameters specified for DTW 
getDtwScore= function(pTf, pTrgt,...){
  res=list(0,0)
  library(dtw)
  res[1]=dtw(pTf, pTrgt,...)$normalizedDistance;
  res[2]=dtw(-1*pTrgt, pTf,... )$normalizedDistance;
  res;
}
