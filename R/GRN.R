

#' find genes expressed dynamically
#'
#' uses slingshot approach
#'
#'
#' @param expDat properly normalized expression matrix
#' @param sampTab sample table that includes dpt_groups, and pseudotime, rownames = cell_id
#' @param path vector of dpt_group names to grossy order/limit cells
#' @param topX number of genes to retur
#'
#' @return top X genes that are dynamically expressedzscore matrix
#' 
#' @export
#'
findDynGenes = function(expDat, sampTab, path){

  ids = vector()
  for(grp in path){
    ids = c(ids, as.vector(sampTab[sampTab$dpt_groups==grp,]$cell_name))
  }


  t1 = sampTab$pseudotime
  names(t1) = as.vector(sampTab$cell_name)

  t1C = t1[ids]
  cat("starting gammma...\n")
  gpChr <-gamFit(expDat[,names(t1C)], rownames(expDat), t1C)
  cells = data.frame(cell_name = names(t1), pseudotime = t1, group = as.vector(sampTab$dpt_groups))
  rownames(cells)= names(t1)
  cells = cells[order(cells$pseudotime),]

  ans <- list( genes = gpChr, cells = cells)
  ans
}





#' @export
getGRNMat = function(expressionMatrix, movingInterval, delay, threshold){
  tree.phate <- phate(expressionMatrix);
  phenotype=data.frame(row.names=row.names(expressionMatrix), pseudo=row.names(expressionMatrix));
  phenotype$pseudo=tree.phate$embedding[,1];
  phenotype$pseudo2=tree.phate$embedding[,2];
  #plot(tree.phate, col = tree.data$branches) #for visual inspection 
  phenotype=phenotype[row.names(phenotype)[order(phenotype$pseudo)],];
  expressionMatrix=expressionMatrix[row.names(phenotype),];
  #smoothing by moving average
  expressionMatrixSmoothed=data.frame(expressionMatrix);
  expressionMatrixSmoothed[expressionMatrixSmoothed==0]=NA;
  for (i in nrow(expressionMatrix):movingInterval){
    expressionMatrixSmoothed[i,]=colMeans(expressionMatrixSmoothed[i:(i-movingInterval),], na.rm = TRUE)
  }
  expressionMatrixSmoothed=expressionMatrixSmoothed[movingInterval:nrow(expressionMatrixSmoothed),];
  scaled=scale(expressionMatrixSmoothed);
  scaled[is.na(scaled)]=0;
  activation=matrix(nrow=ncol(scaled), ncol=ncol(scaled));
  deactivation=matrix(nrow=ncol(scaled), ncol=ncol(scaled));
  for (i in 1:nrow(activation)){
    for (j in 1:ncol(activation)){
      tf=scaled[,i];
      target=scaled[,j]
      tf=filter.fft(tf,BW=0.01, n=10)
      target=filter.fft(target,BW=0.01, n=10)
      activation[i,j]=dtw(Re(tf)[1:(length(tf)-delay)], Re(target)[(delay+1):length(target)])$normalizedDistance;
      deactivation[i,j]=dtw(Re(tf)[1:(length(tf)-delay)], -1*Re(target)[(delay+1):length(target)])$normalizedDistance;
    }
    print(i)
  }
  blah=list(activation, deactivation)
  row.names(activation)=colnames(expressionMatrix)
  colnames(activation)=colnames(expressionMatrix)
  bool=activation;
    for (i in 1:nrow(activation)){
      for (j in 1:ncol(activation)){
        if((activation[i,j])<threshold){
          bool[i,j]=1
        }
        else if((deactivation[i,j])<threshold){
          bool[i,j]=-1
        }
        else {
          bool[i,j]=0
        }
      }
    testGRN=melt(as.matrix(bool))
    testGRN=testGRN[testGRN$value!=0,]
    testGRN$value[testGRN$value==-1]="Inhibition"
    testGRN$value[testGRN$value==1]="Activation"
    testGRN=testGRN[, c(1 ,3, 2)]
    colnames(testGRN)=c("TF", "Interaction", "Target")
    testGRN
  }
}

#' smooths expression
#'
#' smooths expression
#'
#'
#' @param expDat  expression matrix
#' @param cells sample table that includes "cell_name"  "pseudotime" "group"     
#' @param BW bandwith for kernel smoother
#'
#' @return smoothed matrix
#' 
#' @export
#'
grnKsmooth<-function(
  expDat,
  cells,
  BW = .25){

  
   t1 = cells$pseudotime
   names(t1) = as.vector(cells$cell_name)
   sort(t1)
   expDat = expDat[,names(t1)]
   
   ans = matrix(0, nrow=nrow(expDat), ncol=ncol(expDat))
   for(i in seq(nrow(expDat))){
      yy = ksmooth(t1, expDat[i,], kernel="normal", bandwidth = BW)
     ans[i,] = yy$y
   }
   rownames(ans) = rownames(expDat)
   colnames(ans) = colnames(expDat)
   ans
}


#' @export
grnSSF = function(
  expMat, # rows are cells
  cellOrder, # as provided by TI
  movingInterval = 5,
  BW = 0.01,
  fftN = 10){
    
    genes = colnames(expMat)
    expMat = expMat[cellOrder,]

    # smoothing by moving average
    expressionMatrixSmoothed=data.frame(expMat)
    colnames(expressionMatrixSmoothed) = genes
    expressionMatrixSmoothed[expressionMatrixSmoothed==0]=NA
    for (i in nrow(expressionMatrix):movingInterval){
      expressionMatrixSmoothed[i,]=colMeans(expressionMatrixSmoothed[i:(i-movingInterval),], na.rm = TRUE)
    }
    expressionMatrixSmoothed=expressionMatrixSmoothed[movingInterval:nrow(expressionMatrixSmoothed),]
    scaled=scale(expressionMatrixSmoothed)
    scaled[is.na(scaled)]=0

    colnames(scaled) = genes
    cat("FFT\n")
    x = apply(scaled, 2, filter.fft, BW=BW, n = fftN)
    colnames(x) = genes
    list(scaled=scaled, ffTed = Re(x))
}

#' @export
grnDTWx = function(grnFF_res, tfs, delay=5){
 

  ncells = nrow(grnFF_res)
  targets = colnames(grnFF_res)
  ntargets = length(targets)

  tfs = intersect(tfs, targets)
  ntfs = length(tfs)

  activation   = matrix(nrow=ntargets, ncol=ntfs)
  deactivation = matrix(nrow=ntargets, ncol=ntfs)
  row.names(activation) = targets
  colnames(activation)  = tfs

  for (i in 1:ntargets){
    #  cat(targets[i])
      for (j in 1:ntfs){
    #    cat(   " ",tfs[j],"\n")
        activation[i,j]   = dtw( grnFF_res[,j][1:(ncells-delay)], grnFF_res[,i][(delay+1):ncells] )$normalizedDistance
        deactivation[i,j] = dtw( grnFF_res[,j][1:(ncells-delay)], -1*grnFF_res[,i][(delay+1):ncells] )$normalizedDistance
        
      }
  }

  list(act = activation, deact = deactivation)
}




grn_zscores<-function 
(corrVals,
 tfs
  ){
  zscs<-mat_zscores(corrVals);
  gc();
  zscs[,tfs];
}


#' compute context dependent zscores
#'
#' slightly modidied from JJ Faith et al 2007
#' @param corrMat correlation matrix
#'
#' @return matrix of clr zscores
#'
mat_zscores<-function# computes sqrt(zscore_row + zscore_col) .. 
(corrMat ### correlation matrix
){
  corrMat<-abs(corrMat);
  zscs_2<-round(scale(corrMat), 3);
  rm(corrMat);
  gc()
  zscs_2 + t(zscs_2);
}

cn_extractRegsDF<-function# 
(zscores,
 corrMatrix,
 genes,
 threshold,
 removeAuto=TRUE
){
    
  targets<-vector();
  regulators=vector();
  zscoresX<-vector();
  correlations<-vector();
  
  targets<-rep('', 1e6);
  regulators<-rep('', 1e6);
  zscoresX<-rep(0, 1e6);
  correlations<-rep(0, 1e6);
  
  str<-1;
  stp<-1;
  for(target in genes){
    x<-zscores[target,];
    regs<-names(which(x>threshold));
    if(length(regs)>0){
      zzs<-x[regs];
      corrs<-corrMatrix[target,regs];
      ncount<-length(regs);
      stp<-str+ncount-1;
      targets[str:stp]<-rep(target, ncount);
      #    targets<-append(targets,rep(target, ncount));
      regulators[str:stp]<-regs;
      #regulators<-append(regulators, regs);
      #    zscoresX<-append(zscoresX, zzs);
      zscoresX[str:stp]<-zzs;
      correlations[str:stp]<-corrs;
      str<-stp+1;
    }
    #    correlations<-append(correlations, corrs);
  }
  targets<-targets[1:stp];
  regulators<-regulators[1:stp];
  zscoresX<-zscoresX[1:stp];
  correlations<-correlations[1:stp];
  
  
  grn = data.frame(target=targets, reg=regulators, zscore=zscoresX, corr=correlations)
  colnames(grn)[1:2]<-c("TG", "TF");
  
  if(removeAuto){
    tfx = as.vector(grn$TF)
    tgx = as.vector(grn$TG)
    xi = which (tfx != tgx)
    grn = grn[xi,]
  }
  grn

}



#' convert a table to an igraph
#'
#' convert a table to an igraph. This adds an nEnts vertex attribute to count the number of entities in the sub-net
#' 
#' @param grnTab table of TF, TF, maybe zscores, maybe correlations
#' @param simplify false
#' @param directed FALSE,
#' @param weights TRUE
#'
#' @return iGraph object
ig_tabToIgraph<-function#
(grnTab,
 simplify=FALSE,
 directed=FALSE,
 weights=TRUE
){
  
  tmpAns<-as.matrix(grnTab[,c("TF", "TG")]);
  regs<-as.vector(unique(grnTab[,"TF"]));
  ###cat("Length TFs:", length(regs), "\n");
  targs<-setdiff( as.vector(grnTab[,"TG"]), regs);
  
###  cat("Length TGs:", length(targs), "\n");
  myRegs<-rep("Regulator", length=length(regs));
  myTargs<-rep("Target", length=length(targs));
  
  types<-c(myRegs, myTargs);
  verticies<-data.frame(name=c(regs,targs), label=c(regs,targs), type=types);
  
  ### iG<-graph.data.frame(tmpAns,directed=directed,v=verticies);
  iG<-igraph::graph_from_data_frame(tmpAns,directed=directed,v=verticies);
  
  if(weights){
    #E(iG)$weight<-grnTab$weight;    
    E(iG)$weight<-grnTab$zscore;    
  }
  if(directed){
    E(iG)$corr<-sign(grnTab$corr)  
  }
  
  if(simplify){
    iG<-simplify(iG);
  }
  V(iG)$nEnts<-1;
  iG;
}



#' cluster and order genes
#'
#' cluster and order genes
#' 
#' @param expSmooth expression matrix
#' @param dynRes result of running findDynGenes
#' @param pThresh pval threshold
#' @param k number of clusters
#'
#' @return data.frame of dynamically expressed genes, cluster, peakTime, ordered by peaktime
#' @export
#'
caoGenes<-function(
  expSmooth,
  dynRes,
  k=3,
  pThresh=0.01,
  method="kmeans"){

  # Find dyn genes
  genes = dynRes$genes
  genes = names( genes[genes<pThresh] )
  if(length(genes)==0){
    cat("no dynamic genes\n")
    ans = data.frame()
  } 
  else{


    value <- t(scale(t(expSmooth[genes,])))

    # cluster
    ###genedist = utils_myDist(expSmooth[genes,])


    genedist = utils_myDist(value)
    geneTree = hclust(genedist, 'ave')

    if(method=='kmeans'){
      kmeansRes = kmeans(genedist, centers=k)
      geneMods2 = as.character(kmeansRes$cluster)
    }
    else{
      if(method=='dct'){
        geneMods = cutreeDynamicTree(geneTree, deepSplit=FALSE, minModuleSize=50)
        geneMods2 = as.character(geneMods)
      }
      else{
         geneMods = cutree(geneTree, k=k)
        geneMods2 = as.character(geneMods)
      }
    }
   
    names(geneMods2) = genes
  ###   names(geneMods2) = names(geneMods)
    

    geneMods = geneMods2

    # find max peak
    genesPeakTimes = apply(expSmooth[genes,], 1, which.max)


    ### peakTime = names(sort(gemeakTime))

    # order the clusters
    clusterOrder = vector()
    clusterNames = unique(geneMods)
    meanPeakTimes = rep(0, length(clusterNames))
    names(meanPeakTimes) = clusterNames
    for( clu in clusterNames ){
      cat(clu,"\t")
      cgenes = names(which(geneMods==clu))
      cat(length(cgenes),"\t")
      meanPeakTimes[clu] = mean(genesPeakTimes[cgenes])
      cat( meanPeakTimes[clu], "\n")
    }
 
    meanPeakTimes = sort(meanPeakTimes) 
    cluOrdered = names(meanPeakTimes)
    # make the data.frame

    ans<-data.frame()
    for( clu in cluOrdered ){
      cat(clu,"\n")
      cgenes = names(which(geneMods==clu))
      ptX = genesPeakTimes[cgenes]
      genesOrdered = names(sort(ptX))

      tmpAns<-data.frame(gene = genesOrdered, peakTime = ptX, epoch = clu)
      rownames(tmpAns)<-genesOrdered
      ans<-rbind(ans, tmpAns)
    }

  }
  ans
}

#' find TFs from each epoch
#'
#' do this based on # targets, corr, and median weight
#' 
#' @param geneDF from running caoGenes
#' @param grnTab from running cn_extractRegsDF
#' @param topX how man per epoch
#'
#' @return list of vectors with TF names
#' @export
#'
pickExemplars<-function(
  geneDF,
  grnTab,
  topX=3
  ){

  tfs = unique(as.vector(grnTab$TF))
  epochs = unique(as.vector(geneDF$epoch))

  ans = list()
  for(epoch in epochs){

    eGeneDF = geneDF[geneDF$epoch==epoch,]

    eTFs = intersect(tfs, as.vector(eGeneDF$gene))
    tmpScores = unlist(lapply(eTFs, scoreTF, geneDF=eGeneDF, grnTab=grnTab))
    names(tmpScores) = eTFs
    ans[[epoch]] = names(sort(tmpScores, decreasing=T))[1:topX]
  }
  ans
}



scoreTF <- function(
  TF,
  geneDF,
  grnTab){

  grnX = grnTab[grnTab$TF==TF,]
  # assumes no duplocate rows
  rownames(grnX) = as.vector(grnX$TG)

  epoch = as.vector(geneDF[geneDF$gene == TF,]$epoch[1])
  epochGenes = as.vector(geneDF[geneDF$epoch == epoch,]$gene)
 ## outScore = median(grnX$zscore) + nrow(grnX)
  outScore = sum(grnX$zscore * grnX$corr)

  cgenes = intersect(epochGenes, rownames(grnX))
  grnIn = grnX[cgenes,]
  ###inScore = median(grnIn$zscore) + nrow(grnIn)
  inScore = sum(grnIn$zscore * grnIn$corr)

  #list(inScore=inScore, outScore=outScore, epochGenes = epochGenes, grnX= grnX)

  inScore - abs(outScore)
}


#' reconstruct GRN ala CLR
#'
#' reconstruct GRN ala CLR uses minet
#'
#'
#' @param expSm smoothed expression data
#' @param tfs vector of transcription factor names
#' @param zThresh zscore threshold default 2
#'
#' @return data frame of TF TG zscore corr
#' 
#' @export
#'
reconstructGRN <- function(
  expSm,
  tfs,
  zThresh=2){

  ttDat = t(testSm)
  mim <- build.mim(ttDat,estimator="pearson")
  xnet <- clr(mim)
  xcorr = cor(ttDat)
  tfsI = intersect(tfs, colnames(ttDat))

  xnet = xnet[,tfsI]
  xcorr = xcorr[,tfsI]
  cn_extractRegsDF(xnet, xcorr, rownames(expSm), zThresh)
}











