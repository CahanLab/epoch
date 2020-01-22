

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
  #assumes cell_name is rownames of sampTab and colnames of expDat
  sampTab = sampTab[ids,]
  expDat = expDat[,ids]

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
#' @param npoints n.points parameter for kernel smoother
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
      yy = ksmooth(t1, expDat[i,], kernel="normal", bandwidth = BW, n.points=ncol(expDat))
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
#' @param caoGenesRes result of running caoGenesRes data frame with gene peakTime, peakTimeRaw and epoch name
#' @param simplify false
#' @param directed FALSE,
#' @param weights TRUE
#'
#' @return iGraph object
wn_ig_tabToIgraph<-function#
(grnTab,
 caoGenesRes,
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

  epochs  = as.vector(caoGenesRes[c(regs,targs),]$epoch)
  verticies<-data.frame(name=c(regs,targs), label=c(regs,targs), type=types, epoch=epochs);

  
  ### iG<-graph.data.frame(tmpAns,directed=directed,v=verticies);
  iG<-igraph::graph_from_data_frame(tmpAns,directed=directed,v=verticies);
  
  if(weights){
    #E(iG)$weight<-grnTab$weight;    
    ### E(iG)$weight<-grnTab$zscore;   
    E(iG)$weight<-grnTab$graphDist;  
    E(iG)$influ<-grnTab$adjWeight;  
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

findTop <- function(vect, topX=10){
  mean(order(vect, decreasing=T)[1:topX])
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
  expDat, # for alt peakTimes
  dynRes,
  k=3,
  pThresh=0.01,
  method="kmeans"){

  # Find dyn genes
  genes = dynRes$genes
  genes = names( genes[which(genes<pThresh)] )
  if(length(genes)==0){
    cat("no dynamic genes\n")
    ans = data.frame()
  } 
  else{


    value <- t(scale(t(expSmooth[genes,])))
    cat("A\n")
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
        if(method=='pam'){
          geneMods = pam(genedist, k=k, cluster.only=TRUE)
          geneMods2 = as.character(geneMods)
        }
        else{
         geneMods = cutree(geneTree, k=k)
        geneMods2 = as.character(geneMods)
      }
    }}
   
    names(geneMods2) = genes
  ###   names(geneMods2) = names(geneMods)
    

    geneMods = geneMods2

    # find max peak
    genesPeakTimes = apply(expSmooth[genes,], 1, which.max)
    cat("B\n")
    # the expdat won't have been pruned ordered
    expDat = expDat[,colnames(expSmooth)]
    ### genesPeakTimesAlt = apply(expDat[genes,], 1, which.max)
    genesPeakTimesRaw = apply(expDat[genes,], 1, findTop)
    cat("C\n")

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
      ptX_raw = genesPeakTimesRaw[cgenes]
      genesOrdered = names(sort(ptX))

      ### tmpAns<-data.frame(gene = genesOrdered, peakTime = ptX, peakTimeRaw = ptX_raw, epoch = clu)
    
      tmpAns<-data.frame(gene = genesOrdered, peakTime = ptX[genesOrdered], peakTimeRaw = ptX_raw[genesOrdered], epoch = clu)


      rownames(tmpAns)<-genesOrdered
      ans<-rbind(ans, tmpAns)
    }

  }

  # now, re-label the clusters so they make sense
  epochs = unique(as.vector(ans$epoch))
  newepochs = as.vector(ans$epoch)
  new_e = 1
  for( i in epochs ){ # this should go in order they appear
      newepochs[ which(ans$epoch == i ) ] = new_e
      new_e = new_e + 1
  }

  ans$epoch = newepochs
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
  topX=3,
  type='zscore'
  ){

  tfs = unique(as.vector(grnTab$TF))
  epochs = unique(as.vector(geneDF$epoch))

  ans = list()
  for(epoch in epochs){

    eGeneDF = geneDF[geneDF$epoch==epoch,]

    eTFs = intersect(tfs, as.vector(eGeneDF$gene))
    if(length(eTFs)>0){
      if(type=='zscore'){
        tmpScores = unlist(lapply(eTFs, scoreTF, geneDF=eGeneDF, grnTab=grnTab))
      }
      else{
        tmpScores = unlist(lapply(eTFs, scoreTFweight, geneDF=eGeneDF, grnTab=grnTab))
      }
      names(tmpScores) = eTFs
      ans[[epoch]] = names(sort(tmpScores, decreasing=T))[1:topX]
    }
    else{
      ans[[epoch]] <- c()
    }
  }
  ans
}


# compute the consistency * normalized distance between a TF and a TG
conNormDist2 <-function(
  grnTab,
  geneDF,
  tf,
  tg,
  maxDist
  ){

  PkT = geneDF[tf,]$peakTimeRaw - geneDF[tg,]$peakTimeRaw
  regDir = grnTab[grnTab$TF==tf & grnTab$TG==tg,]$corr
  conDir = 1
  
  if( (regDir < 0) & (PkT < 0) ){
    conDir = -1
  }

  if( (regDir > 0) & (PkT > 0) ){
    conDir = -1
  }
  conDir * (abs(PkT)/maxDist )
}

# compute the consistency * normalized distance between a TF and a TG
conNormDist <-function(
  grnTab,
  geneDF,
  maxDist
  ){

  tf = as.vector(grnTab[2])
  tg = as.vector(grnTab[1])

   ### cat(tf, " ", tg,"\n")
  PkT = geneDF[tf,]$peakTimeRaw - geneDF[tg,]$peakTimeRaw
  regDir = grnTab[4]
  conDir = 1
  
  if( (regDir < 0) & (PkT < 0) ){
    conDir = -1
  }

  if( (regDir > 0) & (PkT > 0) ){
    conDir = -1
  }
  conDir * (abs(PkT)/maxDist )
}



distScore <- function(vect, trxLag = 0.05, maxNeg = - 0.2){

#   -2 * (vect-trxLag)**2 + 1

  tmp = dnorm(vect, mean=trxLag, sd = trxLag*2)
  tmp = tmp/max(tmp)
  tmp[which(vect< maxNeg)] <- 0
  tmp


}

grnDistScores <-function(
  grnTab,
  geneDF,
  maxDist
  ){

  normDist = rep(0, nrow(grnTab))
#   weightAdj = rep(1, nrow(grnTab))


    normDist = apply( grnTab, 1, conNormDist, geneDF=geneDF, maxDist=maxDist)
  
  weightAdj = distScore(normDist)

  grnTab <- cbind(grnTab, normDist=normDist)
  grnTab <- cbind(grnTab, weightAdj = weightAdj)
  grnTab
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


scoreTFweight <- function(
  TF,
  geneDF,
  grnTab){

  grnX = grnTab[grnTab$TF==TF,]
  # assumes no duplocate rows
  rownames(grnX) = as.vector(grnX$TG)

  epoch = as.vector(geneDF[geneDF$gene == TF,]$epoch[1])
  epochGenes = as.vector(geneDF[geneDF$epoch == epoch,]$gene)
 ## outScore = median(grnX$zscore) + nrow(grnX)
  outScore = sum(grnX$adjWeight * grnX$corr)

  cgenes = intersect(epochGenes, rownames(grnX))
  grnIn = grnX[cgenes,]
  ###inScore = median(grnIn$zscore) + nrow(grnIn)
  sum(grnIn$adjWeight)# * grnIn$corr)

  sum(grnX$adjWeight)
}


#' reconstruct GRN ala CLR
#'
#' give a lot of information for each TF that can be used as a basis for better and more diverse way of selecting TFs
#'
#'
#' @param grnTab expects that addDistWeight has been run 
#' @param tfs vector of transcription factor names
#' @param zThresh zscore threshold default 2
#'
#' @return data frame of TF TG zscore corr
#' 
#' @export
#'


# return a df of TF, epoch, distToStart, distToPrior, distToNext, weightMean, weightTotal
evalTFs <-function(
  grnTab, # 
  geneDF){


  allTFs = unique(as.vector(grnTab$TF))

  vect_TFs = rep("", length(allTFs))
  vect_n_TFs = rep(0, length(allTFs))
  vect_epochs = rep("", length(allTFs))
  vect_dt_start = rep(0, length(allTFs))
  vect_weightTotal = rep(0, length(allTFs))
  vect_weightMean = rep(0, length(allTFs))
  vect_nTargets = rep(0, length(allTFs))
  vect_peakTimes = rep(0, length(allTFs))

  epochs = unique(geneDF$epoch)
  nTFs = list()
  for(epoch in epochs){
    eGeneDF = geneDF[geneDF$epoch==epoch,]
    nTFs[[epoch]] = length(intersect(rownames(eGeneDF), allTFs))
  }


  for(i in seq(length(allTFs))){
    tf = allTFs[i]
    vect_TFs[i] = tf
    grnX = grnTab[grnTab$TF==tf,]
    # assumes no duplocate rows
    rownames(grnX) = as.vector(grnX$TG)

    epoch = as.vector(geneDF[geneDF$gene == tf,]$epoch[1])
    vect_epochs[i] = epoch
    vect_n_TFs[i] = nTFs[[epoch]]
    vect_weightTotal[i] = sum(grnX$adjWeight)
    vect_weightMean[i]  = mean(grnX$adjWeight)
    vect_nTargets[i] = nrow(grnX)
    vect_peakTimes[i] = as.vector(geneDF[geneDF$gene == tf,]$peakTimeRaw[1])

    eGeneDF = geneDF[geneDF$epoch==epoch,]
    ##eTFs = intersect(rownames(eGeneDF),allTFs) 
    ##eTfDF = eGeneDF[eTFs,]

    ##vect_dt_start[i] = which(rownames(eTfDF)==tf)
    vect_dt_start[i] = which(rownames( eGeneDF)==tf)

  }

  ans = data.frame(TF = vect_TFs, num_TFs = vect_n_TFs, epoch=vect_epochs, distToStart = vect_dt_start, weightTotal = vect_weightTotal, weightMean = vect_weightMean, ntargets = vect_nTargets, peakTime = vect_peakTimes)
  rownames(ans) = vect_TFs
  ans
}



#' reconstruct GRN ala CLR
#'
#' reconstruct GRN ala CLR uses minet
#'
#'
#' @param expDat unsmoothed expression matrix
#' @param tfs vector of transcription factor names
#' @param method either 'pearson' or 'MI' to be used in CLR
#' @param zThresh zscore threshold default 2
#'
#' @return data frame of TF TG zscore corr
#' 
#' @export
#'
reconstructGRN <- function(
  expDat,
  tfs,
  method='pearson',
  zThresh=2){

  if (method=='MI'){
    ttDat = t(expDat)
    mim <- build.mim(as.data.frame(ttDat),estimator="mi.empirical",disc="equalwidth")
    xnet <- clr(mim)
    xcorr = cor(ttDat)
    tfsI = intersect(tfs, colnames(ttDat))

    xnet = xnet[,tfsI]
    xcorr = xcorr[,tfsI]
    cn_extractRegsDF(xnet, xcorr, rownames(expDat), zThresh)

  }else{
    ttDat = t(expDat)
    mim <- build.mim(ttDat,estimator="pearson")
    xnet <- clr(mim)
    xcorr = cor(ttDat)
    tfsI = intersect(tfs, colnames(ttDat))

    xnet = xnet[,tfsI]
    xcorr = xcorr[,tfsI]
    cn_extractRegsDF(xnet, xcorr, rownames(expDat), zThresh)
  }


}

#' reconstruct GRN ala GENIE3
#'
#' reconstruct GRN ala GENIE3 using GENIE3
#'
#'
#' @param expDat unsmoothed expression matrix
#' @param tfs vector of transcription factor names
#' @param weightThresh threshold to consider call positive default 0
#'
#' @return data frame of TF TG zscore corr
#' 
#' @export
#'
reconstructGRN_GENIE3 <- function(
  expDat,
  tfs,
  weightThresh=0){

  require(GENIE3)

  weightmat<-GENIE3(expDat,regulators=tfs)
  gnet<-t(weightmat)
  cn_extractRegsDF(gnet,gnet,rownames(expDat),0)

}


#' add dist and weight
#'
#' add dist and weight to GRN
#'
#'
#' @param grnTab grn table
#' @param geneDF from running caoGenes
#' @param ncells number of cells
#' @param maxNeg threshold for setting adjWeight to 0 
#'
#' @return data frame with normDist and adjWeight added
#' 
#' @export
#'
addDistWeight<-function(
  grnTab,
  geneDF,
  ncells, 
  maxNeg = -0.2){

    normDists  = apply(grnTab, 1, conNormDist, geneDF=geneDF, maxDist=ncells)
    distScores = distScore(normDists, maxNeg = maxNeg)
    grnTab<-cbind(grnTab, adjWeight = grnTab$zscore * distScores)
    grnTab<-cbind(grnTab, normDist = normDists)

    weightToDist = max(grnTab$adjWeight) - grnTab$adjWeight + 0.01
    grnTab<-cbind(grnTab, graphDist = weightToDist)
    grnTab
}


#' find 'best' TFs from a warpNet run
#'
#' just pick 3 TFs, the 1t, the middle, and the last, also allow replacement by n neighbor if higher weight
#' 
#' @param tfTab result of running evalTFs
#' @param nneigh  number of cancdidates per spot to look consider
#' @param weightType weightTotal or weightMean
#'
#' @return vector of TF names
#' @export
#'
pickBestTFs<-function(
  tfTab, # result of running evalTFs
  nneigh=3,
  weightType='weightTotal'
  ){

  epochs = unique(as.vector(tfTab$epoch))

  ans = vector()

  for(epoch in epochs){
    cat("EPOCH:", epoch,"\t")
    xDF = tfTab[tfTab$epoch==epoch,]
    xDF = xDF[order(xDF$distToStart), ]
    cat("(size= ",nrow(xDF),")\n")


    # boundary cases
    ntfs = nrow(xDF)
    if(ntfs <= 3){
      ans = as.vector(xDF$TF)
    }
    else{
      cat("in kmeans\n")
      groupedData = kmeans(xDF$distToStart, centers=3)
      clusters = groupedData$cluster

      cnames = unique(clusters)
      for(cname in cnames){
        cat("subcluster ",cname,"\n")
        subDF = xDF[which(clusters==cname),]
        cat("size of subckuster: ", nrow(subDF),"\n")
        ans = append(ans, as.vector(subDF[which.max(subDF$weightTotal),]$TF))
      }
    }
  }
  unique(ans)
}





#' wrapper function to run all or subset of CLR, Epoch, and GENIE3 without PT-weight. Epoch methods are CLR run on dyngenes.
#'
#' @param expDat unsmoothed expression matrix
#' @param sampTab sample table containing "cell_names" and "dpt_groups","pseudotime" if running Epoch or PT weight
#' @param tfs vector of transcription factor names
#' @param methods which methods to run. Subset of c("CLR_pearson","CLR_MI","Epoch_pearson","Epoch_MI","GENIE3","GENIE3_dyngenes")
#' @param dpt_path gross path used to compute dynamically expressed genes
#' 
#'
#' @return list of unthresholded weight or zscore matrices, smoothed expression matrix, and geneDF
#' 
#' @export
#'
reconstructGRN_across_methods<-function(expDat,sampTab,tfs,methods,dpt_path){
  tfs<-intersect(tfs,rownames(expDat))

  #Initialize
  enet<-NA
  enet_mi<-NA
  xnet<-NA
  xnet_mi<-NA
  gnet<-NA
  gnet_dyn<-NA

  if (methods=='all'){
    methods<-c("CLR_pearson","CLR_MI","Epoch_pearson","Epoch_MI","GENIE3","GENIE3_dyngenes")
  }
  message(paste0("running ",paste(methods,collapse=", ")))

  # CLR - Pearson
  if ('CLR_pearson' %in% methods){
    print("CLR")
    mim<-build.mim(t(expDat),estimator="pearson")
    xnet<-clr(mim)
    xnet<-xnet[,tfs]
  }

  # CLR - MI
  if ('CLR_MI' %in% methods){
    print("CLR-MI")
    mim<-build.mim(as.data.frame(t(expDat)),estimator="mi.empirical",disc="equalwidth")
    xnet_mi<-clr(mim)
    xnet_mi<-xnet_mi[,tfs]
  }

  # Epoch - MI
  if ('Epoch_MI' %in% methods){
    print("Epoch-MI")  
    xdyn <- findDynGenes(expDat, sampTab, dpt_path)
    ccells = xdyn$cells
    
    expSmoothed <- grnKsmooth(expDat, ccells)
    geneDF = caoGenes(expSmoothed, expDat, xdyn, k=3, pThresh=0.01, method='kmeans')
    grnDF <- reconstructGRN(expDat[rownames(geneDF),], tfs, method="MI",zThresh=0)
    grnDF <- addDistWeight(grnDF, geneDF, ncells = ncol(expSmoothed), maxNeg= -0.2)

    enet_mi<-acast(grnDF,TF~TG,value.var="zscore")
    enet_mi[is.na(enet_mi)]<-0
    enet_mi<-t(enet_mi)
  }

  # Epoch - Pearson
  if ('Epoch_pearson' %in% methods){
    print("Epoch-MI")  
    xdyn <- findDynGenes(expDat, sampTab, dpt_path)
    ccells = xdyn$cells
    
    expSmoothed <- grnKsmooth(expDat, ccells)
    geneDF = caoGenes(expSmoothed, expDat, xdyn, k=3, pThresh=0.01, method='kmeans')
    grnDF <- reconstructGRN(expDat[rownames(geneDF),], tfs, zThresh=0)
    grnDF <- addDistWeight(grnDF, geneDF, ncells = ncol(expSmoothed), maxNeg= -0.2)
    
    enet<-acast(grnDF,TF~TG,value.var="zscore")
    enet[is.na(enet)]<-0
    enet<-t(enet)
  }

  # GENIE3, all genes
  if ('GENIE3' %in% methods){
    require(GENIE3)
    print("GENIE3")
    weightmat<-GENIE3(expDat,regulators=tfs)
    gnet<-t(weightmat)
  }

  # GENIE3 on dyn genes
  if ('GENIE3_dyngenes' %in% methods){
    require(GENIE3)
    print("GENIE3-dyn")
    xdyn <- findDynGenes(expDat, sampTab, dpt_path)
    expdyn<-expDat[names(xdyn$genes[which(xdyn$genes < 0.01)]),]
    weightmat<-GENIE3(expdyn,regulators=intersect(tfs,rownames(expdyn)))
    gnet_dyn<-t(weightmat)
  }

  res<-list(enet=enet, enet_mi=enet_mi, xnet=xnet, xnet_mi=xnet_mi, gnet=gnet, gnet_dyn=gnet_dyn)
  res<-res[!is.na(res)]
  
  # one last calculation from initial state to add to res.... all of this needs to be cleaned up because it keeps getting re-run which is unecessary
  xdyn <- findDynGenes(expDat, sampTab, dpt_path)
  expSmoothed <- grnKsmooth(expDat, ccells)
  geneDF = caoGenes(expSmoothed, expDat, xdyn, k=3, pThresh=0.01, method='kmeans')
  geneDF_all = caoGenes(expSmoothed, expDat, xdyn, k=3, pThresh=1000, method='kmeans')

  res$xdyn<-xdyn
  res$geneDF_dyn<-geneDF
  res$geneDF_all<-geneDF_all

  res

}




#' crazy wrapper function to run all combinations of methods
#'
#' @param expDat unsmoothed expression matrix
#' @param sampTab sample table containing "cell_names" and "dpt_groups","pseudotime" if running Epoch or PT weight
#' @param tfs vector of transcription factor names
#' @param dpt_path gross path used to compute dynamically expressed genes
#' 
#'
#' @return list of unthresholded weight or zscore matrices
#' 
#' @export
#'
reconstructGRN_all_methods<-function(expDat,sampTab,tfs,dpt_path){
    tfs<-intersect(tfs,rownames(expDat))

    # GENIE3, all genes
    print("GENIE3")
    require(GENIE3)
    weightmat<-GENIE3(expDat,regulators=tfs)
    gnet<-t(weightmat)

    # CLR - Pearson
    print("CLR")
    mim<-build.mim(t(expDat),estimator="pearson")
    xnet<-clr(mim)
    xnet<-xnet[,tfs]

    # CLR - MI
    print("CLR-MI")
    mim<-build.mim(as.data.frame(t(expDat)),estimator="mi.empirical",disc="equalwidth")
    xnet_mi<-clr(mim)
    xnet_mi<-xnet_mi[,tfs]

    # Epoch - MI
    print("Epoch-MI")
    xdyn <- findDynGenes(expDat, sampTab, dpt_path)
    ccells = xdyn$cells
    expSmoothed <- grnKsmooth(expDat, ccells)
    geneDF = caoGenes(expSmoothed, expDat, xdyn, k=3, pThresh=0.01, method='kmeans')
    grnDF <- reconstructGRN(expDat[rownames(geneDF),], tfs, method="MI",zThresh=0)
    grnDF <- addDistWeight(grnDF, geneDF, ncells = ncol(expSmoothed), maxNeg= -0.2)
    # Just CLR on dyn genes
    enet_mi<-acast(grnDF,TF~TG,value.var="zscore")
    enet_mi[is.na(enet_mi)]<-0
    enet_mi<-t(enet_mi)
    # CLR on dyn genes + PT weighting
    enet_mi_pt<-acast(grnDF,TF~TG,value.var="adjWeight")
    enet_mi_pt[is.na(enet_mi_pt)]<-0
    enet_mi_pt<-t(enet_mi_pt)

    # Epoch - Pearson
    print("Epoch-Pearson")
    ccells = xdyn$cells
    expSmoothed <- grnKsmooth(expDat, ccells)
    geneDF = caoGenes(expSmoothed, expDat, xdyn, k=3, pThresh=0.01, method='kmeans')
    grnDF <- reconstructGRN(expDat[rownames(geneDF),], tfs, zThresh=0)
    grnDF <- addDistWeight(grnDF, geneDF, ncells = ncol(expSmoothed), maxNeg= -0.2)
    # Just CLR on dyn genes
    enet<-acast(grnDF,TF~TG,value.var="zscore")
    enet[is.na(enet)]<-0
    enet<-t(enet)
    # CLR on dyn genes + PT weighting
    enet_pt<-acast(grnDF,TF~TG,value.var="adjWeight")
    enet_pt[is.na(enet_pt)]<-0
    enet_pt<-t(enet_pt)

    # CLR pearson on all genes + PT weighting
    print("CLR-pearson-pt")
    geneDF = caoGenes(expSmoothed, expDat, xdyn, k=3, pThresh=1000, method='kmeans')
    grnDF <- reconstructGRN(expDat[rownames(geneDF),], tfs, zThresh=0)
    grnDF <- addDistWeight(grnDF, geneDF, ncells = ncol(expSmoothed), maxNeg= -0.2)
    xnet_pt<-acast(grnDF,TF~TG,value.var="adjWeight")
    xnet_pt[is.na(xnet_pt)]<-0
    xnet_pt<-t(xnet_pt)

    # CLR MI on all genes + PT weighting
    print("CLR-MI-pt")
    #ccells = xdyn$cells
    #expSmoothed <- grnKsmooth(expDat, ccells)
    geneDF = caoGenes(expSmoothed, expDat, xdyn, k=3, pThresh=1000, method='kmeans')
    grnDF <- reconstructGRN(expDat[rownames(geneDF),], tfs, method="MI",zThresh=0)
    grnDF <- addDistWeight(grnDF, geneDF, ncells = ncol(expSmoothed), maxNeg= -0.2)
    xnet_mi_pt<-acast(grnDF,TF~TG,value.var="adjWeight")
    xnet_mi_pt[is.na(xnet_mi_pt)]<-0
    xnet_mi_pt<-t(xnet_mi_pt)

    # GENIE3 on dyn genes
    print("GENIE3-dyn")
    expdyn<-expDat[names(xdyn$genes[which(xdyn$genes < 0.01)]),]
    weightmat<-GENIE3(expdyn,regulators=intersect(tfs,rownames(expdyn)))
    gnet_dyn<-t(weightmat)

    # GENIE3 on dyn genes + PT weighting
    print("GENIE3-dyn-pt")
    grnDF<-cn_extractRegsDF(gnet_dyn,gnet_dyn,rownames(expdyn),0)
    geneDF = caoGenes(expSmoothed, expDat, xdyn, k=3, pThresh=0.01, method='kmeans')
    grnDF <- addDistWeight(grnDF, geneDF, ncells = ncol(expSmoothed), maxNeg= -0.2)
    gnet_dyn_pt<-acast(grnDF,TF~TG,value.var="adjWeight")
    gnet_dyn_pt[is.na(gnet_dyn_pt)]<-0
    gnet_dyn_pt<-t(gnet_dyn_pt)

    res<-list(enet=enet, enet_mi=enet_mi, enet_mi_pt=enet_mi_pt, enet_pt=enet_pt, xnet=xnet, xnet_mi=xnet_mi, xnet_mi_pt=xnet_mi_pt, xnet_pt=xnet_pt, gnet=gnet, gnet_dyn=gnet_dyn, gnet_dyn_pt=gnet_dyn_pt)

    res
}






