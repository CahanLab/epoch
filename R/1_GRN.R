

# =================== Find dynamically expressed genes =====================


#' find genes expressed dynamically
#'
#' uses slingshot approach
#'
#'
#' @param expDat properly normalized expression matrix
#' @param sampTab sample table that includes pseudotime, rownames = cell_name, and a group column
#' @param path vector of group names to include
#' @param group_column column name in sampTab annotating groups in the path
#' @param pseudotime_column column name in sampTab annotating pseudotime or latent time
#'
#' @return pvals and cell info
#' 
#' @export
#'
findDynGenes<-function(expDat, 
                        sampTab, 
                        path=NULL, 
                        group_column="dpt_groups", 
                        pseudotime_column="pseudotime"){

  sampTab$dpt_groups<-sampTab[,group_column]
  sampTab$pseudotime<-sampTab[,pseudotime_column]
  sampTab$cell_name<-rownames(sampTab)

  if (is.null(path)){
    path<-unique(sampTab[,group_column])
  }

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




# =================== Static GRN reconstruction =====================

#' reconstruct GRN ala CLR
#'
#' reconstruct GRN ala CLR uses minet
#'
#'
#' @param expDat unsmoothed expression matrix
#' @param tfs vector of transcription factor names
#' @param dgenes dynamic genes, found using findDynGenes
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
  dgenes=NULL,
  method='pearson',
  zThresh=2){

  if(!is.null(dgenes)){
    expDat<-expDat[dgenes,]
  }else{
    print("Using all genes. For a proper Epoch run, provide dgenes.")
  }

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
#' @param dgenes dynamic genes, found using findDynGenes
#' @param weightThresh threshold to consider call positive default 0
#'
#' @return data frame of TF TG zscore corr
#' 
#' @export
#'
reconstructGRN_GENIE3 <- function(
  expDat,
  tfs,
  dgenes=NULL,
  weightThresh=0){

  require(GENIE3)

  if(!is.null(dgenes)){
    expDat<-expDat[dgenes,]
  }else{
    print("Using all genes. For a proper Epoch run, provide dgenes.")
  }

  tfs<-intersect(tfs,rownames(expDat))

  weightmat<-GENIE3(expDat,regulators=tfs)
  gnet<-t(weightmat)

  xcorr = cor(t(expDat))
  xnet = xcorr[,tfs]

  cn_extractRegsDF(gnet,xnet,rownames(expDat),0)

}


#' Create GRN table
#'
#' Create GRN table
#'
#'
#' @return data frame of GRN
#'
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






