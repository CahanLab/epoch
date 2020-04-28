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
# reconstructGRN_all_methods<-function(expDat,sampTab,tfs,dpt_path){
#     tfs<-intersect(tfs,rownames(expDat))

#     # GENIE3, all genes
#     print("GENIE3")
#     require(GENIE3)
#     weightmat<-GENIE3(expDat,regulators=tfs)
#     gnet<-t(weightmat)

#     # CLR - Pearson
#     print("CLR")
#     mim<-build.mim(t(expDat),estimator="pearson")
#     xnet<-clr(mim)
#     xnet<-xnet[,tfs]

#     # CLR - MI
#     print("CLR-MI")
#     mim<-build.mim(as.data.frame(t(expDat)),estimator="mi.empirical",disc="equalwidth")
#     xnet_mi<-clr(mim)
#     xnet_mi<-xnet_mi[,tfs]

#     # Epoch - MI
#     print("Epoch-MI")
#     xdyn <- findDynGenes(expDat, sampTab, dpt_path)
#     ccells = xdyn$cells
#     expSmoothed <- grnKsmooth(expDat, ccells)
#     geneDF = caoGenes(expSmoothed, expDat, xdyn, k=3, pThresh=0.01, method='kmeans')
#     grnDF <- reconstructGRN(expDat[rownames(geneDF),], tfs, method="MI",zThresh=0)
#     grnDF <- addDistWeight(grnDF, geneDF, ncells = ncol(expSmoothed), maxNeg= -0.2)
#     # Just CLR on dyn genes
#     enet_mi<-acast(grnDF,TF~TG,value.var="zscore")
#     enet_mi[is.na(enet_mi)]<-0
#     enet_mi<-t(enet_mi)
#     # CLR on dyn genes + PT weighting
#     enet_mi_pt<-acast(grnDF,TF~TG,value.var="adjWeight")
#     enet_mi_pt[is.na(enet_mi_pt)]<-0
#     enet_mi_pt<-t(enet_mi_pt)

#     # Epoch - Pearson
#     print("Epoch-Pearson")
#     ccells = xdyn$cells
#     expSmoothed <- grnKsmooth(expDat, ccells)
#     geneDF = caoGenes(expSmoothed, expDat, xdyn, k=3, pThresh=0.01, method='kmeans')
#     grnDF <- reconstructGRN(expDat[rownames(geneDF),], tfs, zThresh=0)
#     grnDF <- addDistWeight(grnDF, geneDF, ncells = ncol(expSmoothed), maxNeg= -0.2)
#     # Just CLR on dyn genes
#     enet<-acast(grnDF,TF~TG,value.var="zscore")
#     enet[is.na(enet)]<-0
#     enet<-t(enet)
#     # CLR on dyn genes + PT weighting
#     enet_pt<-acast(grnDF,TF~TG,value.var="adjWeight")
#     enet_pt[is.na(enet_pt)]<-0
#     enet_pt<-t(enet_pt)

#     # CLR pearson on all genes + PT weighting
#     print("CLR-pearson-pt")
#     geneDF = caoGenes(expSmoothed, expDat, xdyn, k=3, pThresh=1000, method='kmeans')
#     grnDF <- reconstructGRN(expDat[rownames(geneDF),], tfs, zThresh=0)
#     grnDF <- addDistWeight(grnDF, geneDF, ncells = ncol(expSmoothed), maxNeg= -0.2)
#     xnet_pt<-acast(grnDF,TF~TG,value.var="adjWeight")
#     xnet_pt[is.na(xnet_pt)]<-0
#     xnet_pt<-t(xnet_pt)

#     # CLR MI on all genes + PT weighting
#     print("CLR-MI-pt")
#     #ccells = xdyn$cells
#     #expSmoothed <- grnKsmooth(expDat, ccells)
#     geneDF = caoGenes(expSmoothed, expDat, xdyn, k=3, pThresh=1000, method='kmeans')
#     grnDF <- reconstructGRN(expDat[rownames(geneDF),], tfs, method="MI",zThresh=0)
#     grnDF <- addDistWeight(grnDF, geneDF, ncells = ncol(expSmoothed), maxNeg= -0.2)
#     xnet_mi_pt<-acast(grnDF,TF~TG,value.var="adjWeight")
#     xnet_mi_pt[is.na(xnet_mi_pt)]<-0
#     xnet_mi_pt<-t(xnet_mi_pt)

#     # GENIE3 on dyn genes
#     print("GENIE3-dyn")
#     expdyn<-expDat[names(xdyn$genes[which(xdyn$genes < 0.01)]),]
#     weightmat<-GENIE3(expdyn,regulators=intersect(tfs,rownames(expdyn)))
#     gnet_dyn<-t(weightmat)

#     # GENIE3 on dyn genes + PT weighting
#     print("GENIE3-dyn-pt")
#     grnDF<-cn_extractRegsDF(gnet_dyn,gnet_dyn,rownames(expdyn),0)
#     geneDF = caoGenes(expSmoothed, expDat, xdyn, k=3, pThresh=0.01, method='kmeans')
#     grnDF <- addDistWeight(grnDF, geneDF, ncells = ncol(expSmoothed), maxNeg= -0.2)
#     gnet_dyn_pt<-acast(grnDF,TF~TG,value.var="adjWeight")
#     gnet_dyn_pt[is.na(gnet_dyn_pt)]<-0
#     gnet_dyn_pt<-t(gnet_dyn_pt)

#     res<-list(enet=enet, enet_mi=enet_mi, enet_mi_pt=enet_mi_pt, enet_pt=enet_pt, xnet=xnet, xnet_mi=xnet_mi, xnet_mi_pt=xnet_mi_pt, xnet_pt=xnet_pt, gnet=gnet, gnet_dyn=gnet_dyn, gnet_dyn_pt=gnet_dyn_pt)

#     res
# }

# #' updated crazy wrapper function to run all combinations of methods
# #'
# #' @param expDat unsmoothed expression matrix
# #' @param sampTab sample table containing "cell_names" and "dpt_groups","pseudotime" if running Epoch or PT weight
# #' @param tfs vector of transcription factor names
# #' @param dpt_path gross path used to compute dynamically expressed genes
# #' 
# #'
# #' @return list of unthresholded weight or zscore matrices
# #' 
# #' @export
# #'
# reconstructGRN_all_methods<-function(expDat,sampTab,tfs,dpt_path){
#     tfs<-intersect(tfs,rownames(expDat))

#     # compute dynamic genes
#     xdyn <- findDynGenes(expDat, sampTab, dpt_path)
#     ccells = xdyn$cells
#     expSmoothed <- grnKsmooth(expDat, ccells)

#     #dynamic and non-dynamic geneDF
#     geneDF<-computePT(expSmoothed,expDat,sampTab,xdyn,pThresh=0.01)
#     geneDF_all<-computePT(expSmoothed,expDat,sampTab,xdyn,pThresh=1000)

#     # CLR - Pearson / PT weight
#     print("CLR pearson with and without PT weight")
#     grnDF <- reconstructGRN(expDat[rownames(geneDF_all),], tfs, zThresh=0)
#     grnDF <- weightPT(grnDF,geneDF_all,sampTab)

#     xnet_pt<-acast(grnDF,TF~TG,value.var="adjusted_weight")
#     xnet_pt[is.na(xnet_pt)]<-0
#     xnet_pt<-t(xnet_pt)

#     xnet<-acast(grnDF,TF~TG,value.var="zscore")
#     xnet[is.na(xnet)]<-0
#     xnet<-t(xnet)

#     # CLR - MI / PT weight
#     print("CLR MI with and without PT weight")
#     grnDF <- reconstructGRN(expDat[rownames(geneDF_all),], tfs, method="MI", zThresh=0)
#     grnDF <- weightPT(grnDF,geneDF_all,sampTab)

#     xnet_mi_pt<-acast(grnDF,TF~TG,value.var="adjusted_weight")
#     xnet_mi_pt[is.na(xnet_mi_pt)]<-0
#     xnet_mi_pt<-t(xnet_mi_pt)

#     xnet_mi<-acast(grnDF,TF~TG,value.var="zscore")
#     xnet_mi[is.na(xnet_mi)]<-0
#     xnet_mi<-t(xnet_mi)
    
#     # Epoch - Pearson / PT weight
#     print("Epoch pearson with and without PT weight")
#     grnDF <- reconstructGRN(expDat[rownames(geneDF),], tfs, zThresh=0)
#     grnDF <- weightPT(grnDF,geneDF,sampTab)
#     enet<-acast(grnDF,TF~TG,value.var="zscore")
#     enet[is.na(enet)]<-0
#     enet<-t(enet)

#     enet_pt<-acast(grnDF,TF~TG,value.var="adjusted_weight")
#     enet_pt[is.na(enet_pt)]<-0
#     enet_pt<-t(enet_pt)

#     # Epoch - MI / PT weight
#     print("Epoch MI with and without PT weight")
#     grnDF <- reconstructGRN(expDat[rownames(geneDF),], tfs, method = "MI", zThresh=0)
#     grnDF <- weightPT(grnDF,geneDF,sampTab)
#     enet_mi<-acast(grnDF,TF~TG,value.var="zscore")
#     enet_mi[is.na(enet_mi)]<-0
#     enet_mi<-t(enet_mi)
 
#     enet_mi_pt<-acast(grnDF,TF~TG,value.var="adjusted_weight")
#     enet_mi_pt[is.na(enet_mi_pt)]<-0
#     enet_mi_pt<-t(enet_mi_pt)

#     # GENIE3, all genes
#     print("GENIE3")
#     grnDF <- reconstructGRN_GENIE3(expDat[rownames(geneDF_all),],tfs)
#     gnet<-acast(grnDF,TF~TG,value.var="zscore")
#     gnet[is.na(gnet)]<-0
#     gnet<-t(gnet)

#     # GENIE3 - dyngenes / PT weight
#     print("GENIE3-dyn with and without PT weight")
#     grnDF <- reconstructGRN_GENIE3(expDat[rownames(geneDF),],tfs)
#     grnDF <- weightPT(grnDF,geneDF,sampTab)
#     gnet_dyn<-acast(grnDF,TF~TG,value.var="zscore")
#     gnet_dyn[is.na(gnet_dyn)]<-0
#     gnet_dyn<-t(gnet_dyn)

#     gnet_dyn_pt<-acast(grnDF,TF~TG,value.var="adjusted_weight")
#     gnet_dyn_pt[is.na(gnet_dyn_pt)]<-0
#     gnet_dyn_pt<-t(gnet_dyn_pt)

#     res<-list(enet=enet, enet_mi=enet_mi, enet_mi_pt=enet_mi_pt, enet_pt=enet_pt, xnet=xnet, xnet_mi=xnet_mi, xnet_mi_pt=xnet_mi_pt, xnet_pt=xnet_pt, gnet=gnet, gnet_dyn=gnet_dyn, gnet_dyn_pt=gnet_dyn_pt)

#     res
# }



#' updated updated crazy wrapper function to run all combinations of methods
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

    # compute dynamic genes
    xdyn <- findDynGenes(expDat, sampTab, dpt_path)
    ccells = xdyn$cells
    expSmoothed <- grnKsmooth(expDat, ccells)

    #dynamic and non-dynamic geneDF
    geneDF<-computePT(expSmoothed,expDat,sampTab,xdyn,pThresh=0.01)
    geneDF_all<-computePT(expSmoothed,expDat,sampTab,xdyn,pThresh=1000)

    # CLR - Pearson / PT weight
    print("CLR pearson with and without PT weight")
    grnDF <- reconstructGRN(expDat[rownames(geneDF_all),], tfs, zThresh=0)
    grnDF<-crossweight(grnDF,expSmoothed)

    xnet_pt<-acast(grnDF,TF~TG,value.var="weighted_score")
    xnet_pt[is.na(xnet_pt)]<-0
    xnet_pt<-t(xnet_pt)

    xnet<-acast(grnDF,TF~TG,value.var="zscore")
    xnet[is.na(xnet)]<-0
    xnet<-t(xnet)

    # CLR - MI / PT weight
    print("CLR MI with and without PT weight")
    grnDF <- reconstructGRN(expDat[rownames(geneDF_all),], tfs, method="MI", zThresh=0)
    grnDF<-crossweight(grnDF,expSmoothed)

    xnet_mi_pt<-acast(grnDF,TF~TG,value.var="weighted_score")
    xnet_mi_pt[is.na(xnet_mi_pt)]<-0
    xnet_mi_pt<-t(xnet_mi_pt)

    xnet_mi<-acast(grnDF,TF~TG,value.var="zscore")
    xnet_mi[is.na(xnet_mi)]<-0
    xnet_mi<-t(xnet_mi)
    
    # Epoch - Pearson / PT weight
    print("Epoch pearson with and without PT weight")
    grnDF <- reconstructGRN(expDat[rownames(geneDF),], tfs, zThresh=0)
    grnDF<-crossweight(grnDF,expSmoothed)
    enet<-acast(grnDF,TF~TG,value.var="zscore")
    enet[is.na(enet)]<-0
    enet<-t(enet)

    enet_pt<-acast(grnDF,TF~TG,value.var="weighted_score")
    enet_pt[is.na(enet_pt)]<-0
    enet_pt<-t(enet_pt)

    # Epoch - MI / PT weight
    print("Epoch MI with and without PT weight")
    grnDF <- reconstructGRN(expDat[rownames(geneDF),], tfs, method = "MI", zThresh=0)
    grnDF<-crossweight(grnDF,expSmoothed)
    enet_mi<-acast(grnDF,TF~TG,value.var="zscore")
    enet_mi[is.na(enet_mi)]<-0
    enet_mi<-t(enet_mi)
 
    enet_mi_pt<-acast(grnDF,TF~TG,value.var="weighted_score")
    enet_mi_pt[is.na(enet_mi_pt)]<-0
    enet_mi_pt<-t(enet_mi_pt)

    # GENIE3, all genes
    print("GENIE3")
    grnDF <- reconstructGRN_GENIE3(expDat[rownames(geneDF_all),],tfs)
    gnet<-acast(grnDF,TF~TG,value.var="zscore")
    gnet[is.na(gnet)]<-0
    gnet<-t(gnet)

    # GENIE3 - dyngenes / PT weight
    print("GENIE3-dyn with and without PT weight")
    grnDF <- reconstructGRN_GENIE3(expDat[rownames(geneDF),],tfs)
    grnDF<-crossweight(grnDF,expSmoothed)
    gnet_dyn<-acast(grnDF,TF~TG,value.var="zscore")
    gnet_dyn[is.na(gnet_dyn)]<-0
    gnet_dyn<-t(gnet_dyn)

    gnet_dyn_pt<-acast(grnDF,TF~TG,value.var="weighted_score")
    gnet_dyn_pt[is.na(gnet_dyn_pt)]<-0
    gnet_dyn_pt<-t(gnet_dyn_pt)

    res<-list(enet=enet, enet_mi=enet_mi, enet_mi_pt=enet_mi_pt, enet_pt=enet_pt, xnet=xnet, xnet_mi=xnet_mi, xnet_mi_pt=xnet_mi_pt, xnet_pt=xnet_pt, gnet=gnet, gnet_dyn=gnet_dyn, gnet_dyn_pt=gnet_dyn_pt)

    res
}


# Function for AUPR computation after running reconstructGRN_all_methods
# ATM need to source calcPRs.R... will add this in...
aupr_res2<-function(res,gs){
    
  enet<-res$enet
  enet_mi<-res$enet_mi
  enet_mi_pt<-res$enet_mi_pt
  enet_pt<-res$enet_pt
  xnet<-res$xnet
  xnet_mi<-res$xnet_mi
  xnet_mi_pt<-res$xnet_mi_pt
  xnet_pt<-res$xnet_pt
  gnet<-res$gnet
  gnet_dyn<-res$gnet_dyn
  gnet_dyn_pt<-res$gnet_dyn_pt

  gs_dyn<-gs
  dgenes<-rownames(enet)
  gs_dyn<-gs_dyn[gs_dyn$TF %in% dgenes,]
  gs_dyn<-gs_dyn[gs_dyn$TG %in% dgenes,]
  gs_dyn<-split(gs_dyn$TG,gs_dyn$TF)

  gs<-gs[gs$TF %in% rownames(xnet),]
  gs<-gs[gs$TG %in% rownames(xnet),]
  gs<-split(gs$TG,gs$TF)    

  pr_clr<-totalPR(xnet,gs,tRange=seq(0,30,by=0.001))
  pr_clr_mi<-totalPR(xnet_mi,gs,tRange=seq(0,30,by=0.001))
  pr_clrpt<-totalPR(xnet_pt,gs,tRange=seq(0,30,by=0.001))
  pr_clrpt_mi<-totalPR(xnet_mi_pt,gs,tRange=seq(0,30,by=0.001))

  pr_epoch<-totalPR(enet,gs_dyn,tRange=seq(0,30,by=0.001))
  pr_epochpt<-totalPR(enet_pt,gs_dyn,tRange=seq(0,30,by=0.001))
  pr_epochpt_mi<-totalPR(enet_mi_pt,gs_dyn,tRange=seq(0,30,by=0.001))
  pr_epoch_mi<-totalPR(enet_mi,gs_dyn,tRange=seq(0,30,by=0.001))

  pr_genie<-totalPR(gnet,gs,tRange=seq(0,30,by=0.001))
  pr_geniedyn<-totalPR(gnet_dyn,gs_dyn,tRange=seq(0,30,by=0.001))
  pr_geniedynpt<-totalPR(gnet_dyn_pt,gs_dyn,tRange=seq(0,30,by=0.001))

  pr_clr$method<-"CLR_pearson"
  pr_clr_mi$method<-"CLR_MI"
  pr_clrpt$method<-"CLR_pearson_PT"
  pr_clrpt_mi$method<-"CLR_MI_PT"
  pr_epoch$method<-"Epoch_pearson"
  pr_epochpt$method<-"Epoch_pearson_PT"
  pr_epochpt_mi$method<-"Epoch_MI_PT"
  pr_epoch_mi$method<-"Epoch_MI"
  pr_genie$method<-"GENIE3"
  pr_geniedyn$method<-"GENIE3_dyngenes"
  pr_geniedynpt$method<-"GENIE3_dyngenes_PT"

  totalpr<-rbind(pr_clr,pr_clrpt)
  totalpr<-rbind(totalpr,pr_clr_mi)
  totalpr<-rbind(totalpr,pr_clrpt_mi)
  totalpr<-rbind(totalpr,pr_epoch)
  totalpr<-rbind(totalpr,pr_epochpt)
  totalpr<-rbind(totalpr,pr_epoch_mi)
  totalpr<-rbind(totalpr,pr_epochpt_mi)
  totalpr<-rbind(totalpr,pr_genie)
  totalpr<-rbind(totalpr,pr_geniedyn)
  totalpr<-rbind(totalpr,pr_geniedynpt)

  aupr<-data.frame(method=character(),aupr=numeric(),aupr05_norm=numeric(),aupr10_norm=numeric(),aupr20_norm=numeric())
  for (method in unique(totalpr$method)){
      subset<-totalpr[totalpr$method==method,c("Precision","Recall")]
      subset<-subset[complete.cases(subset),]
      auc<-DescTools::AUC(subset$Recall,subset$Precision, method="trapezoid", subdivisions=10000)

      subset<-subset[subset$Recall<=0.2,]        #none of the parameters in DescTools::AUC are working... do this manually
      auc20<-DescTools::AUC(subset$Recall,subset$Precision, method="trapezoid", subdivisions=10000)

      subset<-subset[subset$Recall<=0.1,]
      auc10<-DescTools::AUC(subset$Recall,subset$Precision, method="trapezoid", subdivisions=10000)

      subset<-subset[subset$Recall<=0.05,]
      auc05<-DescTools::AUC(subset$Recall,subset$Precision, method="trapezoid", subdivisions=10000)
      aupr<-rbind(aupr,data.frame(method=method,aupr=auc,aupr05_norm=auc05/0.05,aupr10_norm=auc10/0.1,aupr20_norm=auc20/0.2))

  }

  res<-list(aupr=aupr, totalpr=totalpr)

}
