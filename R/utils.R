# some custom functions that aren't in any other packages

loadLoomExpDiffMap_UMAP<-function# load a loom object containing expression data
(path,
  cellNameCol='obs_names',
  xname='leiden',
  has_dpt_groups=TRUE,
  has_dm = TRUE
  ){
  lfile <- connect(filename = path, skip.validate=TRUE)
  geneNames<-lfile[["row_attrs"]][["var_names"]][]
  cellNames<-lfile[["col_attrs"]][["obs_names"]][]
  expMat<- t(lfile[["matrix"]][1:length(cellNames),])
  rownames(expMat)<-geneNames
  colnames(expMat)<-cellNames
  
  cluster <- lfile[['col_attrs']][[xname]][]

  pseudotime <- lfile[['col_attrs']][["dpt_pseudotime"]][]
  dpt_groups <- cluster
  if(has_dpt_groups){
    dpt_groups <- lfile[['col_attrs']][["dpt_groups"]][]
  }
  if(has_dm){
    dm1 <- lfile[['col_attrs']][["dm1"]][]
    dm2 <- lfile[['col_attrs']][["dm2"]][]
    dm3 <- lfile[['col_attrs']][["dm3"]][]
  }
  umap1 <- lfile[['col_attrs']][["UMAP1"]][]
  umap2 <- lfile[['col_attrs']][["UMAP2"]][]

  sampTab <- data.frame(cell_name=cellNames, cluster = cluster,  pseudotime= pseudotime, dpt_groups = dpt_groups, umap1=umap1, umap2=umap2, dm1=dm1, dm2=dm2, dm3=dm3)
  
  lfile$close_all()
  rownames(sampTab) = as.vector(sampTab$cell_name)
  list(expDat = expMat, sampTab = sampTab)
}

loadLoomExpDiffMap<-function# load a loom object containing expression data
(path,
  cellNameCol='obs_names',
  xname='leiden',
  has_dpt_groups=TRUE
  ){
  lfile <- connect(filename = path, skip.validate=TRUE)
  geneNames<-lfile[["row_attrs"]][["var_names"]][]
  cellNames<-lfile[["col_attrs"]][["obs_names"]][]
  expMat<- t(lfile[["matrix"]][1:length(cellNames),])
  rownames(expMat)<-geneNames
  colnames(expMat)<-cellNames
  
  cluster <- lfile[['col_attrs']][[xname]][]

  pseudotime <- lfile[['col_attrs']][["dpt_pseudotime"]][]
  dpt_groups <- cluster
  if(has_dpt_groups){
    dpt_groups <- lfile[['col_attrs']][["dpt_groups"]][]
  }
  dpt_order <- lfile[['col_attrs']][["dpt_order"]][]

  sampTab <- data.frame(cell_name=cellNames, cluster = cluster,  pseudotime= pseudotime, dpt_groups = dpt_groups, dpt_order = dpt_order )
  
  lfile$close_all()
  rownames(sampTab) = as.vector(sampTab$cell_name)
  list(expDat = expMat, sampTab = sampTab)
}


loadLoomExpUMAP<-function# load a loom object containing expression data
(path,
  cellNameCol='obs_names',
  xname='louvain',
  phase=TRUE
  ){
  lfile <- connect(filename = path, skip.validate=TRUE)
  geneNames<-lfile[["row_attrs"]][["var_names"]][]
  cellNames<-lfile[["col_attrs"]][["obs_names"]][]
  expMat<- t(lfile[["matrix"]][1:length(cellNames),])
  rownames(expMat)<-geneNames
  colnames(expMat)<-cellNames
  
  timepoint <- lfile[['col_attrs']][['timepoint']][]
  if(phase){
    phase <- lfile[['col_attrs']][['phase']][]
  }
  n_counts <- lfile[['col_attrs']][['n_counts']][]
  n_genes <- lfile[['col_attrs']][['n_genes']][]
  percent_mito <- lfile[['col_attrs']][['percent_mito']][]
  stage_louvain <- lfile[['col_attrs']][[xname]][]
  
  pc1 <- lfile[['col_attrs']][["pca1"]][]
  pc2 <- lfile[['col_attrs']][["pca2"]][]
  pc3 <- lfile[['col_attrs']][["pca3"]][]
  umap1 <- lfile[['col_attrs']][["umap1"]][]
  umap2 <- lfile[['col_attrs']][["umap2"]][]

  if(phase){
    sampTab <- data.frame(cell_name=cellNames, timepoint=timepoint,phase=phase, n_counts=n_counts, n_genes = n_genes, percent_mito=percent_mito, stage_louvain=stage_louvain, pc1=pc1, pc2=pc2, pc3=pc3, umap1=umap1, umap2=umap2)
  }
  else{
    sampTab <- data.frame(cell_name=cellNames, timepoint=timepoint, n_counts=n_counts, n_genes = n_genes, percent_mito=percent_mito, stage_louvain=stage_louvain, pc1=pc1, pc2=pc2,pc3=pc3, umap1=umap1, umap2=umap2)
  }
  lfile$close_all()
  list(expDat = expMat, sampTab = sampTab)
}

# for slingshot
gamFit<-function(expMat,
 genes, # genes to test
  celltime){

  genes2 <- intersect(genes, rownames(expMat))
  # could print out if any missing
  ans <- apply(expMat[genes2,names(celltime)],1,function(z){
     d <- data.frame(z=z, t=celltime)
     tmp <- gam(z ~ lo(celltime), data=d)
    p <- summary(tmp)[4][[1]][1,5]
    p
  })
  ans
}



hm_ti<-function(
  expDat,
  genes,
  grps, ## vector of cellnames -> grp label
  cRow=FALSE,
  cCol=FALSE,
  limits=c(0,10),
  toScale=FALSE,
  fontsize_row=4){
  allgenes<-rownames(expDat)
  missingGenes<-setdiff(genes, allgenes)
  if(length(missingGenes)>0){
    cat("Missing genes: ", paste0(missingGenes, collapse=","), "\n")
    genes<-intersect(genes, allgenes)
  }
  value<-expDat[genes,]
  if(toScale){
      if(class(value)[1]!='matrix'){
        value <- t(scale(Matrix::t(value)))
      }
      else{
        value <- t(scale(t(value)))
      }
    }
  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]
  groupNames<-unique(grps)
  cells<-names(grps)  
##
 ## groupNames<-myGrpSort(grps)
##

  xcol <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(groupNames))
    names(xcol) <- groupNames
    anno_colors <- list(group = xcol)
    xx<-data.frame(group=as.factor(grps))
    rownames(xx)<-cells
   val_col <- colorRampPalette(rev(brewer.pal(n = 12,name = "Spectral")))(25)
   #val_col <- colorRampPalette(brewer.pal(n = 12,name = "Spectral"))(100)
  pheatmap(value, cluster_rows = cRow, cluster_cols = cCol, color=val_col,
        show_colnames = FALSE, annotation_names_row = FALSE,
##        annotation_col = annTab,
        annotation_col = xx,
        annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row=fontsize_row)
}