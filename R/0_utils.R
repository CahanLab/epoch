# some custom functions that aren't in any other packages


#' Loads .obs and .var data from loom file
#'
#' @param path path to loom file
#' @param col_attrs column attributes to include in sample table. By default will include all in loom file.
#' @param obs_names where in loom file observation names are stored
#' @param var_names where in loom file variable names are stored
#'
#' @return list includeing expression matrix and sample table
#' 
#' @export
#'
loadDataFromLoom<-function(path,
                    col_attrs=NULL,
                    obs_names='obs_names',
                    var_names='var_names'
  ){

  lfile <- connect(filename = path, skip.validate=TRUE)

  # set obs_names and var_names
  if ("cell_names" %in% lfile[['col_attrs']]$names){
    obs_names<-"cell_names"
  }
  if ("gene_names" %in% lfile[['row_attrs']]$names){
    var_names<-"gene_names"
  }

  # get column attributes to extract
  if (is.null(col_attrs)){
    col_attrs<-lfile[['col_attrs']]$names
  }
  if (sum(!(col_attrs %in% lfile[['col_attrs']]$names))>0){
    stop("elements in col_attrs not in metadata.")
  }

  sampTab<-lfile$get.attribute.df(attributes = col_attrs, col.names=obs_names)

  geneNames<-lfile[["row_attrs"]][[var_names]][]
  cellNames<-lfile[["col_attrs"]][[obs_names]][]
  expMat<- t(lfile[["matrix"]][,])
  rownames(expMat)<-geneNames
  colnames(expMat)<-cellNames

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


#' smooths expression
#'
#' smooths expression across cells in path
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

  BW = min(BW, (max(cells$pseudotime)-min(cells$pseudotime))/10)
  
   t1 = cells$pseudotime
   names(t1) = as.vector(cells$cell_name)
   sort(t1)
   expDat = expDat[,names(t1)]
   
   ans = matrix(0, nrow=nrow(expDat), ncol=ncol(expDat))
   for(i in seq(nrow(expDat))){
      yy = ksmooth(t1, expDat[i,], kernel="normal", bandwidth = BW, x.points=t1)
      # changed from uniform spacing (n.points=ncol(expDat)) to defining x coordinates
      ans[i,] = yy$y
   }
   rownames(ans) = rownames(expDat)
   colnames(ans) = colnames(expDat)
   ans
}


#' compiles results from findDynGenes run on separate paths
#'
#' compiles results from findDynGenes run on separate paths into a single result
#'
#'
#' @param dynReslist list of results from findDynGenes
#'
#' @return a new dynRes a.k.a. the compiled result of findDynGenes
#' 
#' @export
#'
compileDynGenes<-function(dynReslist){

  gene_df<-as.data.frame(sapply(dynReslist,function(x) x[['genes']]))
  gene_df$min_pval<-apply(gene_df,1,min)
  genes<-c(gene_df$min_pval)
  names(genes)<-rownames(gene_df)

  cell_list<-lapply(dynReslist,function(x) x[['cells']])
  cell_df<-as.data.frame(do.call(rbind,cell_list))
  cell_df<-cell_df[!duplicated(cell_df),]
  cell_df<-cell_df[order(cell_df$pseudotime,decreasing=FALSE),]

  list(genes=genes,cells=cell_df)

}

#' compile epochs on branched trajectories
#'
#' compile epochs based on matching epoch label in result of running assign_epoch
#' 
#' @param epochlist list of results of running assign_epoch on each path
#'
#' @return a compiled list of epoch genes
#' @export
#'
compile_epochs<-function(epochlist){

  compiled_epochs=list()

  for (path in 1:length(epochlist)){
    epochs <- epochlist[[path]]
    epochs$mean_expression<-NULL
    for (epoch in names(epochs)){
      if (epoch %in% names(compiled_epochs)){
        compiled_epochs[[epoch]]<-union(epochs[[epoch]],compiled_epochs[[epoch]])
      }else{
        compiled_epochs[[epoch]]<-epochs[[epoch]]
      }
    }
  }

  means<-data.frame(gene=character(),epoch=character(),mean_expression=numeric())
  for (path in 1:length(epochlist)){
    mean_exp<-epochlist[[path]]$mean_expression
    means<-rbind(means,mean_exp)
  }
  compiled_means<-aggregate(means$mean_expression,by=list(gene=means$gene,epoch=means$epoch),data=means,FUN=mean)
  colnames(compiled_means)<-c("gene","epoch","mean_expression")

  compiled_epochs[["mean_expression"]]<-compiled_means
  compiled_epochs
}








