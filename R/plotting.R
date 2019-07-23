# plotting functions


#' plots results of findDynGenes
#'
#' heatmap
#'
#'
#' @param expDat expression matrix
#' @param dynRes result of running findDynGenes
#' @param topX top x genes to plot
#' '
#' @return heatmaps
#' 
#' @export
#'
hm_dyn_clust<-function(
  expDat,
  dynRes,
  geneAnn,
  limits=c(0,10),
  toScale=FALSE,
  fontsize_row=4){

  allgenes<-rownames(expDat)

  sampTab = dynRes$cells
  t1 = sampTab$pseudotime
  names(t1) = as.vector(sampTab$cell_name)
  grps = as.vector(sampTab$group)
  names(grps) = as.vector(sampTab$cell_name)

  ord1 = sort(t1)

  expDat = expDat[,names(ord1)]
  grps = grps[names(ord1)]


  genes = rownames(geneAnn)

  ### genes = names(genes[genes<threshold])
  # add error check if none pass
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

  xcol <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(groupNames))
    names(xcol) <- groupNames
    anno_colors <- list(group = xcol)
    xx<-data.frame(group=as.factor(grps))
    rownames(xx)<-cells
   val_col <- colorRampPalette(rev(brewer.pal(n = 12,name = "Spectral")))(25)
   #val_col <- colorRampPalette(brewer.pal(n = 12,name = "Spectral"))(100)
   
      pheatmap(value, cluster_cols = FALSE, cluster_rows= FALSE, color=val_col,
        show_colnames = FALSE, annotation_row = geneAnn,
        annotation_col = xx,
        annotation_names_col = FALSE, 
        annotation_names_row = FALSE, 
        annotation_colors = anno_colors, 
        fontsize_row=fontsize_row)
}




#' plots results of findDynGenes
#'
#' heatmap
#'
#'
#' @param expDat expression matrix
#' @param dynRes result of running findDynGenes
#' @param topX top x genes to plot
#' '
#' @return heatmaps
#' 
#' @export
#'
hm_dyn<-function(
  expDat,
  dynRes,
  topX = 25,
  cRow=FALSE,
  cCol=FALSE,
  limits=c(0,10),
  toScale=FALSE,
  fontsize_row=4,
  geneAnn = FALSE){

  allgenes<-rownames(expDat)

  sampTab = dynRes$cells
  t1 = sampTab$pseudotime
  names(t1) = as.vector(sampTab$cell_name)
  grps = as.vector(sampTab$group)
  names(grps) = as.vector(sampTab$cell_name)

  ord1 = sort(t1)

  expDat = expDat[,names(ord1)]
  grps = grps[names(ord1)]


  genes = dynRes$genes
  genes = names(sort(genes, decreasing = FALSE))[1:topX]


  ### genes = names(genes[genes<threshold])
  # add error check if none pass
  missingGenes<-setdiff(genes, allgenes)
  if(length(missingGenes)>0){
    cat("Missing genes: ", paste0(missingGenes, collapse=","), "\n")
    genes<-intersect(genes, allgenes)
  }

  # by default, order the genes by time of peak expression
  peakTime = apply(expDat[genes,], 1, which.max)
  genesOrdered = names(sort(peakTime))

  value<-expDat[genesOrdered,]
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
   if(is.data.frame(geneAnn)){
      cat("right spot\n")
      pheatmap(value, cluster_rows = cRow, cluster_cols = cCol, color=val_col,
        show_colnames = FALSE, annotation_row = geneAnn,
        annotation_col = xx,
        annotation_names_col = FALSE, 
        annotation_colors = anno_colors, 
        fontsize_row=fontsize_row)
    }
    else{
  pheatmap(value, cluster_rows = cRow, cluster_cols = cCol, color=val_col,
        show_colnames = FALSE, annotation_names_row = FALSE,
##        annotation_col = annTab,
        annotation_col = xx,
        annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row=fontsize_row)
}
}





#' @export
plotRegulon <- function(
  expMat,
  tf,
  targets,
  ntrim=0
){

  if(ntrim>0){
    str = ntrim
    stp = nrow(expMat)-ntrim

    expMat = expMat[str:stp,]
  }

  pt = 1:nrow(expMat)
 # xdat = data.frame(pt = pt, expr = c(expMat[,tf], expMat[,target]), type = c(rep("TF", nrow(expMat)), rep("Target", nrow(expMat))))
  gnames = c(tf, targets)
  xdat = gather_(as.data.frame(expMat[,gnames]), "gene", "expression", gnames)
  xdat = cbind(xdat, pt=rep(pt, length(gnames)))

  mygreys = colorRampPalette( brewer.pal(8, "Greys"))(length(targets)+4)[4:(length(targets)+4)]
  mycolors = c("#fc4e2a", mygreys)

  ggplot(xdat, aes(x=pt, y=expression)) + geom_line(aes(colour=gene, group=gene)) + geom_point(aes(colour=gene, group=gene), size=.5) + theme_bw() + scale_colour_manual(values=mycolors)


}


#' @export
plotPairs <- function(
  expMat,
  tf,
  target
){
  pt = 1:nrow(expMat)
  xdat = data.frame(pt = pt, expr = c(expMat[,tf], expMat[,target]), type = c(rep("TF", nrow(expMat)), rep("Target", nrow(expMat))))
  ggplot(xdat, aes(x=pt, y=expr, col=type)) +
  geom_line() + geom_point() + theme_bw()
}




ig_NiceGraph<-function# returns a pretty graph given a grnTab and expression data
(myG, # result of running induced.subgraph(graphX, v=genes), and ig_convertSmall
# expDat,
 grnTab,
 vLabels='Regulator'
 # bold=NULL){
){
  #  myG<-ig_convertSmall(myG);
  # calculate correlations between TF->TG
  lEdges<-length(E(myG));
  myCors<-rep(0, length=lEdges);
  edgeColors<-rep("lightblue", length=lEdges);
  for(i in seq(lEdges)){
    edgeI<-V(myG)$name[get.edge(myG, i)];
    tf<-edgeI[1];
    tg<-edgeI[2];
    #cat(tf, "->", tg,":");
#    myCors[i]<-sign(cor(expDat[tf,], expDat[tg,]));
#    myCors[i]<-sign(corrMat[tf,tg])
    myCors[i]<-sign( grnTab[grnTab$TF==tf & grnTab$TG==tg,]$corr)
    if(myCors[i]>0){
      edgeColors[i]<-"red";
    }
  }
  E(myG)$color<-edgeColors;  
  E(myG)$arrow.size=.5;
  
  ## node colors
  genes<-V(myG)$name;
 ###  geneVals<-apply(expDat[genes,], 1, mean);
  # V(myG)$color<-mp_assignColors(geneVals);
  
  # node size
  outdegree<-degree(myG, mode='out');
  V(myG)$size<-ig_scaleV( outdegree, sf=10, 4);
  
  
  # node labels
  if(length(vLabels)==1){
    if(vLabels=="Target"){
      V(myG)[V(myG)$type=="Regulator"]$label<-'';
    }
    else{
      if(vLabels=="Regulator"){
        V(myG)[V(myG)$type=="Target"]$label<-'';
      }
    }
  }
  else{
    others<-setdiff(V(myG)$name, vLabels);
    xi<-match(others, V(myG)$name);
    V(myG)[xi]$label<-'';
  }
  myG
}
  
  
ig_convertSmall<-function# change igraph attributes so that it is suitable for plotting a small GRN
(iG, ###  in which regulators and targets are already specified
 rCol="#F46D43", # default color for regulator nodes
 tCol="#66C2A5",  # default color for target nodes
 vScale=1,
 exponent=2.6
){
  
  E(iG)$color<-rgb(0,0,.5,alpha=.05);
  
  V(iG)[type=='Regulator']$label<-V(iG)[type=='Regulator']$name;
  #V(iG)[type=='Target']$label<-"";
  V(iG)$size<-ig_scaleV( degree(iG), sf=vScale, minVal=4);
  V(iG)[type=='Regulator']$shape<-"square";
  V(iG)[type=='Regulator']$label.color<-"black";
  V(iG)[type=='Regulator']$label.cex<-.75;
  V(iG)[type=='Target']$label.color<-"blue";
  V(iG)[type=='Target']$label.cex<-.5;
  V(iG)[type=='Target']$shape<-"circle";
  V(iG)[type=='Target']$frame.color=NA;
  V(iG)$label.degree<- -1*pi/2;
  V(iG)$label.dist=.25;
  nRegs<-length(V(iG)[type=='Regulator']$name);
  nTargs<-length(V(iG)[type=='Target']$name);
  #rCols<-unlist(lapply(rep(rCol, nRegs), mp_makeTransparent ) );
  rCols<-rep(rCol, nRegs);
 ## tCols<-unlist(lapply(rep(tCol, nTargs), mp_makeTransparent ) );
  tCols<-rep(tCol, nTargs);
  V(iG)[type=='Regulator']$color<-rCols;
  V(iG)[type=='Target']$color<-tCols;  
  #iG$layout<-layout.fruchterman.reingold(iG, vcount(iG)^exponent);
  iG$layout<-layout.fruchterman.reingold(iG, niter=2000, area=vcount(iG)^3.5, repulserad=vcount(iG)^2.8);  
  iG;
}

 
ig_convertLarge<-function# change igraph attributes so that it is suitable for plotting a small GRN
(iG, ###  in which regulators and targets are already specified
 rCol="#F46D43", # default color for regulator nodes
 tCol="#66C2A5",  # default color for target nodes
 vScale=1,
 exponent=2.6
){
  


  V(iG)[type=='Regulator']$color = adjustcolor("#74a9cf", alpha.f=0.7) # blue
  V(iG)[type=='Target']$color = adjustcolor("black", alpha.f=0.7) 
  
  E(iG)$color<-rgb(0,0,.5,alpha=.05);
  
  V(iG)[type=='Regulator']$label<-V(iG)[type=='Regulator']$name;
  #V(iG)[type=='Target']$label<-"";
  V(iG)$size<-ig_scaleV( degree(iG), sf=vScale, minVal=1);
  #deg = degree(iG, mode="all")
  #V(net)$size <- deg*3
 # V(iG)$size<- deg*3
  V(iG)[type=='Regulator']$shape<-"square";
  V(iG)[type=='Regulator']$frame.color<-adjustcolor("black", alpha.f=0.5)
  ###V(iG)[type=='Regulator']$frame.color<-"white"#adjustcolor("black", alpha.f=0.5)
  V(iG)[type=='Regulator']$label.color<-"black";
  V(iG)[type=='Regulator']$label.cex<-.75;

  V(iG)[type=='Target']$label.color<-"blue";
  V(iG)[type=='Target']$label.cex<-.5;
  V(iG)[type=='Target']$shape<-"circle";
###  V(iG)[type=='Target']$frame.color=NA;
  V(iG)[type=='Target']$label = ''

  E(iG)$width <- 0.5
  E(iG)$color <- adjustcolor("Grey", alpha.f=0.5)
  E(iG)$arrow.mode = 0

  E(iG)[corr<0]$color <- adjustcolor("black", alpha.f=0.5)
  E(iG)[corr>0]$color <- adjustcolor("tomato", alpha.f=0.5)


  V(iG)$label.degree<- -1*pi/2;
  V(iG)$label.dist=.25;
  nRegs<-length(V(iG)[type=='Regulator']$name);
  nTargs<-length(V(iG)[type=='Target']$name);
  #rCols<-unlist(lapply(rep(rCol, nRegs), mp_makeTransparent ) );
  rCols<-rep(rCol, nRegs);
 ## tCols<-unlist(lapply(rep(tCol, nTargs), mp_makeTransparent ) );
  tCols<-rep(tCol, nTargs);
  #V(iG)[type=='Regulator']$color<-rCols;
  #V(iG)[type=='Target']$color<-tCols;  
  #iG$layout<-layout.fruchterman.reingold(iG, vcount(iG)^exponent);
  #iG$layout<-layout.fruchterman.reingold(iG, niter=2000, area=vcount(iG)^3.5, repulserad=vcount(iG)^2.8);  
  iG$layout <- layout_with_lgl(iG)
  iG;
}



ig_convertTFs<-function# change igraph attributes so that it is suitable for plotting a GRN of only regs
(iG, ###  in which regulators and targets are already specified
 rCol="#F46D43", # default color for regulator nodes
 tCol="#66C2A5",  # default color for target nodes
 vScale=1,
 exponent=2.6
){
  
  E(iG)$color<-rgb(0,0,.5,alpha=.05);
  
  iG<-delete.vertices(iG, V(iG)[type=='Target']$name);
  V(iG)[type=='Regulator']$label<-V(iG)[type=='Regulator']$name;
  #V(iG)[type=='Target']$label<-"";
  V(iG)$size<-ig_scaleV( degree(iG), sf=vScale, minVal=1);
  V(iG)[type=='Regulator']$shape<-"square";
  V(iG)[type=='Regulator']$label.color<-"black";
  V(iG)[type=='Regulator']$label.cex<-.75;
  nRegs<-length(V(iG)[type=='Regulator']$name);


  rCols<-rep(rCol, nRegs);
  V(iG)[type=='Regulator']$color<-rCols;
  iG$layout<-layout.fruchterman.reingold(iG, niter=500, area=vcount(iG)^4, repulserad=vcount(iG)^exponent);  
  iG;
}



ig_scaleV<-function# return a vector of scaled sizes for a vector of verticies
(vals, # values associated with verticies.  Can be number of sub-net nodes (members) or degree (number of edges))
 sf=5, # scaling factor, so sf=5 means that the maximum vertix will have cex=5)
 minVal=2
){
  vals<-vals-min(vals); # shift such that m
  vals<-vals/max(vals); # scale such that range vals == 0,1
  minVal+(vals*sf); # scale to sf
}


# make a graph of the TFs, top targets, selecting only top XX targets each
ig_exemplars<-function(
  grnTab,
  geneDF, 
  tfList, # named list 
  topX = 5,
  posOnly=TRUE
  ){

  if(posOnly){
    grnTab = grnTab[grnTab$corr>0,]
  }
  grnTmp = data.frame()
  tfs = unlist(tfList)
  for(tf in tfs){
    grnX = grnTab[grnTab$TF==tf,]
    if(nrow(grnX)<topX){
      grnTmp = rbind(grnTmp, grnX)
    }
    else{
      grnX = grnX[order(grnX$zscore, decreasing=TRUE),][1:topX,]
      grnTmp = rbind(grnTmp, grnX)
    }
  }

  ig_tabToIgraph(grnTmp, directed=T)

}





ig_convertMedium<-function# change igraph attributes so that it is suitable for plotting a medium GRN
(iG, ###  in which regulators and targets are already specified
 rCol="#74a9cf", # default color for regulator nodes blue
 tCol="#66C2A5",  # default color for target nodes
 vScale=1
){
  


  V(iG)[type=='Regulator']$color = adjustcolor("#74a9cf", alpha.f=0.7) # blue
  V(iG)[type=='Target']$color = adjustcolor("black", alpha.f=0.7) 
  
  E(iG)$color<-rgb(0,0,.5,alpha=.05);
  
  V(iG)[type=='Regulator']$label<-V(iG)[type=='Regulator']$name;
  V(iG)[type=='Target']$label<-"";
  V(iG)$size<-ig_scaleV( degree(iG), sf=vScale, minVal=1);
  #deg = degree(iG, mode="all")
  #V(net)$size <- deg*3
 # V(iG)$size<- deg*3
  V(iG)[type=='Regulator']$shape<-"square";
  V(iG)[type=='Regulator']$frame.color<-adjustcolor("black", alpha.f=0.5)
  ###V(iG)[type=='Regulator']$frame.color<-"white"#adjustcolor("black", alpha.f=0.5)
  V(iG)[type=='Regulator']$label.color<-"black";
  V(iG)[type=='Regulator']$label.cex<-.75;

  V(iG)[type=='Target']$label.color<-"blue";
  V(iG)[type=='Target']$label.cex<-.5;
  V(iG)[type=='Target']$shape<-"circle";
###  V(iG)[type=='Target']$frame.color=NA;
  V(iG)[type=='Target']$label <-V(iG)[type=='Target']$name;

  E(iG)$width <- 0.5

  E(iG)[corr<0]$color <- adjustcolor("black", alpha.f=0.5)
  E(iG)[corr>0]$color <- adjustcolor("tomato", alpha.f=0.5)


  V(iG)$label.degree<- -1*pi/2;
  V(iG)$label.dist=.25;
   E(iG)$arrow.mode = 0

  iG;
}







