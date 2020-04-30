# plotting functions

#' quick plot of dynamic networks
#'
#'
#' @param grn the result of running epochGRN
#' @param epochs the result of running assign_epochs
#' @param tfs tfs
#' @param only_TFs plot only TF network
#' 
#' @return 
#' 
#' @export
#'
plot_dynamic_network<-function(grn,epochs,tfs,only_TFs=TRUE,order=NULL){

  g<-list()

  if (!is.null(order)){
    grn<-grn[order]
  }

  for (i in 1:length(grn)){
    df<-grn[[i]]
    
    if (only_TFs){
      df<-df[df$TG %in% tfs,]
    }

    net<-graph_from_data_frame(df[,c("TF","TG")],directed=FALSE)
    tfnet<-ggnetwork(net,layout="fruchtermanreingold",cell.jitter=0)
    
    g[[i]]<-ggplot()+
      geom_edges(data=tfnet,aes(x=x, y=y, xend=xend, yend=yend),size=0.75,curvature=0.1, alpha=.6)+
      geom_nodes(data=tfnet,aes(x=x, y=y, xend=xend, yend=yend),size=5,alpha=.5)+
      #geom_nodelabel_repel(data=tfnet,aes(x=x, y=y, label=vertex.names),size=6, color="#8856a7")+
      geom_nodelabel_repel(data=tfnet,aes(x=x, y=y, label=vertex.names),size=2.5, color="#5A8BAD")+
      theme_blank()+theme(legend.position="none")+
      ggtitle(names(grn)[i])

  }

  do.call(grid.arrange,g)

}

#' quick plot of top regulators in dynamic networks
#'
#'
#' @param grn the result of running epochGRN
#' @param gene_ranks the result of running compute_pagerank
#' @param tfs tfs
#' @param numTopTFs number of top regulators to plot
#' @param numTargets number of top targets to plot for each regulator
#' @param only_TFs plot only TF network
#' @param order which epochs or transitions to plot
#' 
#' @return 
#' 
#' @export
#'
plot_top_regulators<-function(grn,gene_ranks,tfs,numTopTFs=5, numTargets=5, only_TFs=TRUE,order=NULL){

  g<-list()

  if (!is.null(order)){
    grn<-grn[order]
  }

  for (i in 1:length(grn)){
    epoch<-names(grn)[i]
    df<-grn[[epoch]]
    
    if (only_TFs){
      df<-df[df$TG %in% tfs,]
    }

    rank<-gene_ranks[[epoch]]
    topregs<-rownames(rank[rank$is_regulator==TRUE,])[1:numTopTFs]
    df<-df[df$TF %in% topregs,]

    topdf<-data.frame(TF=character(),TG=character())
    for (reg in topregs){
      targets<-as.character(df[df$TF==reg,"TG"])
      rank_targets<-rank[targets,]
      rank_targets<-rank_targets[order(rank_targets$page_rank,decreasing=TRUE),]
      num<-numTargets
      if (numTargets>length(targets)){
        num<-length(targets)
      }

      toptargets<-rownames(rank_targets)[1:num]
      add<-df[df$TF==reg,c("TF","TG")]
      add<-add[add$TG %in% toptargets,]

      topdf<-rbind(topdf,add)
    }

    net<-graph_from_data_frame(topdf[,c("TF","TG")],directed=FALSE)
    tfnet<-ggnetwork(net,layout="fruchtermanreingold",cell.jitter=0)
    tfnet$is_regulator<-as.character(tfnet$vertex.names %in% topregs)
    
    cols<-c("TRUE"="#8C4985","FALSE"="darkgray")

    g[[i]]<-ggplot()+
      geom_edges(data=tfnet,aes(x=x, y=y, xend=xend, yend=yend),size=0.75,curvature=0.1, alpha=.6)+
      geom_nodes(data=tfnet,aes(x=x, y=y, xend=xend, yend=yend,color=is_regulator),size=6,alpha=.8)+
      scale_color_manual(values=cols)+
      #geom_nodelabel_repel(data=tfnet,aes(x=x, y=y, label=vertex.names),size=6, color="#8856a7")+
      geom_nodelabel_repel(data=tfnet,aes(x=x, y=y, label=vertex.names),size=2.5, color="#5A8BAD")+
      theme_blank()+theme(legend.position="none")+
      ggtitle(names(grn)[i])

  }

  do.call(grid.arrange,g)

}

#' quick plot of top regulators given targets in dynamic networks based on reconstruction weight
#'
#'
#' @param grn the result of running epochGRN
#' @param targets targets
#' @param weight_column column name containing reconstruction weights to use
#' @param numTopRegulators number of top regulators to plot
#' @param order which epochs or transitions to plot
#' 
#' @return 
#' 
#' @export
#'
plot_targets_with_top_regulators<-function(grn,targets,weight_column="zscore",gene_ranks=NULL,numTopRegulators=5,order=NULL){
  g<-list()

  if (!is.null(order)){
    grn<-grn[order]
  }

  for (i in 1:length(grn)){
    epoch<-names(grn)[i]
    df<-grn[[epoch]]

    #look for targets in epoch GRN
    tgs<-as.character(df[df$TG %in% targets,"TG"])
    if (length(tgs)==0){
      next
    }

    #find top regulators for each target
    edges_to_keep<-data.frame(TF=character(),TG=character())
    if(weight_column=="page_rank"){
      for (tg in tgs){
        if (is.null(gene_ranks)){
          stop("Need to supply gene_ranks.")
        }
        rank<-gene_ranks[[epoch]]
        regs_of_targets<-as.character(df[df$TG==tg,"TF"])
        rank_regs<-rank[regs_of_targets,]
        rank_regs<-rank_regs[order(rank_regs$page_rank,decreasing=TRUE),]
        top_regs<-rownames(rank_regs)[1:numTopRegulators]

        edges<-df[df$TG==tg,]
        edges<-edges[edges$TF %in% top_regs,c("TF","TG")]

        edges_to_keep<-rbind(edges_to_keep,data.frame(TF=as.character(edges$TF),TG=as.character(edges$TG)))

        #random little hack to fix stupid issue with ggnetwork...
        if(nrow(edges_to_keep)==1){
          edges_to_keep<-rbind(edges_to_keep,data.frame(TF=NA,TG=NA))
        }

      }
    }else{
      for (tg in tgs){
        edges<-df[df$TG==tg,]
        edges<-edges[order(edges[,weight_column],decreasing=TRUE),]
        edges<-edges[1:numTopRegulators,c("TF","TG")]

        edges_to_keep<-rbind(edges_to_keep,data.frame(TF=as.character(edges$TF),TG=as.character(edges$TG)))
      }
    }
    #plot 
    net<-graph_from_data_frame(edges_to_keep[,c("TF","TG")],directed=FALSE)
    tfnet<-ggnetwork(net,layout="fruchtermanreingold",cell.jitter=0)
    tfnet$type<-"tf"
    tfnet$type[tfnet$vertex.names %in% targets]<-"tg"

    tfnet<-tfnet[!(is.na(tfnet$vertex.names)),]

    cols<-c("tf"="darkgray","tg"="#F5663A")

    g[[i]]<-ggplot()+
      geom_edges(data=tfnet,aes(x=x, y=y, xend=xend, yend=yend),size=0.75,curvature=0.1, alpha=.6)+
      geom_nodes(data=tfnet,aes(x=x, y=y, xend=xend, yend=yend,shape=type,color=type),size=6,alpha=.8)+
      #scale_color_manual(values=c("#F5663A","#8C4985"),labels=c("tg","tf"))+
      scale_color_manual(values=cols)+
      geom_nodes(data=tfnet[tfnet$vertex.names %in% targets,],aes(x=x,y=y,xend=xend,yend=yend),shape=24,size=6,color="black",stroke=0.25)+
      #geom_nodelabel_repel(data=tfnet,aes(x=x, y=y, label=vertex.names),size=3, color="#8856a7")+
      geom_nodelabel_repel(data=tfnet,aes(x=x, y=y, label=vertex.names),size=2.5, color="#5A8BAD")+
      theme_blank()+theme(legend.position="none")+
      ggtitle(names(grn)[i])

  }
  g[sapply(g,is.null)]<-NULL
  do.call(grid.arrange,g)

}

#' quick plot of top regulators given targets in dynamic networks based on reconstruction weight, colored by expression and interaction type
#'
#'
#' @param grn the result of running epochGRN
#' @param targets targets
#' @param epochs result of running assign_epochs
#' @param weight_column column name containing reconstruction weights to use
#' @param numTopRegulators number of top regulators to plot
#' @param order which epochs or transitions to plot
#' @param fixed_layout whether or not to fix node positions across epoch networks
#' @param declutter if TRUE, will only label nodes with active interactions in given network
#' 
#' @return 
#' 
#' @export
#'
plot_targets_with_top_regulators_detail<-function(grn,targets,epochs,weight_column="zscore",gene_ranks=NULL,numTopRegulators=5,order=NULL,fixed_layout=TRUE,declutter=TRUE){
  g<-list()

  mean_expression<-epochs$mean_expression

  if (!is.null(order)){
    grn<-grn[order]
  }

  #----------- Extract graphs for each known target network------------
  ktgraph<-list()
  for (i in 1:length(grn)){
    epoch<-names(grn)[i]
    df<-grn[[epoch]]

    #look for targets in epoch GRN
    tgs<-as.character(df[df$TG %in% targets,"TG"])
    if (length(tgs)==0){
      next
    }

    df$interaction<-"activation"
    df$interaction[df$corr<0]<-"repression"
    #find top regulators for each target
    edges_to_keep<-data.frame(TF=character(),TG=character())
    if(weight_column=="page_rank"){
      for (tg in tgs){
        if (is.null(gene_ranks)){
          stop("Need to supply gene_ranks.")
        }
        rank<-gene_ranks[[epoch]]
        regs_of_targets<-as.character(df[df$TG==tg,"TF"])
        rank_regs<-rank[regs_of_targets,]
        rank_regs<-rank_regs[order(rank_regs$page_rank,decreasing=TRUE),]
        top_regs<-rownames(rank_regs)[1:numTopRegulators]

        edges<-df[df$TG==tg,]
        edges<-edges[edges$TF %in% top_regs,c("TF","TG","interaction")]

        edges_to_keep<-rbind(edges_to_keep,data.frame(TF=as.character(edges$TF),TG=as.character(edges$TG),interaction=as.character(edges$interaction)))

        #random little hack to fix stupid issue with ggnetwork...
        if(nrow(edges_to_keep)==1){
          edges_to_keep<-rbind(edges_to_keep,data.frame(TF=NA,TG=NA))
        }

      }
    }else{
      for (tg in tgs){
        edges<-df[df$TG==tg,]
        edges<-edges[order(edges[,weight_column],decreasing=TRUE),]
        edges<-edges[1:numTopRegulators,c("TF","TG","interaction")]

        edges_to_keep<-rbind(edges_to_keep,data.frame(TF=as.character(edges$TF),TG=as.character(edges$TG),interaction=as.character(edges$interaction)))
      }
    
    }

    ktgraph[[epoch]]<-edges_to_keep

  }

  # ---------- compute node coordinates if fixed_layout==TRUE -------------
  if (fixed_layout){
    # # using sna package
    # # aggregate epoch networks
    # agg<-dplyr::bind_rows(grn)[,c("TF","TG")]
    # agg<-agg[!duplicated(agg),]
    # agg<-as_adjacency_matrix(graph_from_data_frame(agg,directed=FALSE))
    # # assign layout-- fruchterman-reingold 
    # layout<-sna::gplot.layout.fruchtermanreingold(as.matrix(agg),layout.par = c())

    # using igraph
    # aggregate epoch networks
    agg<-dplyr::bind_rows(ktgraph)[,c("TF","TG","interaction")]
    agg<-agg[!duplicated(agg),]
    agg<-graph_from_data_frame(agg,directed=FALSE)
    agg<-delete_vertices(agg,v=V(agg)$name[is.na(V(agg)$name)])

    # assign layout -- fruchterman-reingold
    layout<-layout_with_fr(agg)
    rownames(layout)<-V(agg)$name

  }


  # -------------PLOT each network--------------
  for (i in 1:length(grn)){
    epoch<-names(grn)[i]
    edges_to_keep<-ktgraph[[epoch]]

    # Covert to igraph object
    if (fixed_layout){
      # make network
      net<-graph_from_data_frame(edges_to_keep[,c("TF","TG","interaction")],directed=FALSE)
      # add vertices (no edges) that aren't in epoch network
      addvtcs<-V(agg)$name[!(V(agg)$name %in% V(net)$name)]
      net<-add_vertices(net,length(addvtcs),attr=list(name=addvtcs))

    }else{
      net<-graph_from_data_frame(edges_to_keep[,c("TF","TG","interaction")],directed=FALSE)
    }

    # remove nodes with name NA
    net<-delete_vertices(net,v=V(net)$name[is.na(V(net)$name)])
    net<-delete_vertices(net,v=V(net)$name[V(net)$name=="NA"])      # stupid hack

    # add expression values as node attribute
    expression_from<-mean_expression[mean_expression$epoch==strsplit(epoch,split="..",fixed=TRUE)[[1]][1],]
    expression_to<-mean_expression[mean_expression$epoch==strsplit(epoch,split="..",fixed=TRUE)[[1]][2],]
    if(strsplit(epoch,split="..",fixed=TRUE)[[1]][1] == strsplit(epoch,split="..",fixed=TRUE)[[1]][2]){
      # epoch (non-transition) network
      V(net)$expression<-expression_from$mean_expression[match(V(net)$name,expression_from$gene)]
    }else{
      #transition network expression of TF from source epoch, expression of target from target epoch
      V(net)$expression<-ifelse(V(net)$name %in% edges_to_keep$TF,expression_from$mean_expression[match(V(net)$name,expression_from$gene)],expression_to$mean_expression[match(V(net)$name,expression_to$gene)])
    }

    # convert to ggnetwork object
    if (fixed_layout){
      # order layout
      layout_ordered<-layout[V(net)$name,]
      tfnet<-ggnetwork(net,layout=layout_ordered,cell.jitter=0)
      }else{
        tfnet<-ggnetwork(net,layout="fruchtermanreingold",cell.jitter=0)
      }
    
    # specify if node is known regulator
    tfnet$type<-"regulator"
    tfnet$type[tfnet$vertex.names %in% targets]<-"known target"
    tfnet$type<-factor(tfnet$type,levels=c("regulator","known target"))

    tfnet<-tfnet[!(is.na(tfnet$vertex.names)),]

    cols<-c("activation"="blue","repression"="red")

    g[[i]]<-ggplot()+
      geom_edges(data=tfnet,aes(x=x, y=y, xend=xend, yend=yend,color=interaction),size=0.75,curvature=0.1, alpha=.6)+
      geom_nodes(data=tfnet,aes(x=x, y=y, xend=xend, yend=yend,shape=type,size=expression,alpha=expression),color="black")+
      scale_color_manual(values=cols)+
      #geom_nodes(data=tfnet[tfnet$vertex.names %in% targets,],aes(x=x,y=y,xend=xend,yend=yend),shape=24,size=6,color="black",stroke=0.25)+
      #geom_nodelabel_repel(data=tfnet,aes(x=x, y=y, label=vertex.names),size=3, color="#8856a7")+
      theme_blank()+
      ggtitle(names(grn)[i])

    if(declutter){
      keep<-union(edges_to_keep$TF,edges_to_keep$TG)
      g[[i]]<-g[[i]]+geom_nodelabel_repel(data=tfnet[tfnet$vertex.names %in% keep,],aes(x=x, y=y, label=vertex.names),size=2.5, color="#5A8BAD")
    }else{
      g[[i]]<-g[[i]]+geom_nodelabel_repel(data=tfnet,aes(x=x, y=y, label=vertex.names),size=2.5, color="#5A8BAD")
    }

    common_legend<-get_legend(g[[i]])

    g[[i]]<-g[[i]]+theme(legend.position="none")

  }

  # ------- Facet Plots --------
  g[sapply(g,is.null)]<-NULL

  g$legend<-common_legend
  do.call(grid.arrange,g)

}

#helpful piece of code to extract a legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
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
hm_dyn_clust<-function(
  expDat,
  dynRes,
  geneAnn,
  row_cols,
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
#    anno_colors <- list(group = xcol)
    names(row_cols) = unique(as.vector(geneAnn$epoch))
    anno_colors <- list(group = xcol, epoch=row_cols )
    xx<-data.frame(group=as.factor(grps))
    rownames(xx)<-cells

    geneX = as.data.frame(geneAnn[,"epoch"])
    rownames(geneX) = rownames(geneAnn)
    colnames(geneX) = "epoch"

  val_col <- colorRampPalette(rev(brewer.pal(n = 12,name = "Spectral")))(25)
   
   
      pheatmap(value, cluster_cols = FALSE, cluster_rows= FALSE, color=val_col,
        show_colnames = FALSE, annotation_row = geneX,
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
        annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row=fontsize_row, border_color=NA)
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







