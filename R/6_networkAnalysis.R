# =================== Functions to analyze networks =====================

# =================== Functions to select effector targets =====================
load_SP_effectors<-function(path){

  cwd = getwd()
  setwd(path)
  cmd = paste0("ls *.tsv")
  fnames = system(cmd, intern=T)
  
  ans = list()
  for(fname in fnames){
    newname = strsplit(fname, ".", fixed=T)[[1]][1]
    cat(newname,"\n")
  #   fname = paste0(path, "/", fname)
     xtemp = read.csv(fname, sep="\t", as.is=TRUE)
     ans[[newname]] = xtemp
   }
   setwd(cwd)
   ans
}

# Updated scoring and target selection
# score_targets adds a column "mean_score" which is just copied from chip-atlas computed average column
# 				adds a column returning the maximum score for a target across samples
#				adds a column returning the percent frequency of a hit amongst samples based on >threshold
score_targets<-function(aList,threshold=50){
	# remove STRING score
	aList<-lapply(aList,function(x){x$STRING<-NULL;x})

	# set rownames if not already set
	aList<-lapply(aList,function(x){if(!is.null(x$Target_genes)){rownames(x)<-x$Target_genes;x$Target_genes<-NULL};x})

	# add mean_score column
	aList<-lapply(aList,function(x){x$mean_score<-x[,grepl("Average",colnames(x))];x})
	# remove 'average' column
	aList<-lapply(aList,function(x){x<-x[,!grepl("Average",colnames(x))]})

	# add max_score
	aList<-lapply(aList,function(x){x$max_score<-apply(x,1,max);x})

	# add percent_freq
	aList<-lapply(aList,function(x){columns<-colnames(x)[!(colnames(x) %in% c("Target_genes","STRING","mean_score","max_score"))];
									n_pos<-apply(x[,columns],1,function(y){sum(y>threshold)});
									x$percent_freq<-n_pos/length(columns);
									x})

	aList

}

# aList the result of running score_targets
# by_rank if TRUE will return top n_targets targets; if FALSE will use thresholds
# column column name used in either ranking or thresholding the targets
# n_targets the number of targets to return for each effector if by_rank=TRUE
find_targets<-function(aList,column="max_score",by_rank=FALSE,n_targets=2000,threshold=NULL){

	# note: binding scores are -log10(q-value/FDR) from MACS2

	if (by_rank){
		aList<-lapply(aList,function(x){x<-x[order(x[,column],decreasing=TRUE),];
										n_targets<-min(n_targets,nrow(x));
										x[1:n_targets,]})
	}else if (!is.null(threshold)){
		aList<-lapply(aList,function(x){x[(x[,column]>threshold),]})
	}else{
		stop("Insufficient parameters.")
	}

	res<-lapply(aList,function(x){rownames(x)})

	res

}


# =================== Functions to find and quantify average subnetwork or module expression across time =====================


# Function to look at average module (group of genes) expression
mean_module_expression<-function(expX,module_list){
    res<-data.frame(matrix(ncol=ncol(expX),nrow=length(module_list)))
    colnames(res)<-colnames(expX)
    rownames(res)<-names(module_list)
    for (m in names(module_list)){
        genes<-intersect(module_list[[m]],rownames(expX))
        exp_sub<-expX[genes,]
        if(length(genes)==1){
            res[m,]<-expX[genes,]
        }else{
            res[m,]<-colMeans(exp_sub)
        }
    }
    res
}


mean_subnetwork_expression<-function(expX,community_list){
    # get list of subnetwork names
    subnets<-c()
    for (e in names(community_list)){
        subnets<-c(subnets,paste(e,unique(as.character(community_list[[e]]$communities)),sep="_"))
    }

    # compute subnetwork expression
    res<-data.frame(matrix(ncol=ncol(expX),nrow=length(subnets)))
    colnames(res)<-colnames(expX)
    rownames(res)<-subnets
    for (n in subnets){
        e<-strsplit(n,"_")[[1]][1]
        c<-strsplit(n,"_")[[1]][2]

        genes<-community_list[[e]]$gene[as.character(community_list[[e]]$communities)==as.character(c)]
        exp_sub<-expX[as.character(genes),]
        
        if(length(genes)==1){
            res[n,]<-expX[as.character(genes),]
        }else{
            res[n,]<-colMeans(exp_sub)
        }
    }

    res
}


subnets<-function(df,tfs,weight_column,tfonly=TRUE){
    if (tfonly){
        df<-df[df$TG %in% tfs,]
    }    

    df<-df[,c("TF","TG",weight_column,"interaction")]
    colnames(df)<-c("TF","TG","weight","interaction")

    net<-graph_from_data_frame(df,directed=FALSE)

    c<-as.data.frame(as.table(membership(cluster_louvain(net))))
    colnames(c)<-c("gene","communities")

    c
}



# =================== Functions to analyze flattened networks =====================

# Flatten networks
# function for within an epoch network, or within a static network
flatten_network<-function(grnDF,communities,community_column="communities",weight_column="weighted_score"){
    # match TFs/TGs to their communities
    grnDF$TG_module<-communities[,community_column][match(grnDF$TG,communities$gene)]
    grnDF$TF_module<-communities[,community_column][match(grnDF$TF,communities$gene)]
    grnDF$TG_module[is.na(grnDF$TG_module)]<- grnDF$TF_module[is.na(grnDF$TG_module)]

    grnDF$interaction<-"activation"
    grnDF$interaction[grnDF$corr<0]<-"repression"

    flatDF<-grnDF[,c("TF_module","TG_module",weight_column,"interaction")]
    
    # make edge weight negative for repressive edges
    flatDF$edge_score<-flatDF[,weight_column]
    flatDF$edge_score[flatDF$interaction=="repression"]<-flatDF$edge_score[flatDF$interaction=="repression"]*-1
    
    # aggregate edges between modules by sum (repressive edges will just be negative)
    flatDF<-flatDF[,c("TF_module","TG_module","edge_score")]
    flatDF<-aggregate(.~TF_module + TG_module,data=flatDF,FUN=sum)

    # normalize edge weight by possible number of edges
    n_members<-data.frame(table(communities[,community_column]))
    flatDF$TG_members<-n_members$Freq[match(as.character(flatDF$TG_module),as.character(n_members$Var1))]
    flatDF$TF_members<-n_members$Freq[match(as.character(flatDF$TF_module),as.character(n_members$Var1))]        

    flatDF$normalized_edge_score<-flatDF$edge_score/(flatDF$TF_members*flatDF$TG_members)
    flatDF$edge_length<-1-flatDF$normalized_edge_score
    
    flatDF[,c("TF_module","TG_module","edge_score","normalized_edge_score","edge_length")]
    
}


# modnet the result from running flatten_network
# module the module name at the end of the paths
find_paths_to<-function(modnet,module){
    modnet_ig<-igraph::graph_from_data_frame(modnet,directed=TRUE)
    mods<-V(modnet_ig)$name
    mods<-mods[mods!=module]

    res<-distances(modnet_ig,to=module,weights=E(modnet_ig)$edge_length,mode="out")
    colnames(res)<-c("path_length")

    paths<-c()
    for (mod in rownames(res)){
        path<-igraph::shortest_paths(modnet_ig,from=mod,to=module,mode="out",weights=E(modnet_ig)$edge_length)$vpath[[1]]$name
        paths<-c(paths,paste(path,collapse="--"))
    }

    res<-cbind(res,data.frame(path=paths))

    res$path[is.infinite(res$path_length)]<-NA

    res
}


# returns rough roots in the network
# rough roots selected as those connected to most number of nodes
rough_hierarchy<-function(modnet){
    if (class(modnet)!="igraph"){
        modnet<-igraph::graph_from_data_frame(modnet,directed=TRUE)
    }

    distances<-data.frame(distances(modnet,mode="out"))
    distances$num_connected<-rowSums(distances!=Inf)

    roots<-rownames(distances)[distances$num_connected==max(distances$num_connected)]
    list(roots=roots,num_paths=max(distances$num_connected))
}




# =================== Functions to perform shortest path analyses =====================


# Function to return shortest path from 1 TF to 1 TG in static network
# grnDF a static network table
# from a TF
# to a TG
# weight_column column name in grnDF with edge weights
static_shortest_path<-function(grnDF,from,to,weight_column="weighted_score",compare_to_average=FALSE){
    require(igraph)
    # compute relative edge lengths
    grnDF$normalized_score<-grnDF[,weight_column]/max(grnDF[,weight_column])
    grnDF$edge_length<-1-grnDF$normalized_score

    # find shortest paths
    ig<-igraph::graph_from_data_frame(grnDF[,c("TF","TG","edge_length","corr")],directed=TRUE)
    path<-igraph::shortest_paths(ig,from=from,to=to,mode="out",output="both",weights=E(ig)$edge_length)

    vpath<-path$vpath[[1]]$name
    epath<-path$epath[[1]]

    if(length(epath)==0){
        stop("No path")
    }

    # compute distance
    distance<-igraph::distances(ig,v=from,to=to,weights=E(ig)$edge_length,mode="out")[1,1]

    if (compare_to_average){
    	avg_path_length<-igraph::mean_distance(ig,directed=TRUE)
    }

    # dealing with repressive interactions
    # if path has even number of repressive interactions then TF on = target turned on
    # else, if odd, TF on = target off
    nrep<-sum(path$epath[[1]]$corr<0)
    
    if (compare_to_average){
    	list(path=vpath,distance=distance,distance_over_average=distance/avg_path_length,action=ifelse((sum(path$epath[[1]]$corr<0)%%2)==0,1,-1))
    }else{
    	list(path=vpath,distance=distance,action=ifelse((sum(path$epath[[1]]$corr<0)%%2)==0,1,-1))

    }
}


# Dynamic version of ^^
# there is probably a better way of doing this that takes in to account the dynamic aspect of the network
# for now, this function will merge networks and run shortest path on merged network..
dynamic_shortest_path<-function(grn,from,to,weight_column="weighted_score",compare_to_average=FALSE){

    # merge
    grnDF<-do.call("rbind",grn)
    grnDF<-grnDF[!duplicated(grnDF[,c("TF","TG")]),]
    
    res<-static_shortest_path(grnDF,from,to,weight_column,compare_to_average)
    res
}


# quick function to loop through ^^ for multiple TFs and targets. Returns a data frame
#
# grn dynamic network (list)
# from a vector of TFs
# to a vector of targets
dynamic_shortest_path_multiple<-function(grn,from,to,weight_column="weighted_score"){

	res<-data.frame(from=character(),to=character(),path=character(),distance=numeric(),action=numeric())

	for (tf in from){
		for (target in to){
			tryCatch({
				path<-dynamic_shortest_path(grn,from=tf,to=target,weight_column=weight_column,compare_to_average=FALSE)
				res<-rbind(res,data.frame(from=tf,to=target,path=paste(path[['path']],collapse="--"),distance=path[["distance"]],action=path[['action']]))
				},error=function(e){})
		}
	}

	grnDF<-do.call("rbind",grn)
    grnDF<-grnDF[!duplicated(grnDF[,c("TF","TG")]),]
    grnDF$normalized_score<-grnDF[,weight_column]/max(grnDF[,weight_column])
    grnDF$edge_length<-1-grnDF$normalized_score
    ig<-igraph::graph_from_data_frame(grnDF[,c("TF","TG","edge_length","corr")],directed=TRUE)

	avg_path_length<-igraph::mean_distance(ig,directed=TRUE)
	res$distance_over_average<-res$distance/avg_path_length

	res
}

# adds an extra column to the result of dynamic_shortest_path_multiple that predicts overall action
# based on correaltion between "from" and "to"
# (i.e. a check on the existing odd/even number of repressive edges)
# spDF result of running dynamic_shortest_path_multiple
# expMat 
cor_and_add_action<-function(spDF,expMat){
	genes<-union(spDF$from,spDF$to)
	expMat<-expMat[genes,]

	corrs<-apply(spDF[,c("from","to")],1,function(x){cor(expMat[x[1],],expMat[x[2],])})
	spDF$action_by_corr<-corrs
	spDF$action_by_corr<-ifelse(spDF$action_by_corr<0,-1,1)

	spDF
}


# =================== Functions to perform reachability analyses =====================

# a very basic reachability function....
# max dist numeric value, otherwise "less_than_mean"
# from is a either a character or vector of characters
static_reachability<-function(grnDF,from,max_dist="less_than_mean",weight_column="weighted_score",tfs=NULL,tf_only=FALSE){
    require(igraph)
    # compute relative edge lengths
    grnDF$normalized_score<-grnDF[,weight_column]/max(grnDF[,weight_column])
    grnDF$edge_length<-1-grnDF$normalized_score

    if (tf_only){
    	if(is.null(tfs)){
    		stop("Supply TFs.")
    	}else{
    		grnDF<-grnDF[grnDF$TG %in% tfs,]
    	}
    }

    ig<-igraph::graph_from_data_frame(grnDF[,c("TF","TG","edge_length","corr")],directed=TRUE)

    if (max_dist=="less_than_mean"){
    	avg_path_length<-igraph::mean_distance(ig,directed=TRUE)
    	max_dist<-avg_path_length
    }

    from<-intersect(from,union(grnDF$TF,grnDF$TG))

    # compute distances
    distance_mat<-igraph::distances(ig,v=from,weights=E(ig)$edge_length,mode="out")

    reachable<-lapply(1:nrow(distance_mat),function(x) {t<-colnames(distance_mat)[(distance_mat[x,]<max_dist)];
    													t<-t[t!=rownames(distance_mat)[x]]; t})
    names(reachable)<-rownames(distance_mat)

    reachable
  
}

# the dynamic version of ^^ .. again, initial simplified version
dynamic_reachability<-function(grn,from,max_dist="less_than_mean",weight_column="weighted_score",tfs=NULL,tf_only=FALSE){

    # merge
    grnDF<-do.call("rbind",grn)
    grnDF<-grnDF[!duplicated(grnDF[,c("TF","TG")]),]
    
    res<-static_reachability(grnDF,from=from,max_dist=max_dist,weight_column=weight_column,tfs=tfs,tf_only=tf_only)
    res
}




# =================== Useful plotting functions =====================
# function to plot heatmap with pre-split expX
# expList list of expression matrices
heatmap_by_treatment_group<-function(expList,sampTab,pseudotime_column="pseudotime",toScale=T,limits=c(0,5),smooth=TRUE,order_by="WAG",thresh_on=0.02,fontSize=8,anno_colors=NULL){
    if (class(expList)!="list"){
        expList<-list(A=expList)
    }
    expList<-lapply(expList,function(x){st<-sampTab[colnames(x),];
                                        st<-st[order(st[,pseudotime_column],decreasing=FALSE),];
                                        x[,rownames(st)]})
    if (smooth){
        expList<-lapply(expList,function(x){st<-data.frame(cell_name=rownames(sampTab[colnames(x),]),pseudotime=sampTab[colnames(x),pseudotime_column]);
                                            rownames(st)<-st$cell_name;
                                            grnKsmooth(x,st,BW=0.05)})
    }


    expList<-lapply(expList,function(x){x[rowSums(is.na(x))==0,]})

    # For this function, order rows by time of hitting some expression threshold
    peakTime = apply(expList[[order_by]], 1, function(x){ifelse(any(x>thresh_on),mean(which(x>thresh_on)[1:10]),length(x)+1)})
    peakTime[is.na(peakTime)]<-ncol(expList[[order_by]])+1

    #peakTime = apply(expList[[order_by]], 1, which.max)
    genesOrdered = names(sort(peakTime))    

    expX<-do.call(cbind,expList)
    sampTab<-sampTab[colnames(expX),]

    sampTab$cell_name<-rownames(sampTab)
    if ("epoch" %in% colnames(sampTab)){
        col_ann<-sampTab[,c("treatment",pseudotime_column,"epoch")]
    }else{
        col_ann<-sampTab[,c("treatment",pseudotime_column)]
    }

    value<-expX[genesOrdered,]
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

    val_col<-colorRampPalette(rev(brewer.pal(n = 12,name = "Spectral")))(25)

    # find gaps rows and columns
    gaps<-sapply(expList,ncol)
    gaps<-cumsum(gaps)
    gaps<-gaps[-length(gaps)]

    if(is.null(anno_colors)){
    	pheatmap(value,cluster_rows=F,cluster_cols=F, show_colnames=F,annotation_col=col_ann,border_color=NA,gaps_col=gaps,annotation_names_row=F,,fontsize=fontSize)
    }else{
    	pheatmap(value,cluster_rows=F,cluster_cols=F, show_colnames=F,annotation_colors = anno_colors,annotation_col=col_ann,border_color=NA,gaps_col=gaps,annotation_names_row=F,,fontsize=fontSize)

    }


}


plot_heatmap_by_treatment<-function(expList,sampTab,pseudotime_column="latent_time",toScale=T,limits=c(0,5),smooth=TRUE,anno_colors=NULL,fontSize=8){
    if (class(expList)!="list"){
        expList<-list(A=expList)
    }

    expList<-lapply(expList,function(x){st<-sampTab[colnames(x),];
                                        st<-st[order(st[,pseudotime_column],decreasing=FALSE),];
                                        x[,rownames(st)]})
    if (smooth){
        expList<-lapply(expList,function(x){st<-data.frame(cell_name=rownames(sampTab[colnames(x),]),pseudotime=sampTab[colnames(x),pseudotime_column]);
                                            rownames(st)<-st$cell_name;
                                            grnKsmooth(x,st,BW=0.03)})
    }

    # For this function, order rows by time of peak expression in WAG (first element in list)
    # chunked by community epochs
    rows<-data.frame(subnet=rownames(expList[[1]]),network=unlist(lapply(strsplit(rownames(expList[[1]]),"_"),'[[',1)))
    genesOrdered<-c()
    for (n in unique(rows$network)){
        peakTime = apply(expList[[1]][startsWith(rownames(expList[[1]]),n),],1,function(x){mean(order(x,decreasing=TRUE)[1:30])})
        genesOrdered<-c(genesOrdered,names(sort(peakTime)))
    }

    expX<-do.call(cbind,expList)
    sampTab<-sampTab[colnames(expX),]

    sampTab$cell_name<-rownames(sampTab)
    col_ann<-sampTab[,c("treatment",pseudotime_column)]

    row_ann<-data.frame(subnet=rownames(expX),network=unlist(lapply(strsplit(rownames(expX),"_"),'[[',1)))
    rownames(row_ann)<-row_ann$subnet
    row_ann$subnet<-NULL
    row_ann$network<-factor(row_ann$network,levels=unique(row_ann$network))

    value<-expX[genesOrdered,]
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

    val_col<-colorRampPalette(rev(brewer.pal(n = 12,name = "Spectral")))(25)

    # find gaps rows and columns
    gaps<-sapply(expList,ncol)
    gaps<-cumsum(gaps)
    gaps<-gaps[-length(gaps)]

    rowgaps<-table(row_ann$network)
    rowgaps<-cumsum(rowgaps)
    rowgaps<-rowgaps[-length(rowgaps)]

    if(is.null(anno_colors)){
    	pheatmap(value,cluster_rows=F,cluster_cols=F, show_colnames=F,annotation_col=col_ann,annotation_row=row_ann,border_color=NA,gaps_col=gaps,gaps_row=rowgaps,annotation_names_row=F,fontsize=fontSize)
    }else{
    	pheatmap(value,cluster_rows=F,cluster_cols=F, show_colnames=F,annotation_colors = anno_colors,annotation_col=col_ann,annotation_row=row_ann,border_color=NA,gaps_col=gaps,gaps_row=rowgaps,annotation_names_row=F,fontsize=fontSize)

    }
}





plot_reachability_community_coverage<-function(reachableList,communities){

	# get rid of duplicates in reachableList, in case there are any
	reachableList<-lapply(reachableList,function(x){x[!duplicated(x)]})

	res<-data.frame(matrix(ncol=length(reachableList),nrow=length(communities)))
	colnames(res)<-names(reachableList)
	rownames(res)<-names(communities)

	for(m in names(communities)){
		pct_covered<-sapply(reachableList,function(x){sum(x %in% communities[[m]])/length(communities[[m]])})
		res[m,]<-pct_covered
	}

	res




}






order_genes<-function(exp,sampTab,pseudotime_column,smooth=TRUE){
	st<-sampTab[colnames(exp),]
    st<-st[order(st[,pseudotime_column],decreasing=FALSE),]
    exp<-exp[,rownames(st)]

	if (smooth){
        st<-data.frame(cell_name=rownames(st),pseudotime=st[,pseudotime_column]);
        rownames(st)<-st$cell_name;
        exp<-grnKsmooth(exp,st,BW=0.05)
    }

    peakTime = apply(exp, 1, which.max)
	genesOrdered = names(sort(peakTime))
	genesOrdered

}

