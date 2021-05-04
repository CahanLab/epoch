# =================== Functions to compare networks =====================


#' Function to bootstrap Epoch reconstruction 
#'
#' @param expX genes-by-cells expression matrix
#' @param sampTab sample table
#' @param tfs TFs
#' @param n_reconstructions the number of networks to reconstruct
#' @param ncells_per_sample the number of cells to sample for each reconstruction
#' @param pThresh p-value threshold for selecting dynamic genes
#' @param zThresh z-score threshold in static reconstruction
#' @param pseudotime_column column name in sampTab with pseudotime annotation
#' @param group_column column name in sampTab for restricting reconstruction
#' @param method CLR method, either "pearson" or "MI"
#' @param crossweight whether or not to perform crossweighting
#' @param limit_to vector of genes. If not NULL, will skip finding dynamic genes and use genes in limit_to for reconstruction
#'
#' @return list of grnDFs
#' 
#' @export
#'
sample_and_epoch_reconstruct<-function(expX,
								sampTab,
								tfs,
								n_reconstructions=10,
								ncells_per_sample=400,
								pThresh=0.05,
								zThresh=5,
								pseudotime_column="dpt_pseudotime",
								group_column="leiden",
								method="pearson",
								crossweight=FALSE,
								limit_to=NULL){

	res<-list()

	for (i in 1:n_reconstructions){
		print(i)
		# sample
		sampTab_subset<-sampTab[sample(rownames(sampTab),ncells_per_sample,replace=FALSE),]
		expX_subset<-expX[,rownames(sampTab_subset)]
		expX_subset<-expX_subset[rowSums(expX_subset)!=0,]
		tfs_subset<-intersect(tfs,rownames(expX_subset))

		# reconstruct
		if (is.null(limit_to)){
		  xdyn<-findDynGenes(expX_subset,sampTab_subset,unique(sampTab_subset[,group_column]),group_column=group_column,pseudotime_column=pseudotime_column)
		  dgenes<-names(xdyn$genes)[xdyn$genes<pThresh]
		  dgenes<-dgenes[!is.na(dgenes)]
		}else{
		  dgenes<-limit_to
		  dgenes<-intersect(dgenes,rownames(expX_subset))
		}
		
		grnDF<-reconstructGRN(expX_subset[dgenes,],tfs_subset,method=method,zThresh=zThresh)

		if (crossweight){
			sampTab_subset<-sampTab_subset[order(sampTab_subset[,pseudotime_column],decreasing=FALSE),]
			grnDF<-crossweight(grnDF,expX_subset[,rownames(sampTab_subset)])
		}

		res[[i]]<-grnDF
	}

	res

}

# =================== Differential Network Functions =====================
# grnDFs list of grnDF
edge_uniqueness<-function(grnDFs,tfs,weight_column){
    
    # all genes
    genes<-unlist(lapply(grnDFs, function(x) union(x$TG,x$TF)),use.names=FALSE)
    genes<-unique(genes)
    tfs<-intersect(genes,tfs)

    # cast grnDFs into adjacency matrices -- have to do this rather than just 
    # using dataframes... roundabout, but allows us to keep track of 0-weighted edges too...
    adj_list <- sapply(names(grnDFs),function(x) NULL)

    for (df in names(adj_list)){
        graph <- graph_from_data_frame(grnDFs[[df]][,c("TF","TG",weight_column)],directed=TRUE)
        addvtcs <- setdiff(genes,V(graph)$name)
        graph <- add_vertices(graph,length(addvtcs),name=addvtcs)
        adj <- as_adjacency_matrix(graph,attr=weight_column)
        
        # prune matrix
        adj <- adj[tfs,]
        adj <- t(as.matrix(adj))

        adj_list[[df]]<-adj

    }
    
    # for each edge, compute max-min edge weight (highest scoring treatment/condition - lowest).
    # This should pull out the most variable edges by magnitude.
    # Will need to check how many of these are due to nodes not in network.
    # compare against random?

    score_df<-edge_rank(adj_list)
    score_df

}

edge_rank<-function(grnMats){
    full_df<-as.data.frame(as.table(as.matrix(grnMats[[1]])))
    full_df<-full_df[,1:2]
    colnames(full_df)<-c("TG","TF")

    for (df in names(grnMats)){
        add <- as.data.frame(as.table(as.matrix(grnMats[[df]])))
        colnames(add) <- c("TG","TF",df)
        
        full_df<-merge(full_df,add,by=c("TG","TF"))
    }

    # compute min and max scores; diff score
    full_df$min<-apply(full_df[,!names(full_df) %in% c("TG","TF")],1,FUN=min)
    full_df$max<-apply(full_df[,!names(full_df) %in% c("TG","TF")],1,FUN=max)

    full_df$diff <- full_df$max - full_df$min

    # order
    full_df<-full_df[order(full_df$diff,decreasing=TRUE),]

    full_df
}


# edgeDF is the result of running 'edge_uniqueness' followed by split into dynamic network
# type DE_on vs DE_off network == "on" or "off"
dynamic_difference_network<-function(edgeDF, epochs, condition, type, diff_thresh=3, condition_thresh=6){
    epochs$mean_expression<-NULL

    conditions<-names(edgeDF[[1]])[!(names(edgeDF[[1]]) %in% c("TG","TF","min","max","diff"))]
    if (!(condition %in% conditions)){
        stop("condition not in network.")
    }
    diffnet<-lapply(edgeDF, function(x) x[,c("TG","TF",condition,"diff")])
    diffnet<-lapply(diffnet, function(x) x[x$diff>diff_thresh,])
    
    if (type=="on"){
    	diffnet<-Map(function(x,y) x[x$TF %in% y,],diffnet,epochs)
        diffnet<-lapply(diffnet, function(x) x[x[,condition]>=condition_thresh,])
    }else if (type=="off"){
        diffnet<-lapply(diffnet, function(x) x[x[,condition]<condition_thresh,])
    }else{
        stop("type should be 'on' or 'off'.") 
    }

    diffnet
}

diffnet_community_detection<-function(diffnet,method="louvain",use_weights=FALSE,weight_column=NULL){
    # if directed=TRUE, remember to flip TG and TF in diffnet dfs
    graphs<-lapply(diffnet,function(x) {g<-igraph::graph_from_data_frame(x,directed=FALSE);g})
    
    if (use_weights){
        weights<-igraph::edge_attr(x,weight_column)
    }else{
        weights<-NA
    }

    if (method=="louvain"){
    	communities<-lapply(graphs,function(x) {c<-igraph::cluster_louvain(x,weights=weights);c})
    }

    communities

}


# useful function to add correlation-based interaction type to diffnet after using dynamic_difference_network
add_type<-function(diffnet,type,grnDF_on,grnDF_offlist){
    diffnet<-diffnet[sapply(diffnet,function(x) dim(x)[1])>0]
    fun<-function(x,y){merge(x,y,by=c("TG","TF"),all.x=TRUE)}
    if(type=="on"){
        added<-lapply(diffnet,fun,grnDF_on[,c("TG","TF","corr")])
        added<-lapply(added,transform,interaction=ifelse(corr<0,"repression","activation"))
    }else if (type=="off"){
        added<-diffnet
        for (i in 1:length(grnDF_offlist)){
            added<-lapply(added,fun,grnDF_offlist[[i]][,c("TG","TF","corr")])
            #added<-lapply(added,function(x) {colnames(x)[colnames(x)=="corr"]<-paste0("othernet",i); x})
        }
        added<-lapply(added,function(x){
            cols<-colnames(x)[grepl("corr",colnames(x))];
            x$interaction<-NA;
            x[(rowSums((x[,cols]>=0) | (is.na(x[,cols])),na.rm=TRUE))==length(cols),"interaction"]<-"activation";
            x[(rowSums((x[,cols]<0) | (is.na(x[,cols])),na.rm=TRUE))==length(cols),"interaction"]<-"repression";
            x
        })
    }else{
        stop("invalid type.")
    }

    added
}


# =================== General, Broad comparison functions=====================


compute_frobenius_distance<-function(netlist1,netlist2,weight_column="zscore",compare_within_netlist1=TRUE,compare_within_netlist2=TRUE){

	# between networks in netlist1
	if (compare_within_netlist1){
		df1<-data.frame(t(combn(1:length(netlist1),2)))
		score<-c()
		for (i in 1:nrow(df1)){
			net1<-netlist1[[df1$X1[i]]]
			net2<-netlist1[[df1$X2[i]]]

			net1<-acast(net1,TF~TG,value.var=weight_column)
			net1[is.na(net1)]<-0
			net1<-t(net1)

			net2<-acast(net2,TF~TG,value.var=weight_column)
			net2[is.na(net2)]<-0
			net2<-t(net2)

			rows<-union(rownames(net1),rownames(net2))
			cols<-union(colnames(net1),colnames(net2))

			# theres a better solution than this.... 
			missing<-setdiff(rows,rownames(net1))
			addnet1<-matrix(0,nrow=length(missing),ncol=ncol(net1))
			rownames(addnet1)<-missing
			colnames(addnet1)<-colnames(net1)
			net1<-rbind(net1,addnet1)

			missing<-setdiff(rows,rownames(net2))
			addnet2<-matrix(0,nrow=length(missing),ncol=ncol(net2))
			rownames(addnet2)<-missing
			colnames(addnet2)<-colnames(net2)
			net2<-rbind(net2,addnet2)

			missing<-setdiff(cols,colnames(net1))
			addnet1<-matrix(0,nrow=nrow(net1),ncol=length(missing))
			rownames(addnet1)<-rownames(net1)
			colnames(addnet1)<-missing
			net1<-cbind(net1,addnet1)

			missing<-setdiff(cols,colnames(net2))
			addnet2<-matrix(0,nrow=nrow(net2),ncol=length(missing))
			rownames(addnet2)<-rownames(net2)
			colnames(addnet2)<-missing
			net2<-cbind(net2,addnet2)

			# order rows and columns, take difference, then norm
			net1<-net1[rows,cols]
			net2<-net2[rows,cols]

			diff<-net1-net2
			score<-c(score,norm(diff,"F"))
		}

		df1$score<-score
	}

	# between networks in netlist2
	if (compare_within_netlist2){
		df2<-data.frame(t(combn(1:length(netlist1),2)))
		score<-c()

		for (i in 1:nrow(df2)){
			net1<-netlist2[[df2$X1[i]]]
			net2<-netlist2[[df2$X2[i]]]

			net1<-acast(net1,TF~TG,value.var=weight_column)
			net1[is.na(net1)]<-0
			net1<-t(net1)

			net2<-acast(net2,TF~TG,value.var=weight_column)
			net2[is.na(net2)]<-0
			net2<-t(net2)

			rows<-union(rownames(net1),rownames(net2))
			cols<-union(colnames(net1),colnames(net2))

			# theres a better solution than this.... 
			missing<-setdiff(rows,rownames(net1))
			addnet1<-matrix(0,nrow=length(missing),ncol=ncol(net1))
			rownames(addnet1)<-missing
			colnames(addnet1)<-colnames(net1)
			net1<-rbind(net1,addnet1)

			missing<-setdiff(rows,rownames(net2))
			addnet2<-matrix(0,nrow=length(missing),ncol=ncol(net2))
			rownames(addnet2)<-missing
			colnames(addnet2)<-colnames(net2)
			net2<-rbind(net2,addnet2)

			missing<-setdiff(cols,colnames(net1))
			addnet1<-matrix(0,nrow=nrow(net1),ncol=length(missing))
			rownames(addnet1)<-rownames(net1)
			colnames(addnet1)<-missing
			net1<-cbind(net1,addnet1)

			missing<-setdiff(cols,colnames(net2))
			addnet2<-matrix(0,nrow=nrow(net2),ncol=length(missing))
			rownames(addnet2)<-rownames(net2)
			colnames(addnet2)<-missing
			net2<-cbind(net2,addnet2)

			# order rows and columns, take difference, then norm
			net1<-net1[rows,cols]
			net2<-net2[rows,cols]

			diff<-net1-net2
			score<-c(score,norm(diff,"F"))
		}

		df2$score<-score
	}

	# between networks in netlist1 and netlist2
	df3<-data.frame(expand.grid(1:length(netlist1),1:length(netlist2)))
	colnames(df3)<-c("X1","X2")
	score<-c()

	for (i in 1:nrow(df3)){
		net1<-netlist1[[df3$X1[i]]]
		net2<-netlist2[[df3$X2[i]]]

		net1<-acast(net1,TF~TG,value.var=weight_column)
		net1[is.na(net1)]<-0
		net1<-t(net1)

		net2<-acast(net2,TF~TG,value.var=weight_column)
		net2[is.na(net2)]<-0
		net2<-t(net2)

		rows<-union(rownames(net1),rownames(net2))
		cols<-union(colnames(net1),colnames(net2))

		# theres a better solution than this.... 
		missing<-setdiff(rows,rownames(net1))
		addnet1<-matrix(0,nrow=length(missing),ncol=ncol(net1))
		rownames(addnet1)<-missing
		colnames(addnet1)<-colnames(net1)
		net1<-rbind(net1,addnet1)

		missing<-setdiff(rows,rownames(net2))
		addnet2<-matrix(0,nrow=length(missing),ncol=ncol(net2))
		rownames(addnet2)<-missing
		colnames(addnet2)<-colnames(net2)
		net2<-rbind(net2,addnet2)

		missing<-setdiff(cols,colnames(net1))
		addnet1<-matrix(0,nrow=nrow(net1),ncol=length(missing))
		rownames(addnet1)<-rownames(net1)
		colnames(addnet1)<-missing
		net1<-cbind(net1,addnet1)

		missing<-setdiff(cols,colnames(net2))
		addnet2<-matrix(0,nrow=nrow(net2),ncol=length(missing))
		rownames(addnet2)<-rownames(net2)
		colnames(addnet2)<-missing
		net2<-cbind(net2,addnet2)

		# order rows and columns, take difference, then norm
		net1<-net1[rows,cols]
		net2<-net2[rows,cols]

		diff<-net1-net2
		score<-c(score,norm(diff,"F"))
	}

	df3$score<-score

	#summarize
	if (compare_within_netlist1){
		df1$group<-"netlist1"
	}
	if (compare_within_netlist2){
		df2$group<-"netlist2"
	}
	df3$group<-"cross_comparison"

	if (compare_within_netlist1 & compare_within_netlist2){
		df<-rbind(df1,df2)
		df<-rbind(df,df3)
	}else if (compare_within_netlist2){
		df<-rbind(df2,df3)
	}else if (compare_within_netlist1){
		df<-rbind(df1,df3)
	}else {
		df<-df3
	}
	
	df

}

# compute Jaccard index of top regulators
# right now only PageRank supported
compute_JI_topregs<-function(netlist1,netlist2,n_regs=15,method="pagerank",weight_column="zscore",compare_within_netlist1=TRUE,compare_within_netlist2=TRUE){
	# between networks in netlist1
	if (compare_within_netlist1){
    df1<-data.frame(t(combn(1:length(netlist1),2)))
  	score<-c()
  
  	for (i in 1:nrow(df1)){
  		net1<-netlist1[[df1$X1[i]]]
  		net2<-netlist1[[df1$X2[i]]]
  
  		net1_ranks<-compute_pagerank(list(X=net1),weight_column=weight_column)
  		net2_ranks<-compute_pagerank(list(X=net2),weight_column=weight_column)
  
  		net1_topregs<-net1_ranks$X$gene[1:n_regs]
  		net2_topregs<-net2_ranks$X$gene[1:n_regs]
  
  		ji<-length(intersect(net1_topregs,net2_topregs))/length(union(net1_topregs,net2_topregs))
  
  		score<-c(score,ji)
  	}
  
  	df1$jaccard<-score
	}

	# between networks in netlist2
	if (compare_within_netlist2){
    df2<-data.frame(t(combn(1:length(netlist2),2)))
  	score<-c()
  
  	for (i in 1:nrow(df2)){
  		net1<-netlist2[[df2$X1[i]]]
  		net2<-netlist2[[df2$X2[i]]]
  
  		net1_ranks<-compute_pagerank(list(X=net1),weight_column=weight_column)
  		net2_ranks<-compute_pagerank(list(X=net2),weight_column=weight_column)
  
  		net1_topregs<-net1_ranks$X$gene[1:n_regs]
  		net2_topregs<-net2_ranks$X$gene[1:n_regs]
  
  		ji<-length(intersect(net1_topregs,net2_topregs))/length(union(net1_topregs,net2_topregs))
  
  		score<-c(score,ji)
  	}
  
  	df2$jaccard<-score
	}

	# between networks in netlist1 and netlist2
	df3<-data.frame(expand.grid(1:length(netlist1),1:length(netlist2)))
	colnames(df3)<-c("X1","X2")
	score<-c()

	for (i in 1:nrow(df3)){
		net1<-netlist1[[df3$X1[i]]]
		net2<-netlist2[[df3$X2[i]]]

		net1_ranks<-compute_pagerank(list(X=net1),weight_column=weight_column)
		net2_ranks<-compute_pagerank(list(X=net2),weight_column=weight_column)

		net1_topregs<-net1_ranks$X$gene[1:n_regs]
		net2_topregs<-net2_ranks$X$gene[1:n_regs]

		ji<-length(intersect(net1_topregs,net2_topregs))/length(union(net1_topregs,net2_topregs))

		score<-c(score,ji)
	}

	df3$jaccard<-score

	#summarize
	if (compare_within_netlist1){
	  df1$group<-"netlist1"
	}
	if (compare_within_netlist2){
	  df2$group<-"netlist2"
	}
	df3$group<-"cross_comparison"
	
	if (compare_within_netlist1 & compare_within_netlist2){
	  df<-rbind(df1,df2)
	  df<-rbind(df,df3)
	}else if (compare_within_netlist2){
	  df<-rbind(df2,df3)
	}else if (compare_within_netlist1){
	  df<-rbind(df1,df3)
	}else {
	  df<-df3
	}
	
	df

}

JI_across_topregs<-function(netlist1,netlist2,n_regs=3:15,func="mean",method="pagerank",weight_column="zscore",compare_within_netlist1=TRUE,compare_within_netlist2=TRUE){
  if (func=="mean"){
    res<-data.frame(group=character(),jaccard=numeric(),n_regs=numeric())
    for (i in n_regs){
      ji<-compute_JI_topregs(netlist1,netlist2,n_regs=i,method=method,weight_column=weight_column,compare_within_netlist1 = compare_within_netlist1,compare_within_netlist2 = compare_within_netlist2)
      ji<-ji[,c("group","jaccard")]
      ji<-aggregate(.~group,ji,mean)
      ji$n_regs<-i
      
      res<-rbind(res,ji)
    }
    
  }
  res
}


tally_topregs<-function(netlist1,netlist2,n_regs=15,method="pagerank",weight_column="zscore"){

	# in netlist1
	topregs_1<-c()
	for (i in 1:length(netlist1)){
		net1<-netlist1[[i]]
		net1_ranks<-compute_pagerank(list(X=net1),weight_column=weight_column)
		net1_topregs<-net1_ranks$X$gene[1:n_regs]

		topregs_1<-c(topregs_1,net1_topregs)
	}

	# in netlist2
	topregs_2<-c()
	for (i in 1:length(netlist2)){
		net2<-netlist2[[i]]
		net2_ranks<-compute_pagerank(list(X=net2),weight_column=weight_column)
		net2_topregs<-net2_ranks$X$gene[1:n_regs]

		topregs_2<-c(topregs_2,net2_topregs)
	}

	# summarize
	topregs_1<-table(topregs_1)
	topregs_2<-table(topregs_2)

	list(netlist1=topregs_1, netlist2=topregs_2)

}

average_pagerank_ranks<-function(netlist1,weight_column="zscore",max_rank=NULL){

	if (is.null(names(netlist1))){
		names(netlist1)<-paste0("X",seq(1:length(netlist1)))
	}

	if (length(netlist1)==1){
		net_ranks<-compute_pagerank(netlist1, weight_column=weight_column)
		net_ranks<-lapply(net_ranks, function(x) x[x$is_regulator,])
		net_ranks<-lapply(net_ranks, function(x) cbind(x,rank=1:nrow(x)))
		net_ranks<-lapply(net_ranks, function(x) x[,c("gene","rank")])
		res<-net_ranks[[1]]
		colnames(res)<-c("gene","mean_rank")
		return(res)
	}

	net_ranks<-compute_pagerank(netlist1,weight_column=weight_column)
	net_ranks<-lapply(net_ranks, function(x) x[x$is_regulator,])
	net_ranks<-lapply(net_ranks, function(x) cbind(x,rank=1:nrow(x)))
	net_ranks<-lapply(net_ranks, function(x) x[,c("gene","rank")])
	net_ranks<-lapply(net_ranks, function(x) as.data.frame(t(x)))

	if (is.null(max_rank)){
		max_rank<-max(unlist(lapply(net_ranks,function(x) nrow(x))))
	}

	ranks<-plyr::rbind.fill(net_ranks,fill=NA)

	# there's a better way to do this... stupid hacky dataframe manipulating... could've just used merge() and set all=TRUE
	ranks<-ranks[seq(2,nrow(ranks),2),]

	ranks<-data.frame(sapply(ranks,function(x) as.numeric(as.character(x))))
	ranks[is.na(ranks)]<-max_rank
	avgranks<-colMeans(ranks)

	res<-data.frame(gene=names(avgranks),mean_rank=avgranks)
	res<-res[order(res$mean_rank,decreasing=FALSE),]

	res

}


# function that takes a network as grnDF, returns weighted betweennes of TFs
compute_betweenness<-function(grnDF,tfs=NULL,weight_column="zscore",normalized=TRUE){
	if (is.null(tfs)){
		tfs<-unique(grnDF$TF)
	}

	grnDF<-grnDF[,c("TF","TG",weight_column)]
	colnames(grnDF)<-c("TF","TG","weight")
	g<-graph_from_data_frame(grnDF,directed=FALSE)

	betweenness<-betweenness(g,tfs,directed=FALSE,normalized=normalized)
	betweenness

}

# function that takes a list of network dataframes
# returns long dataframe listing network, betweenness of each TF, degree of each TF
biglist_compute_betweenness_degree<-function(netlist,weight_column="zscore",normalized=TRUE){

	res<-data.frame(betweenness=numeric(),degree=numeric(),network=character(),TF=character())

	if (is.null(names(netlist))){
		names(netlist)<-seq(1:length(netlist))
	}

	for (name in names(netlist)){
		net<-netlist[[name]]
		g<-graph_from_data_frame(net[,c("TF","TG")],directed=FALSE)
		betweenness<-compute_betweenness(net,weight_column=weight_column,normalized=normalized)
		degree<-degree(g,unique(net$TF),mode="all",normalized=normalized)

		df<-cbind(betweenness,degree=degree[names(betweenness)])
		df<-as.data.frame(df)
		df$network<-name
		df$TF<-rownames(df)
		
		res<-rbind(res,df)

	}

	res

}


# =================== Useful Plotting Functions =====================
# plotting the diffnet
plot_dyn_diffnet<-function(grn,tfs,only_TFs=TRUE,order=NULL){
  g<-list()

  if (!is.null(order)){
    grn<-grn[order]
  } 

  for (i in 1:length(grn)){
    df<-grn[[i]]
    
    if (only_TFs){
      df<-df[df$TG %in% tfs,]
    }
    if (nrow(df)==0){
    	next
    }

    net<-graph_from_data_frame(df[,c("TF","TG","interaction")],directed=FALSE)
    layout<-layout_with_fr(net)
    rownames(layout)<-V(net)$name
    layout_ordered<-layout[V(net)$name,]
    tfnet<-ggnetwork(net,layout=layout_ordered,cell.jitter=0)
    tfnet$is_regulator<-as.character(tfnet$name %in% tfs)

    cols<-c("activation"="blue","repression"="red")
    g[[i]]<-ggplot()+
      geom_edges(data=tfnet,aes(x=x, y=y, xend=xend, yend=yend, color=interaction),size=0.75,curvature=0.1, alpha=.6)+
      scale_color_manual(values=cols)+
      geom_nodes(data=tfnet,aes(x=x, y=y, xend=xend, yend=yend),color="darkgray",size=6,alpha=.5)+
      geom_nodes(data=tfnet[tfnet$is_regulator=="TRUE",],aes(x=x, y=y, xend=xend, yend=yend),color="#8C4985",size=6,alpha=.8)+
      geom_nodelabel_repel(data=tfnet,aes(x=x, y=y, label=name),size=2.5, color="#5A8BAD")+
      theme_blank()+
      ggtitle(names(grn)[i])

    g[[i]]<-g[[i]]+theme(legend.position="none")
  }
  g<-g[!sapply(g,is.null)]
  do.call(grid.arrange,g)
}


# same as above, but will also do community detection and betweenness computation, and color/fade accordingly
plot_diffnet_detail<-function(grn,tfs,only_TFs=TRUE,order=NULL,weight_column="zscore",compute_betweenness=TRUE){
  g<-list()

  if (!is.null(order)){
    grn<-grn[order]
  } 

  for (i in 1:length(grn)){
    df<-grn[[i]]
    
    if (only_TFs){
      df<-df[df$TG %in% tfs,]
    }
    if (nrow(df)==0){
    	next
    }

    df<-df[,c("TF","TG",weight_column,"interaction")]
    colnames(df)<-c("TF","TG","weight","interaction")

    net<-graph_from_data_frame(df,directed=FALSE)

    if(compute_betweenness){
	    b<-betweenness(net,directed=FALSE,normalized=TRUE)
	    b<-as.data.frame(b)
	    colnames(b)<-"betweenness"
	    b$gene<-rownames(b)
	    b<-b[,c("gene","betweenness")]
	    c<-as.data.frame(as.table(membership(cluster_louvain(net))))
	    colnames(c)<-c("gene","communities")
	    vtx_features<-merge(b,c,by="gene",all=TRUE)
	}else{
		c<-as.data.frame(as.table(membership(cluster_louvain(net))))
	    colnames(c)<-c("gene","communities")
	    vtx_features<-c
	}

    layout<-layout_with_fr(net)
    rownames(layout)<-V(net)$name
    layout_ordered<-layout[V(net)$name,]
    tfnet<-ggnetwork(net,layout=layout_ordered,cell.jitter=0)
    tfnet$is_regulator<-as.character(tfnet$name %in% tfs)
    if(compute_betweenness){
    	tfnet$betweenness<-vtx_features$betweenness[match(tfnet$name,vtx_features$gene)]
    }
    tfnet$communities<-as.factor(vtx_features$communities[match(tfnet$name,vtx_features$gene)])

    cols<-c("activation"="blue","repression"="red")
    num_cols2<-length(unique(tfnet$communities))
    if (num_cols2<=8){
      cols2<-brewer.pal(num_cols2,"Set1")
    }else{
      cols2<-colorRampPalette(brewer.pal(8,"Set1"))(num_cols2)
    }
    names(cols2)<-unique(tfnet$communities)
    cols<-c(cols,cols2)

    if(compute_betweenness){
	    g[[i]]<-ggplot()+
	      geom_edges(data=tfnet,aes(x=x, y=y, xend=xend, yend=yend, color=interaction),size=0.75,curvature=0.1, alpha=.6)+
	      scale_color_manual(values=cols)+
	      geom_nodes(data=tfnet,aes(x=x, y=y, xend=xend, yend=yend),color="darkgray",size=6,alpha=.5)+
	      geom_nodes(data=tfnet[tfnet$is_regulator=="TRUE",],aes(x=x, y=y, xend=xend, yend=yend,color=communities, alpha=betweenness+.5),size=6)+
	      geom_nodelabel_repel(data=tfnet,aes(x=x, y=y, label=name),size=2.5, color="#5A8BAD")+
	      theme_blank()+
	      ggtitle(names(grn)[i])
  	}else{
  		g[[i]]<-ggplot()+
	      geom_edges(data=tfnet,aes(x=x, y=y, xend=xend, yend=yend, color=interaction),size=0.75,curvature=0.1, alpha=.6)+
	      scale_color_manual(values=cols)+
	      geom_nodes(data=tfnet,aes(x=x, y=y, xend=xend, yend=yend),color="darkgray",size=6,alpha=.5)+
	      geom_nodes(data=tfnet[tfnet$is_regulator=="TRUE",],aes(x=x, y=y, xend=xend, yend=yend,color=communities),alpha=.5,size=6)+
	      geom_nodelabel_repel(data=tfnet,aes(x=x, y=y, label=name),size=2.5, color="#5A8BAD")+
	      theme_blank()+
	      ggtitle(names(grn)[i])
  	}

    g[[i]]<-g[[i]]+theme(legend.position="none")
  }
  g<-g[!sapply(g,is.null)]
  do.call(grid.arrange,g)
}

