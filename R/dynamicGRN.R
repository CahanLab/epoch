
# first attempt at assigning genes to epochs as time periods (rather than clustering genes in epochs)
# method is one of "cell_order", "pseudotime" or "group"

#' Assigns genes to epochs
#'
#' Assigns genes to epochs
#' 
#' @param expSmoothed smoothed expression matrix
#' @param expDat expression matrix
#' @param sampTab sample table
#' @param dynRes result of running findDynGenes
#' @param num_epochs number of epochs, ignored if method = "group"
#' @param method How to divide time into epochs -- one of "cell_order", "pseudotime", or "group"
#' @param pThresh pval threshold
#' @param toScale whether or not to scale the data
#'
#' @return epochs a list detailing genes active in each epoch
#' @export
#'
assign_epochs<-function(expSmoothed,expDat,dynRes,num_epochs=3,method="pseudotime",pThresh=0.01,toScale=FALSE){

	# limit exp to dynamically expressed genes
	exp<-expSmoothed[names(dynRes$genes[dynRes$genes<pThresh]),]
	# scale the data if needed
	if (toScale){
		if(class(exp)[1]!='matrix'){
        	exp <- t(scale(Matrix::t(exp)))
      	}else{
        	exp <- t(scale(t(exp)))
      	}
	}

	navg<-ceiling(ncol(exp)*0.05)

	# compute thresholds for each gene
	# for each gene, order expression --- compute top as average of top 5%, bottom as average of bottom 5% 
	# set threshold (active/inactive) as midpoint (or maybe 1/3 in case gene is expressed at different levels) between top and bottom
	thresholds<-data.frame(gene=rownames(exp),thresh=rep(0,length(rownames(exp))))
	rownames(thresholds)<-thresholds$gene
	for (gene in rownames(exp)){
		profile<-exp[gene,][order(exp[gene,],decreasing=FALSE)]
		bottom<-mean(profile[1:navg])
		top<-mean(profile[(length(profile)-navg):length(profile)])

		thresh<-((top-bottom)*0.33) + bottom
		thresholds[gene,"thresh"]<-thresh
	}

	# order cells in exp along pseudotime-- cells ordered in dynRes
	t1 = dynRes$cells$pseudotime
   	names(t1) = as.vector(dynRes$cells$cell_name)
   	sort(t1, decreasing=FALSE)
	exp<-exp[,names(t1)]

	mean_expression<-data.frame(gene=character(),epoch=numeric(),mean_expression=numeric())
	if (method=="cell_order"){
		# Divide epochs evenly by cell order number; assumes even sampling of cells across pseudotime

		# create empty list to contain active genes
		epoch_names<-paste0(rep("epoch",num_epochs),seq(1:num_epochs))
		epochs<-vector("list",length(epoch_names))
		names(epochs)<-epoch_names

		# determine activity based on average expression in each epoch
		chunk_size<-floor(ncol(exp)/num_epochs)
		for (i in 1:length(epochs)){
			chunk<-exp[,((i-1)*chunk_size):(i*chunk_size)]
			chunk_df<-data.frame(means=rowMeans(chunk))
			chunk_df<-cbind(chunk_df,thresholds)
			chunk_df$active<-(chunk_df$means>=chunk_df$thresh)

			epochs[[i]]<-rownames(chunk_df[chunk_df$active,])

			mean_expression<-rbind(mean_expression,data.frame(gene=rownames(chunk),epoch=rep(i,length(rownames(chunk))),mean_expression=rowMeans(chunk)))

		}

	}else if (method=="pseudotime"){
		# Divide epochs by pseudotime

		# create empty list to contain active genes
		epoch_names<-paste0(rep("epoch",num_epochs),seq(1:num_epochs))
		epochs<-vector("list",length(epoch_names))
		names(epochs)<-epoch_names
		
		# determine activity based on average expression in each epoch
		ptmax<-max(dynRes$cells$pseudotime)
		ptmin<-min(dynRes$cells$pseudotime)
		chunk_size<-(ptmax-ptmin)/num_epochs
		for (i in 1:length(epochs)){
			lower_bound<-ptmin+((i-1)*chunk_size)
			upper_bound<-ptmin+(i*chunk_size)
			chunk_cells<-rownames(dynRes$cells[dynRes$cells$pseudotime>=lower_bound & dynRes$cells$pseudotime<=upper_bound,])
			chunk<-exp[,chunk_cells]

			chunk_df<-data.frame(means=rowMeans(chunk))
			chunk_df<-cbind(chunk_df,thresholds)
			chunk_df$active<-(chunk_df$means>=chunk_df$thresh)

			epochs[[i]]<-rownames(chunk_df[chunk_df$active,])

			mean_expression<-rbind(mean_expression,data.frame(gene=rownames(chunk),epoch=rep(i,length(rownames(chunk))),mean_expression=rowMeans(chunk)))

		}


	}else if (method=="group"){
		# Assign epochs based on group in dynRes$cells$group; ingores num_epochs

		# create empty list to contain active genes
		epoch_names<-unique(dynRes$cells$group)
		epochs<-vector("list",length(epoch_names))
		names(epochs)<-epoch_names

		# determine activity based on average expression in each epoch
		for (epoch in names(epochs)){
			chunk_cells<-rownames(dynRes$cells[dynRes$cells$group==epoch,])
			chunk<-exp[,chunk_cells]

			chunk_df<-data.frame(means=rowMeans(chunk))
			chunk_df<-cbind(chunk_df,thresholds)
			chunk_df$active<-(chunk_df$means>=chunk_df$thresh)

			epochs[[epoch]]<-rownames(chunk_df[chunk_df$active,])

			mean_expression<-rbind(mean_expression,data.frame(gene=rownames(chunk),epoch=rep(i,length(rownames(chunk))),mean_expression=rowMeans(chunk)))

		}

	}else{
		stop('invalid method.')
	}

	mean_expression$epoch<-paste0("epoch",mean_expression$epoch)
	epochs$mean_expression<-mean_expression
	epochs

}

#' Divides grnDF into epochs, filters interactions between genes not in same or consecutive epochs 
#'
#' @param grnDF result of GRN reconstruction
#' @param epochs result of running assign_epochs
#' @param epoch_network dataframe outlining higher level epoch connectivity (i.e. epoch transition network). If NULL, will assume epochs is ordered linear trajectory
#'
#' @return list of GRNs across epochs and transitions
#' @export
#'
epochGRN<-function(grnDF, epochs, epoch_network=NULL){
	
	epochs$mean_expression<-NULL
	all_dyngenes<-unique(unlist(epochs,use.names = FALSE))

	# assign epoch_network if NULL
	if (is.null(epoch_network)){
		
		epoch_network<-data.frame(from=character(), to=character())
		for (i in 1:(length(names(epochs))-1)){
			df<-data.frame(from=c(names(epochs)[i]),to=c(names(epochs)[i+1]))
			epoch_network<-rbind(epoch_network,df)
		}

	}

	# add in self loops to epoch_network (i.e. for interactions ocurring inside same epoch)
	epoch_network<-rbind(epoch_network,data.frame(from=names(epochs),to=names(epochs)))

	# build network for each epoch and transition
	# this will inherently filter out interactions not belonging to same or consecutive epochs (i.e. anything not in epoch_network)
	epoch_network$name<-paste(epoch_network[,1],epoch_network[,2],sep="..")
	GRN<-vector("list",nrow(epoch_network))
	names(GRN)<-epoch_network$name

	print(epoch_network)

	for (t in 1:nrow(epoch_network)){
		from<-as.character(epoch_network[t,1])
		to<-as.character(epoch_network[t,2])

		# TF must be active in source epoch
		temp<-grnDF[grnDF$TF %in% epochs[[from]],]

		# For transition network
		if (from!=to){
			# target turns on: target is not active in source epoch but is active in target epoch

			# target turns off: target is active in source epoch but is not active in target epoch

			# remove any other interaction (i.e. interactions that are constant -- target on in both epochs or off in both epochs)
			remove_tgs_in_both<-intersect(epochs[[to]],epochs[[from]])
			remove_tgs_in_neither<-intersect(setdiff(all_dyngenes,epochs[[from]]),setdiff(all_dyngenes,epochs[[to]]))

			temp<-temp[!(temp$TG %in% remove_tgs_in_both),]
			temp<-temp[!(temp$TG %in% remove_tgs_in_neither),]

		}
		# Else, For epoch network (non-transition network), keep all interactions as long as TF is active

		GRN[[epoch_network[t,"name"]]]<-temp

	}
	
	GRN

}


#' Function to compute page rank of TF+target networks
#'
#' @param dynnet result of dynamic GRN reconstruction
#' @param weight_column name of column in dynnet to weight edges
#' @param directed_graph if GRN is directed or not
#'
#' @return list of dataframes of active genes in each epoch and transition, ranked by page rank
#' @export
#'
compute_pagerank<-function(dynnet,weight_column="zscore",directed_graph=FALSE){

  ranks<-sapply(names(dynnet), function (x) NULL)
  for (net in names(dynnet)){
    df<-dynnet[[net]]
    
    df<-df[,c("TF","TG",weight_column)]
    colnames(df)<-c("TF","TG","weight")
    
    tfnet<-graph_from_data_frame(df,directed=directed_graph)
    
    pagerank<-data.frame(page_rank(tfnet,directed=directed_graph)$vector)
    colnames(pagerank)<-c("page_rank")
    pagerank$gene<-rownames(pagerank)
    pagerank<-pagerank[,c("gene","page_rank")]
    pagerank<-pagerank[order(pagerank$page_rank,decreasing=TRUE),]

    pagerank$is_regulator<-FALSE
    pagerank$is_regulator[pagerank$gene %in% unique(df$TF)]<-TRUE
    
    ranks[[net]]<-pagerank
  }
  
  ranks
}


#' Function to compute page rank of TF-only networks
#'
#' @param dynnet result of dynamic GRN reconstruction
#' @param tfs TFs
#' @param weight_column name of column in dynnet to weight edges
#' @param directed_graph if GRN is directed or not
#'
#' @return list of dataframes of active TFs in each epoch and transition, ranked by page rank
#' @export
tfnet_pagerank<-function(dynnet,tfs,weight_column="zscore",directed_graph=FALSE){

  tfranks<-sapply(names(dynnet), function (x) NULL)
  for (net in names(dynnet)){
    df<-dynnet[[net]]
    #limit dynamic network to TFs
    df<-df[df$TG %in% tfs,]
    
    df<-df[,c("TF","TG",weight_column)]
    colnames(df)<-c("TF","TG","weight")
    
    tfnet<-graph_from_data_frame(df,directed=directed_graph)
    
    pagerank<-data.frame(page_rank(tfnet,directed=directed_graph)$vector)
    colnames(pagerank)<-c("page_rank")
    pagerank$TF<-rownames(pagerank)
    pagerank<-pagerank[,c("TF","page_rank")]
    pagerank<-pagerank[order(pagerank$page_rank,decreasing=TRUE),]
    
    tfranks[[net]]<-pagerank
  }
  
  tfranks
}













