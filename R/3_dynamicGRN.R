
# =================== Functions to define, assign, and extract the dynamic network =====================


#' Define epochs
#'
#' Define epochs
#' 
#' @param dynRes results of running findDynGenes, or a list of results of running findDynGenes per path. If list, names should match names of expDat.
#' @param expDat genes-by-cells expression matrix, or a list of expression matrices per path. If list, names should match names of dynRes.
#' @param method method to define epochs. Either "pseudotime", "cell_order", "group", "con_similarity", "kmeans", "hierarchical"
#' @param num_epochs number of epochs to define. Ignored if epoch_transitions, pseudotime_cuts, or group_assignments are provided.
#' @param pseudotime_cuts vector of pseudotime cutoffs. If NULL, cuts are set to max(pseudotime)/num_epochs.
#' @param group_assignments a list of vectors where names(assignment) are epoch names, and vectors contain groups belonging to corresponding epoch
#'
#' @return updated list of dynRes with epoch column included in dynRes$cells
#' @export
#'
define_epochs<-function(dynRes,
						expDat,
						method="pseudotime",
						num_epochs=2,
						pseudotime_cuts=NULL,
						group_assignments=NULL,
						pThresh_dyn=0.05,
						winSize=2){

	# put dynRes and expDat into list, if not already
	# for cleaner code later on... just put everything is in list.
	if(class(dynRes[[1]])!="list"){
		dynRes<-list(dynRes)
		expDat<-list(expDat)
	}

	new_dynRes<-list()
	for (path in 1:length(dynRes)){
		if (!is.null(names(dynRes))){
			path<-names(dynRes)[path]
		}
		path_dyn<-dynRes[[path]]
		path_dyn$cells$epoch<-NA

		# Define epochs by cell order: equal number of cells in each epoch
		if (method=="cell_order"){
			t1 = path_dyn$cells$pseudotime
   			names(t1) = as.vector(path_dyn$cells$cell_name)
   			sort(t1, decreasing=FALSE)
   			chunk_size<-floor(length(t1)/num_epochs)

   			for (i in 1:num_epochs){
   				if (i==num_epochs){
   					cells_in_epoch<-names(t1)[(1+((i-1)*chunk_size)):length(t1)]
   				}else{
   					cells_in_epoch<-names(t1)[(1+((i-1)*chunk_size)):(i*chunk_size)]
   				}
   				path_dyn$cells$epoch[path_dyn$cells$cell_name %in% cells_in_epoch]<-paste0("epoch",i)
   			}
   			
		}

		# Define epochs by pseudotime: use pseudotime cutoffs
		if (method=="pseudotime"){
			if (is.null(pseudotime_cuts)){
				# if pseudotime_cuts NULL, then split pseudotime evenly (note: this is not the same as using cell_order)
				pseudotime_cuts<-seq(min(path_dyn$cells$pseudotime),max(path_dyn$cells$pseudotime),(max(path_dyn$cells$pseudotime)-min(path_dyn$cells$pseudotime))/num_epochs)
				pseudotime_cuts<-pseudotime_cuts[-length(pseudotime_cuts)]
				pseudotime_cuts<-pseudotime_cuts[-1]
			}

			path_dyn<-split_epochs_by_pseudotime(path_dyn,pseudotime_cuts)
		}


		# Define epochs by group: use clusters to separate cells into epochs
		if (method=="group"){
			if (is.null(group_assignments)){
				stop("Must provide group_assignments for group method.")
			}

			path_dyn<-split_epochs_by_group(path_dyn,group_assignments)
		}

		# define by consecutive similarity
		if (method=="con_similarity"){
			cuts = find_cuts_by_similarity(expDat[[path]], path_dyn, winSize=winSize, pThresh_dyn=pThresh_dyn)
			path_dyn<-split_epochs_by_pseudotime(path_dyn,cuts)
		}

		# define by kmeans
		if (method=="kmeans"){
			cuts = find_cuts_by_clustering(expDat[[path]], path_dyn, num_epochs=num_epochs, method="kmeans", pThresh_dyn=pThresh_dyn)
			path_dyn<-split_epochs_by_pseudotime(path_dyn,cuts)
		}

		# define by hierarchical clustering
		if (method=="hierarchical"){
			cuts = find_cuts_by_clustering(expDat[[path]], path_dyn, num_epochs=num_epochs, method="hierarchical", pThresh_dyn=pThresh_dyn)
			path_dyn<-split_epochs_by_pseudotime(path_dyn,cuts)
		}


		new_dynRes[[path]]<-path_dyn

	}

	if (length(new_dynRes)==1){
		new_dynRes<-new_dynRes[[1]]
	}
	new_dynRes

}

#' Splits data into epochs
#'
#' Splits data into epochs by assigning cells to epochs
#' 
#' @param dynRes result of running findDynGenes or compileDynGenes
#' @param cuts vector of pseudotime cutoffs
#' @param epoch_names names of resulting epochs, must have length of length(cuts)+1
#'
#' @return updated dynRes with epoch column included in dynRes$cells
#' @export
#'
split_epochs_by_pseudotime<-function(dynRes,cuts,epoch_names=NULL){

	sampTab<-dynRes$cells

	if (max(cuts)>max(sampTab$pseudotime)){
		stop("Cuts must be within pseudotime.")
	}

	if (!is.null(epoch_names) & (length(epoch_names)!=length(cuts)+1)){
		stop("Length of epoch_names must be equal to 1+length(cuts).")
	}

	if (is.null(epoch_names)){
		epoch_names<-paste0("epoch",(1:(length(cuts)+1)))
	}

	cuts<-c(-0.1,cuts,max(sampTab$pseudotime))
	sampTab$epoch<-NA
	for (i in 2:length(cuts)){
		sampTab$epoch[(cuts[i-1]<sampTab$pseudotime) & (sampTab$pseudotime<=cuts[i])]<-epoch_names[i-1]
	}

	dynRes$cells<-sampTab
	dynRes

}

#' Splits data into epochs manually
#'
#' Splits data into epochs given group assignment
#' 
#' @param dynRes result of running findDynGenes or compileDynGenes
#' @param assignment a list of vectors where names(assignment) are epoch names, and vectors contain groups belonging to corresponding epoch
#'
#' @return updated dynRes with epoch column included in dynRes$cells
#' @export
#'
split_epochs_by_group<-function(dynRes,assignment){

	sampTab<-dynRes$cells
	sampTab$epoch<-NA

	for (e in names(assignment)){
		sampTab$epoch[sampTab$group %in% assignment[[e]]]<-e
	}

	dynRes$cells<-sampTab
	dynRes

}


# #' Returns cuts to define epochs
# #'
# #' Returns cuts to define epochs via sliding window comparison
# #' 
# #' @param expDat genes-by-cells expression matrix
# #' @param dynRes result of running findDynGenes or define_epochs
# #' @param winSize number of cells to each side to compare each cell to
# #' @param pThresh_dyn pval threshold if gene is dynamically expressed
# #'
# #' @return vector of  pseudotimes at which to cut data into epochs
# #' @export
# #'
# find_cuts_by_similarity<-function(
# 	expDat,
# 	dynRes,
# 	winSize=5,
# 	limit_to=NULL,
# 	pThresh_dyn=0.05){
	
#   # limit exp to dynamically expressed genes
# 	expDat<-expDat[names(dynRes$genes[dynRes$genes<pThresh_dyn]),]
# 	#cat(nrow(expDat),"\n")
	
# 	if (!is.null(limit_to)){
#     	limit_to<-intersect(limit_to,rownames(expDat))
#     	expDat<-expDat[limit_to,]
#     }

# 	# compute PCC between all cells -- just easier this way
# 	xcorr = cor(expDat[,rownames(dynRes$cells)])
# 	xdist = 1 - xcorr
# 	mydists = rep(1, nrow(xdist)-1)
# 	pShift = rep(0, nrow(xdist)-1)
# 	end = nrow(xdist)
# 	for(i in 2 : end-1){
#   		xi = seq(i-1 : max(1, i-winSize))
#   		yi = seq(i+1 : min(end, i+winSize))
#   		choice = which.max( c(mean(xdist[i,xi]), mean(xdist[i,yi])))
#   		pShift[i] = c(0,1)[choice]
#   		mydists[i] = mean(xdist[i,xi]) - mean(xdist[i,yi])
# 	}
# 	# the first one is bogus
# 	cuts_index = which(pShift==1)
# 	cuts_index = cuts_index[2:length(cuts_index)]
# 	cat("Cut points: ",dynRes$cells[cuts_index,]$pseudotime, "\n")

# 	# It is kinda confusing, but we want the pt of the cell just prior to the determined cutpoint
# 	# But, because of the way that this is indexed, we don't need to adjust anything
# 	dynRes$cells[cuts_index,]$pseudotime
# }


# same thing as original find_cuts_by_similarity, but code is easier to understand

#' Returns cuts to define epochs
#'
#' Returns cuts to define epochs via sliding window comparison
#' 
#' @param expDat genes-by-cells expression matrix
#' @param dynRes result of running findDynGenes or define_epochs
#' @param winSize number of cells to each side to compare each cell to
#' @param limit_to vector of genes on which to base epoch cuts, for example, limiting to TFs
#' @param pThresh_dyn pval threshold if gene is dynamically expressed
#'
#' @return vector of  pseudotimes at which to cut data into epochs
#' @export
#'
find_cuts_by_similarity<-function(
    expDat,
    dynRes,
    winSize=2,
    limit_to=NULL,
    pThresh_dyn=0.05){
    
  	# limit exp to dynamically expressed genes
    expDat<-expDat[names(dynRes$genes[dynRes$genes<pThresh_dyn]),]
    #cat(nrow(expDat),"\n")
    
    if (!is.null(limit_to)){
    	limit_to<-intersect(limit_to,rownames(expDat))
    	expDat<-expDat[limit_to,]
    }

    # compute PCC between all cells -- just easier this way
    xcorr = cor(expDat[,rownames(dynRes$cells)])
    pShift = rep(0, nrow(xcorr)-1)
    end = nrow(xcorr)

    for (i in 2:end-1){
    	past = seq(i-1 : max(1, i-winSize))
    	pastandfuture = seq(i+1 : min(end, i+winSize))
    	closer_to = which.max( c(mean(xcorr[i,past]), mean(xcorr[i,pastandfuture])))
    	pShift[i] = c(0,1)[closer_to]
    }

    consecutive_diffs<-diff(pShift)
	cuts_index = which(consecutive_diffs==1)

    #cat("Cut indicies: ",cuts_index, "\n")
    cat("Cut points: ",dynRes$cells[cuts_index,]$pseudotime, "\n")

    # It is kinda confusing, but we want the pt of the cell just prior to the determined cutpoint
    # But, because of the way that this is indexed, we don't need to adjust anything
    dynRes$cells[cuts_index,]$pseudotime
}


#' Returns cuts to define epochs
#'
#' Returns cuts to define epochs via clustering
#' 
#' @param expDat genes-by-cells expression matrix
#' @param dynRes result of running findDynGenes or define_epochs
#' @param num_epochs the number of epochs
#' @param limit_to vector of genes on which to base epoch cuts, for example, limiting to TFs
#' @param method what clustering method to use, either 'kmeans' or 'hierarchical'
#' @param pThresh_dyn pval threshold if gene is dynamically expressed
#'
#' @return vector of  pseudotimes at which to cut data into epochs
#' @export
#'
find_cuts_by_clustering<-function(expDat, 
									dynRes, 
									num_epochs,
									limit_to=NULL,
									method="kmeans",
									pThresh_dyn=0.05){

	expDat<-expDat[names(dynRes$genes[dynRes$genes<pThresh_dyn]),]
    #cat(nrow(expDat),"\n")

    if (!is.null(limit_to)){
    	limit_to<-intersect(limit_to,rownames(expDat))
    	expDat<-expDat[limit_to,]
    }

    if (method=="kmeans"){
    	clustering<-kmeans(t(expDat),num_epochs,iter.max=100)$cluster
    	
    	cuts<-c()
    	for (cluster in unique(clustering)){
    		cells<-names(clustering)[clustering==cluster]
    		cuts<-c(cuts,max(dynRes$cells$pseudotime[rownames(dynRes$cells) %in% cells]))
    	}
    	cuts<-sort(cuts)
    	cuts<-cuts[1:length(cuts)-1]
    }else if (method=="hierarchical"){
    	clustering<-hclust(dist(t(expDat)),method="centroid")
    	clusterCut <- cutree(clustering, 3)
    	cuts<-c()
    	for (cluster in unique(clusterCut)){
    		cells<-names(clusterCut)[clusterCut==cluster]
    		cuts<-c(cuts,max(dynRes$cells$pseudotime[rownames(dynRes$cells) %in% cells]))
    	}
    	cuts<-sort(cuts)
    	cuts<-cuts[1:length(cuts)-1]
    }else{
    	stop("method must be either 'kmeans' or 'hierarchical'.")
    }

    cat("Cut points: ",cuts, "\n")
    cuts
}




#' Assigns genes to epochs
#'
#' Assigns genes to epochs
#' 
#' @param expDat genes-by-cells expression matrix 
#' @param dynRes individual path result of running define_epochs
#' @param method method of assigning epoch genes, either "active_expression" (looks for active expression in epoch) or "DE" (looks for differentially expressed genes per epoch)
#' @param pThresh_dyn pval threshold if gene is dynamically expressed
#' @param pThresh_DE pval if gene is differentially expressed. Ignored if method is active_expression.
#' @param active_thresh value between 0 and 1. Percent threshold to define activity
#' @param toScale whether or not to scale the data
#' @param forceGenes whether or not to rescue orphan dyanmic genes, forcing assignment into epoch with max expression.
#'
#' @return epochs a list detailing genes active in each epoch
#' @export
#'
assign_epochs<-function(expDat,
						dynRes,
						method='active_expression',
						pThresh_dyn=0.05,
						pThresh_DE=0.05,
						active_thresh=0.33,
						toScale=FALSE,
						forceGenes=TRUE){

	if (active_thresh<0 | active_thresh>1){
		stop("active_thresh must be between 0 and 1.")
	}

	# limit exp to dynamically expressed genes
	exp<-expDat[names(dynRes$genes[dynRes$genes<pThresh_dyn]),]
	# scale the data if needed
	if (toScale){
		if(class(exp)[1]!='matrix'){
        	exp <- t(scale(Matrix::t(exp)))
      	}else{
        	exp <- t(scale(t(exp)))
      	}
	}

	# Assign epochs based on assignment in dynRes$cells$epoch
	# create empty list to contain active genes
	epoch_names<-unique(dynRes$cells$epoch)
	epochs<-vector("list",length(epoch_names))
	names(epochs)<-epoch_names

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

		thresh<-((top-bottom)*active_thresh) + bottom
		thresholds[gene,"thresh"]<-thresh
	}

	mean_expression<-data.frame(gene=character(),epoch=numeric(),mean_expression=numeric())
	if (method=="active_expression"){
		# determine activity based on average expression in each epoch
		for (epoch in names(epochs)){
			chunk_cells<-rownames(dynRes$cells[dynRes$cells$epoch==epoch,])
			chunk<-exp[,chunk_cells]

			chunk_df<-data.frame(means=rowMeans(chunk))
			chunk_df<-cbind(chunk_df,thresholds)
			chunk_df$active<-(chunk_df$means>=chunk_df$thresh)

			epochs[[epoch]]<-rownames(chunk_df[chunk_df$active,])

			mean_expression<-rbind(mean_expression,data.frame(gene=rownames(chunk),epoch=rep(epoch,length(rownames(chunk))),mean_expression=rowMeans(chunk)))
		}

	}else{
		# determine which genes are differentially expressed per epoch

		# DESeq will require raw count data
		# for (epoch in names(epochs)){
		# 	st<-dynRes$cells
		# 	st$epoch[st$epoch!=epoch]<-"background"
		# 	exp<-exp[,rownames(st)]

		# 	cds<-newCountDataSet(exp,st$epoch)
		# 	cds<-estimateSizeFactors(cds)
		# 	cds<-tryCatch({estimateDispersions(cds)},error=function(e){print("using fitType=local"); estimateDispersions(cds,fitType='local')})
		# 	diffres<-nbinomTest(cds,epoch,"background")

  #   		epochs[[epoch]]<-diffres$id[diffres$padj<pThresh_DE]

		# 	chunk_cells<-rownames(dynRes$cells[dynRes$cells$epoch==epoch,])
		# 	chunk<-exp[,chunk_cells]
		# 	mean_expression<-rbind(mean_expression,data.frame(gene=rownames(chunk),epoch=rep(epoch,length(rownames(chunk))),mean_expression=rowMeans(chunk)))
		# }

		# Data is scaled and log-normalized. Use t-test here
		for (epoch in names(epochs)){
			chunk_cells<-rownames(dynRes$cells[dynRes$cells$epoch==epoch,])
			chunk<-exp[,chunk_cells]
			background<-exp[,!(colnames(exp) %in% chunk_cells)]

			diffres<-data.frame(gene=character(),mean_diff=double(),pval=double())
			for (gene in rownames(exp)){
				t<-t.test(chunk[gene,],background[gene,])
				ans<-data.frame(gene=gene,mean_diff=(t$estimate[1]-t$estimate[2]),pval=t$p.value)
				diffres<-rbind(diffres,ans)
			}

 			diffres$padj<-p.adjust(diffres$pval,method="BH")
    		diffres<-diffres[diffres$mean_diff>0,]				#Filter for genes that are on

    		epochs[[epoch]]<-diffres$gene[diffres$padj<pThresh_DE]

    		# Filter for genes that are actively expressed
    		chunk_df<-data.frame(means=rowMeans(chunk))
			chunk_df<-cbind(chunk_df,thresholds)
			chunk_df$active<-(chunk_df$means>=chunk_df$thresh)

			epochs[[epoch]]<-intersect(epochs[[epoch]],rownames(chunk_df[chunk_df$active,]))

			mean_expression<-rbind(mean_expression,data.frame(gene=rownames(chunk),epoch=rep(epoch,length(rownames(chunk))),mean_expression=rowMeans(chunk)))
		}


	}

	# some dynamically assigned genes may not be assigned to any epoch
	# if forceGenes, then assign these orphan genes to the epoch in which they have max expression
	if(forceGenes){
		assignedGenes = unique(unlist(epochs))
		orphanGenes = setdiff(rownames(exp), assignedGenes)
		cat("There are ", length(orphanGenes), " orphan genes\n")
		for(oGene in orphanGenes){
			xdat = mean_expression[mean_expression$gene==oGene,]
			ep = xdat[which.max(xdat$mean_expression),]$epoch
			epochs[[ep]] = append(epochs[[ep]], oGene)
		}
	}

	epochs$mean_expression<-mean_expression

	epochs

}

#' Assigns genes to epochs just based on which mean is maximal
#'
#' @export
#'
assign_epochs_simple<-function(expDat,dynRes,num_epochs=3,pThresh=0.01,toScale=FALSE){

	# limit exp to dynamically expressed genes
	exp<-expDat[names(dynRes$genes[dynRes$genes<pThresh]),]
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
	
	# Divide epochs by pseudotime

	# create empty list to contain active genes
	epoch_names<-paste0(rep("epoch",num_epochs),seq(1:num_epochs))
	epochs<-vector("list",length(epoch_names))
	names(epochs)<-epoch_names
	
	# determine activity based on average expression in each epoch
	ptmax<-max(dynRes$cells$pseudotime)
	ptmin<-min(dynRes$cells$pseudotime)
	chunk_size<-(ptmax-ptmin)/num_epochs

	cellsEps = rep("", length(names(t1)))
	names(cellsEps) = names(t1)

	for (i in 1:length(epochs)){
		lower_bound<-ptmin+((i-1)*chunk_size)
		upper_bound<-ptmin+(i*chunk_size)
		chunk_cells<-rownames(dynRes$cells[dynRes$cells$pseudotime>=lower_bound & dynRes$cells$pseudotime<=upper_bound,])
		chunk<-exp[,chunk_cells]

		chunk_df<-data.frame(means=rowMeans(chunk))
		chunk_df<-cbind(chunk_df,thresholds)
		chunk_df$active<-(chunk_df$means>=chunk_df$thresh)

		epochs[[i]]<-rownames(chunk_df[chunk_df$active,])
		genesPeakTimes = apply(chunk, 1, which.max)
		gpt = as.vector(dynRes$cells[chunk_cells,][genesPeakTimes,]$pseudotime)

		mean_expression<-rbind(mean_expression,data.frame(gene=rownames(chunk),epoch=rep(i,length(rownames(chunk))),mean_expression=rowMeans(chunk), peakTime = gpt))
		cellsEps[chunk_cells] = epoch_names[i]
	}

	# assign genes to epochs
	genes = unique(as.vector(mean_expression$gene))
	cat("n genes: ",length(genes),"\n")
		eps = rep("", length(genes))
		geneEpPT = rep(0, length(genes))
		epMean = rep(0, length(genes))

		names(eps) = genes
		names(geneEpPT) = genes
		names(epMean) = genes
		for(gene in genes){
		x = mean_expression[mean_expression$gene==gene,]
		xi = which.max(x$mean_expression)
		eps[[gene]] = as.vector(x[xi,]$epoch)
		geneEpPT[[gene]] = as.vector(x[xi,]$peakTime)
		epMean[[gene]] = max(x$mean_expression)
		}

		geneDF = data.frame(gene=genes, epoch=eps, peakTime = geneEpPT, epMean = epMean, pval=dynRes$genes[genes])
		cells2 = dynRes$cells[names(t1),]
		cells2 = cbind(cells2, epoch=cellsEps)
		#cells2 = 'a'

#	mean_expression$epoch<-paste0("epoch",mean_expression$epoch)
#	epochs$mean_expression<-mean_expression
#	epochs
	list(genes=geneDF, cells=cells2)
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




# =================== Old function to cluster and order genes =====================
# useful for quick and simple epoch assignment, take the place of define_epochs and assign_epochs

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




# =================== Functions to assess importance of transcription factors =====================



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


#' Function to compute betweenness and degree
#'
#' @param dynnet result of dynamic GRN reconstruction
#' @param weight_column name of column in dynnet to weight edges
#' @param directed_graph if GRN is directed or not
#'
#' @return list of dataframes of active TFs in each epoch and transition, ranked by betweenness*degree
#' @export
compute_betweenness_degree<-function(dynnet,weight_column="zscore",directed_graph=FALSE){

  ranks<-sapply(names(dynnet), function (x) NULL)
  for (net_name in names(dynnet)){
    df<-dynnet[[net_name]]
    
    df<-df[,c("TF","TG",weight_column)]
    colnames(df)<-c("TF","TG","weight")

    net<-graph_from_data_frame(df,directed=directed_graph)
    b<-betweenness(net,directed=directed_graph,normalized=TRUE)
    b<-as.data.frame(b)

    degree<-degree(net,mode="all",normalized=TRUE)
    degree<-data.frame(degree[rownames(b)])

	b<-cbind(b,degree)
	colnames(b)<-c("betweenness","degree")
	b$temp<-b$betweenness*b$degree

	b<-b[order(b$temp,decreasing=TRUE),]
	b$temp<-NULL

	b$gene<-rownames(b)
	b<-b[,c("gene","betweenness","degree")]

	b$is_regulator<-FALSE
    b$is_regulator[b$gene %in% unique(df$TF)]<-TRUE

    ranks[[net_name]]<-b
  }
  
  ranks

}












