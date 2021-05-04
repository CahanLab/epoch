#' Perform crossweighting
#'
#' @param grnTab GRN dataframe, the result of running reconstructGRN or reconstructGRN_GENIE3
#' @param expDat genes-by-cells expression matrix
#' @param xdyn result of running findDynGenes
#' @param lag lag window on which to run cross-correlation. Cross-correlaiton computed from -lag to +lag.
#' @param min minimum of weighting window. Edges with offsets (or absolute offsets if symmetric_filter=TRUE) less than min will not be negatively weighted.
#' @param max maximum of weighting window. Edges with offsets (or absolute offsets if symmetric_filter=TRUE) greater than max will have weights set to 0.
#' @param symmetric_filter whether or not to employ a symmetric weight scheme. If true, absolute offset is used in place of offset.
#' @param filter_thresh after crossweighting, edges with weights less than filter_thresh will be set to 0.
#'
#' @return grnTab with offset and weighted_score added
#' 
#' @export
#'
crossweight<-function(grnTab,
					expDat,
					xdyn,
					lag=floor(ncol(expDat)/5),
					min=ceiling(ncol(expDat)/50),
					max=floor(ncol(expDat)/12),
					symmetric_filter=FALSE,
					filter_thresh=0){

	# order expDat
	expDat<-expDat[,rownames(xdyn$cells)]

	grnTab$TG<-as.character(grnTab$TG)
	grnTab$TF<-as.character(grnTab$TF)

	offset<-apply(grnTab,1,cross_corr,expDat=expDat,lag=lag)
	grnTab$offset<-offset

	weighted_score<-c()
	for (i in 1:nrow(grnTab)){
		new<-score_offset(grnTab$zscore[i],grnTab$offset[i],min=min,max=max,symmetric_filter=symmetric_filter)
		weighted_score<-c(weighted_score,new)
	}

	grnTab$weighted_score<-weighted_score

	grnTab<-grnTab[grnTab$weighted_score>filter_thresh,]

	grnTab
}


cross_corr<-function(grn_row,expDat,lag){

	tg<-grn_row[1]
	tf<-grn_row[2]

	x<-ccf(expDat[tf,],expDat[tg,],lag,pl=FALSE)

	df<-data.frame(lag=x$lag,cor=abs(x$acf))
    df<-df[order(df$cor,decreasing=TRUE),]
    offset<-mean(df$lag[1:ceiling((2/3)*lag)])

    offset

}

score_offset<-function(score,offset,min=2,max=20,symmetric_filter=FALSE){

	if(symmetric_filter){
		offset<-abs(offset)
	}

	if (offset<=min){
		res<-score
	}else if (offset>=max){
		res<-0
	}else{
		# linear weighting scheme according to y=(-x/(max-min))+1
		weight<-(-offset/(max-min))+1
		res<-score*weight
	}

	res
}


# estimates min and max values for crossweighting
# for now assumes uniform cell density across pseudotime/only considers early time
# this needs to be refined if it's to be useful...
#' @export
crossweight_params<-function(expDat,xdyn,pseudotime_min=0.005,pseudotime_max=0.01){

	expDat<-expDat[,rownames(xdyn$cells)]
	ncells<-nrow(xdyn$cells)
	min<-nrow(xdyn$cells[xdyn$cells$pseudotime<pseudotime_min,])
	max<-nrow(xdyn$cells[xdyn$cells$pseudotime<pseudotime_max,])

	params<-list(min=min,max=max)

}







