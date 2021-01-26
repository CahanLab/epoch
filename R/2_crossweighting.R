#' compute cross correlation
#'
#' @param grnTab
#' @param expSmoothed
#' @param lag
#'
#' @return grnTab with max corr and offset added
#' 
#' @export
#'
crossweight<-function(grnTab,
					expSmoothed,
					lag=floor(ncol(expSmoothed)/5),
					min=ceiling(ncol(expSmoothed)/50),
					max=floor(ncol(expSmoothed)/12),
					filter_thresh=0){


	grnTab$TG<-as.character(grnTab$TG)
	grnTab$TF<-as.character(grnTab$TF)

	offset<-apply(grnTab,1,cross_corr,expSmoothed=expSmoothed,lag=lag)
	grnTab$offset<-offset

	weighted_score<-c()
	for (i in 1:nrow(grnTab)){
		new<-score_offset(grnTab$zscore[i],grnTab$offset[i],min=min,max=max)
		weighted_score<-c(weighted_score,new)
	}

	grnTab$weighted_score<-weighted_score

	grnTab<-grnTab[grnTab$weighted_score>filter_thresh,]

	grnTab
}


cross_corr<-function(grn_row,expSmoothed,lag){

	tg<-grn_row[1]
	tf<-grn_row[2]

	x<-ccf(expSmoothed[tf,],expSmoothed[tg,],lag,pl=FALSE)

	df<-data.frame(lag=x$lag,cor=abs(x$acf))
    df<-df[order(df$cor,decreasing=TRUE),]
    offset<-mean(df$lag[1:ceiling((2/3)*lag)])

    offset

}

score_offset<-function(score,offset,min=2,max=20){

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










