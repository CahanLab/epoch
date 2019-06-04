getGRNMat = function(expressionMatrix, movingInterval, delay, threshold){
  tree.phate <- phate(expressionMatrix);
  phenotype=data.frame(row.names=row.names(expressionMatrix), pseudo=row.names(expressionMatrix));
  phenotype$pseudo=tree.phate$embedding[,1];
  phenotype$pseudo2=tree.phate$embedding[,2];
  #plot(tree.phate, col = tree.data$branches) #for visual inspection 
  phenotype=phenotype[row.names(phenotype)[order(phenotype$pseudo)],];
  expressionMatrix=expressionMatrix[row.names(phenotype),];
  #smoothing by moving average
  expressionMatrixSmoothed=data.frame(expressionMatrix);
  expressionMatrixSmoothed[expressionMatrixSmoothed==0]=NA;
  for (i in nrow(expressionMatrix):movingInterval){
    expressionMatrixSmoothed[i,]=colMeans(expressionMatrixSmoothed[i:(i-movingInterval),], na.rm = TRUE)
  }
  expressionMatrixSmoothed=expressionMatrixSmoothed[movingInterval:nrow(expressionMatrixSmoothed),];
  scaled=scale(expressionMatrixSmoothed);
  scaled[is.na(scaled)]=0;
  activation=matrix(nrow=ncol(scaled), ncol=ncol(scaled));
  deactivation=matrix(nrow=ncol(scaled), ncol=ncol(scaled));
  for (i in 1:nrow(activation)){
    for (j in 1:ncol(activation)){
      tf=scaled[,i];
      target=scaled[,j]
      tf=filter.fft(tf,BW=0.01, n=10)
      target=filter.fft(target,BW=0.01, n=10)
      activation[i,j]=dtw(Re(tf)[1:(length(tf)-delay)], Re(target)[(delay+1):length(target)])$normalizedDistance;
      deactivation[i,j]=dtw(Re(tf)[1:(length(tf)-delay)], -1*Re(target)[(delay+1):length(target)])$normalizedDistance;
    }
    print(i)
  }
  blah=list(activation, deactivation)
  row.names(activation)=colnames(expressionMatrix)
  colnames(activation)=colnames(expressionMatrix)
  bool=activation;
    for (i in 1:nrow(activation)){
      for (j in 1:ncol(activation)){
        if((activation[i,j])<threshold){
          bool[i,j]=1
        }
        else if((deactivation[i,j])<threshold){
          bool[i,j]=-1
        }
        else {
          bool[i,j]=0
        }
      }
    testGRN=melt(as.matrix(bool))
    testGRN=testGRN[testGRN$value!=0,]
    testGRN$value[testGRN$value==-1]="Inhibition"
    testGRN$value[testGRN$value==1]="Activation"
    testGRN=testGRN[, c(1 ,3, 2)]
    colnames(testGRN)=c("TF", "Interaction", "Target")
    testGRN
  }
}

grnSSF = function(
  expMat, # rows are cells
  cellOrder, # as provided by TI
  movingInterval = 5,
  BW = 0.01,
  fftN = 10){
    
    genes = colnames(expMat)
    expMat = expMat[cellOrder,]

    # smoothing by moving average
    expressionMatrixSmoothed=data.frame(expMat)
    colnames(expressionMatrixSmoothed) = genes
    expressionMatrixSmoothed[expressionMatrixSmoothed==0]=NA
    for (i in nrow(expressionMatrix):movingInterval){
      expressionMatrixSmoothed[i,]=colMeans(expressionMatrixSmoothed[i:(i-movingInterval),], na.rm = TRUE)
    }
    expressionMatrixSmoothed=expressionMatrixSmoothed[movingInterval:nrow(expressionMatrixSmoothed),]
    scaled=scale(expressionMatrixSmoothed)
    scaled[is.na(scaled)]=0

    colnames(scaled) = genes
    cat("FFT\n")
    x = apply(scaled, 2, filter.fft, BW=BW, n = fftN)
    colnames(x) = genes
    list(scaled=scaled, ffTed = Re(x))
}


grnDTWx = function(grnFF_res, tfs, delay=5){
 

  ncells = nrow(grnFF_res)
  targets = colnames(grnFF_res)
  ntargets = length(targets)

  tfs = intersect(tfs, targets)
  ntfs = length(tfs)

  activation   = matrix(nrow=ntargets, ncol=ntfs)
  deactivation = matrix(nrow=ntargets, ncol=ntfs)
  row.names(activation) = targets
  colnames(activation)  = tfs

  for (i in 1:ntargets){
    #  cat(targets[i])
      for (j in 1:ntfs){
    #    cat(   " ",tfs[j],"\n")
        activation[i,j]   = dtw( grnFF_res[,j][1:(ncells-delay)], grnFF_res[,i][(delay+1):ncells] )$normalizedDistance
        deactivation[i,j] = dtw( grnFF_res[,j][1:(ncells-delay)], -1*grnFF_res[,i][(delay+1):ncells] )$normalizedDistance
        
      }
  }

  list(act = activation, deact = deactivation)
}


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



