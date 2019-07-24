# warpnet
gene regulatory network reconstruction from scRNA-Seq using dynamic time warp


## Introduction
warpNet: reconstruct GRNs from scRNA-Seq data useing trajectory inference and dynamic time warp.

## Data

This data is unpublished data from our lab. It is of muscle development at e12.5. It has already been normalized, and the varying genes have been identified. It has also been clustered, and analyzed with RNA Velocity (Velocyto). Here is what the data look like after applying diffusion map:

<img src="img/musclePCA_072219.png">

And here is RNA Velocity result:

<img src="img/muscleVelocyto_072219.png">

And here is the ordering of the cells based on diffusion pseudotime (dpt):

<img src="img/muscleDPT_072219.png">

Now, let's use warpNet to reconstruct the GRNs that underpin this trajectory.


### Walk thru

## Set up
```R

library(igraph)
library(qgraph)
library(loomR)
library(gam)

library(devtools)
install_github("pcahan1/singleCellNet")

library(singleCellNet)

install_github("pcahan1/warpnet")
library(warpnet)

library(minet)

mydate<-utils_myDate()

```

## Load data
```R
# Data is in the R package data folder
pathToWarpNet = "../"
mmTFs<-utils_loadObject( paste0(pathToWarpNet, "data/mmTFs.rda") )

list12<-loadLoomExpUMAP( paste0(pathToWarpNet, "data/adMuscle_E12_DPT_071919.loom", xname='leiden', has_dpt_groups=FALSE)
expDat<-list12[['expDat']]
sampTab<-list12[['sampTab']] # holds the clustering and pca/umap info

grps<-as.vector(sampTab$cluster)
names(grps)<-as.vector(sampTab$cell_name)
```

## Find dynamically expressed genes
```R

# PCA was done in scanpy/python and stored in sampTab
system.time(xdyn <- findDynGenes(expDat, sampTab, c("0","1")))

starting gammma...
   user  system elapsed 
 14.195   1.693  15.972
```

## Smooth expression
```R
ccells = xdyn$cells
system.time(expSmoothed <- grnKsmooth(expDat, ccells))
  user  system elapsed 
  0.832   0.013   0.850 

 dim(expSmoothed)
[1] 2074  262
```

## Cluster genes into epochs
```R

geneDF = caoGenes(expSmoothed, xdyn, k=3, pThresh=0.01, method='kmeans')
gdfForHM = as.data.frame(geneDF[,"epoch"])
rownames(gdfForHM) = rownames(geneDF)
hm_dyn_clust(testSm, xdyn, geneAnn= gdfForHM, toScale=TRUE)
```

<img src="img/heatmapDynGenes_072219.png">

## Reconstruct GRN
 NB: only using dynamically expressed genes for GRN
```R
system.time(grnDF <- reconstructGRN(expSmoothed[rownames(geneDF),], mmTFs, zThresh=3))
   user  system elapsed 
  0.215   0.023   0.238

dim(grnDF)
[1] 1803    4
```

## Add PT-based weight to GRN

normDist is difference in peakTimes between TF and TG. Negative if the sign of the difference is inconsistent with predicted regulatory influence based on CLR network (from pearson)
```R
system.time(grnDF <- addDistWeight(grnDF, geneDF, ncells = ncol(expSmoothed), maxNeg= -0.2))
   user  system elapsed 
  0.183   0.004   0.186 

 grnDF[1:3,]
    TG    TF   zscore      corr adjWeight   normDist
1 Eya1 Meis1 3.022976 0.9911989  2.692854 0.09809160
2 Eya1  Pcna 3.245824 0.9913648  3.075591 0.08282443
3  Msc Dnmt1 3.006199 0.9881022  2.647884 0.10038168

```

## Need a function that will give a lot of information for each TF that can be used as a basis for better and more diverse way of selecting TFs
Here is a start to do so
```R
tfTab = evalTFs(grnDF, geneDF)

 tfTab[order(tfTab$weightMean, decreasing=T),][1:10,]
                 TF num_TFs epoch distToStart weightTotal weightMean
Zeb1           Zeb1      15     2           3    76.14665   3.626031
Hist1h2ak Hist1h2ak      36     1         199    58.97373   3.103880
Hist1h2ae Hist1h2ae      36     1         193    20.54976   2.935680
Myod1         Myod1      15     2          82   115.62926   2.890731
Smyd1         Smyd1      25     3         202   115.27090   2.881773
H2afx         H2afx      36     1         191    44.35505   2.772191
Hmgb2         Hmgb2      36     1         177    65.87164   2.744652
Mef2a         Mef2a      15     2         106    91.15749   2.604500
Six2           Six2      36     1          99    20.82408   2.603010
Prox1         Prox1      25     3         141    27.47881   2.498073
```


## Score TFs per epoch

```R
 pickExemplars(  geneDF,grnDF, topX=5, type='weight')
$`1`
[1] "H2afz" "Pax7"  "Meis1" "Hmgb2" "Msc"  

$`2`
[1] "Myod1"  "Zbtb18" "Mef2a"  "Zeb1"   "Tmf1"  

$`3`
[1] "Smyd1" "Foxd3" "Rxrg"  "Mef2c" "Klf5" 
```


## Plot these + top 5 positive regulons
```R
iG_x = ig_exemplars(grnDF, geneDF, topTFs, topX=5) 
x2 = ig_convertMedium(iG_x, vScale=2)
e = get.edgelist(x2, names=FALSE)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(x2))
plot(x2, layout=l)
```

<img src="img/smallGRN_072419.png">

