# mESC directed differentiation GRN reconstruction


## Introduction
Here we apply Epoch to reconstruct the dynamic GRN underlying in vitro directed differentiation of mESCs and in vivo gastrulation.

For a more in depth look at Epoch, and to read about how we applied it to elucidate signaling-induced GRN topology changes in early mouse embryonic stem cell (ESC) directed differentiation, check out our [preprint here](https://www.biorxiv.org/content/10.1101/2021.05.06.443021v2).


1. [Dynamic GRN Reconstruction](#reconstruction)
		+ [Mesoderm network reconstruction](#mesoderm)
		+ [Endoderm network reconstruction](#endoderm)
		+ [Neuroectoderm network reconstruction](#neuroectoderm)
		+ [in vivo network reconstruction](#invivo)


## Installing Epoch
```R

devtools::install_github("pcahan1/epoch")
library(epoch)

```

## Dynamic GRN Reconstruction with Epoch <a name="reconstruction"></a>

##### The data
```R
early_diff <- loadDataFromLoom("20191205multiseq_scvelo_20201020.loom")

early_diff_exp<-early_diff[['expDat']]
early_diff_st<-early_diff[['sampTab']]

head(early_diff_st[,c("timepoint","treatment","leiden_refined","latent_time")])

#                  timepoint treatment leiden_refined latent_time
# AAACCCAAGGTGATCG         3       WAG              0  0.22996202
# AAACCCAAGTCACAGG         2       pre              3  0.12755865
# AAACCCAAGTGCAGGT         3        WA            0,4  0.08814099
# AAACCCACACTGAATC         4       WAB              2  0.58602914
# AAACCCACATCGGTTA         0       pre              1  0.14227916
# AAACCCATCAATCCAG         4       WAB              3  0.33609538

mmTFs<-utils_loadObject("data/mmTFs.rda")
mmTFs<-intersect(rownames(early_diff_exp),mmTFs)

```
The sample table (early_diff_st) contains cluster annotation (specified here as "leiden_refined") and pseudotime annotation (specified here as "latent_time").


##### Reconstruct the mesoderm network <a name="mesoderm"></a>
  
Begin by limiting the data to the mesoderm lineage.
```R
p1_st<-early_diff_st[early_diff_st$leiden_refined %in% c('1','3','0','8'),]
p1_exp<-early_diff_exp[,rownames(p1_st)]

```
  
Find dynamic genes and reconstruct
```R
# Find dynamic genes
p1dyn<-findDynGenes(p1_exp, p1_st, c('1','3','0','8'),group_column="leiden_refined",pseudotime_column="latent_time")

pThresh<-0.05
p1dgenes<-names(p1dyn$genes)[p1dyn$genes<pThresh]

# Reconstruct and perform crossweighting. For now, we will keep threshold at 0 (i.e. don't filter)
p1grnDF<-reconstructGRN(p1_exp,mmTFs,p1dgenes,method="MI",zThresh=0)
p1grnDF<-crossweight(p1grnDF,p1_exp,p1dyn,filter_thresh=0)

```
  
Smooth expression and plot a nice heatmap
```R
p1smooth<-grnKsmooth(p1_exp, p1dyn$cells, BW=.05)

plot_heatmap<-function(xdyn,expSmoothed,tfs){
    tfstoplot<-intersect(names(xdyn$genes[xdyn$genes<pThresh]),tfs)
    dynTFs<-xdyn
    dynTFs$genes<-dynTFs$genes[names(dynTFs$genes) %in% tfstoplot]
    hm_dyn(expSmoothed,dynTFs,topX=500, toScale=T)
}

plot_heatmap(p1dyn,p1smooth,mmTFs)

```
  
Based on heatmap, define epochs and extract the dynamic network
```R
p1dyn<-define_epochs(p1dyn,p1smooth,method="pseudotime",num_epochs=3,pseudotime_cuts=c(.15,.25))
p1epochs<-assign_epochs(p1_exp,p1dyn,method="DE")
p1dynamic_grn<-epochGRN(p1grnDF, p1epochs)

```
  
p1dynamic_grn now contains the dynamic mesoderm network; The same process and parameters were used to reconstruct treatment-specific networks, but limited as following:
* WAG-treatment network: reconstruction limited to "pre" and "WAG" treated cells
* WAB-treatment network: reconstruction limited to "pre" and "WAB" treated cells
* WAN-treatment network: reconstruction limited to "pre" and "WAN" treated cells
* WA-treatment network: reconstruction limited to "pre" and "WA" treated cells
  


##### Reconstruct the endoderm network <a name="endoderm"></a>
  
Begin by limiting the data to the endoderm lineage.
```R
p2_st<-early_diff_st[early_diff_st$leiden_refined %in% c('1','3','0','7,1','7,0'),]
p2_exp<-early_diff_exp[,rownames(p2_st)]

```
  
Find dynamic genes and reconstruct
```R
# Find dynamic genes
p2dyn<-findDynGenes(p2_exp, p2_st, c('1','3','0','7,1','7,0'),group_column="leiden_refined",pseudotime_column="latent_time")

pThresh<-0.05
p2dgenes<-names(p2dyn$genes)[p2dyn$genes<pThresh]

# Reconstruct and perform crossweighting. For now, we will keep threshold at 0 (i.e. don't filter)
p2grnDF<-reconstructGRN(p2_exp,mmTFs,p2dgenes,method="MI",zThresh=0)
p2grnDF<-crossweight(p2grnDF,p2_exp,p2dyn,filter_thresh=0)

```
  
Smooth expression and plot a nice heatmap
```R
p2smooth<-grnKsmooth(p2_exp, p2dyn$cells, BW=.05)

plot_heatmap<-function(xdyn,expSmoothed,tfs){
    tfstoplot<-intersect(names(xdyn$genes[xdyn$genes<pThresh]),tfs)
    dynTFs<-xdyn
    dynTFs$genes<-dynTFs$genes[names(dynTFs$genes) %in% tfstoplot]
    hm_dyn(expSmoothed,dynTFs,topX=500, toScale=T)
}

plot_heatmap(p2dyn,p2smooth,mmTFs)

```
  
Based on heatmap, define epochs and extract the dynamic network
```R
p2dyn<-define_epochs(p2dyn,p2smooth,method="pseudotime",num_epochs=3,pseudotime_cuts=c(.15,.25))
p2epochs<-assign_epochs(p2_exp,p2dyn,method="DE")
p2dynamic_grn<-epochGRN(p2grnDF, p2epochs)

```
  
p2dynamic_grn now contains the dynamic endoderm network. ; The same process and parameters were used to reconstruct treatment-specific networks, but limited as following:
* WAG-treatment network: reconstruction limited to "pre" and "WAG" treated cells
* WAB-treatment network: reconstruction limited to "pre" and "WAB" treated cells
* WAN-treatment network: reconstruction limited to "pre" and "WAN" treated cells
* WA-treatment network: reconstruction limited to "pre" and "WA" treated cells
  


##### Reconstruct the neuroectoderm network <a name="neuroectoderm"></a>
  
Begin by limiting the data to the neuroectoderm lineage.
```R
p3_st<-early_diff_st[early_diff_st$leiden_refined %in% c('1','3','0','0,4','0,2','0,3','7,2','2','4','5','9'),]
p3_exp<-early_diff_exp[,rownames(p3_st)]

```
  
Find dynamic genes and reconstruct
```R
# Find dynamic genes
p3dyn<-findDynGenes(p3_exp, p3_st, c('1','3','0','0,4','0,2','0,3','7,2','2','4','5','9'),group_column="leiden_refined",pseudotime_column="latent_time")

pThresh<-0.05
p3dgenes<-names(p3dyn$genes)[p3dyn$genes<pThresh]

# Reconstruct and perform crossweighting. For now, we will keep threshold at 0 (i.e. don't filter)
p3grnDF<-reconstructGRN(p3_exp,mmTFs,p3dgenes,method="MI",zThresh=0)
p3grnDF<-crossweight(p3grnDF,p3_exp,p3dyn,filter_thresh=0)

```
  
Smooth expression and plot a nice heatmap
```R
p3smooth<-grnKsmooth(p3_exp, p3dyn$cells, BW=.05)

plot_heatmap<-function(xdyn,expSmoothed,tfs){
    tfstoplot<-intersect(names(xdyn$genes[xdyn$genes<pThresh]),tfs)
    dynTFs<-xdyn
    dynTFs$genes<-dynTFs$genes[names(dynTFs$genes) %in% tfstoplot]
    hm_dyn(expSmoothed,dynTFs,topX=500, toScale=T)
}

plot_heatmap(p3dyn,p3smooth,mmTFs)

```
  
Based on heatmap, define epochs and extract the dynamic network
```R
p3dyn<-define_epochs(p3dyn,p3smooth,method="pseudotime",num_epochs=3,pseudotime_cuts=c(.15,.25))
p3epochs<-assign_epochs(p3_exp,p3dyn,method="DE")
p3dynamic_grn<-epochGRN(p3grnDF, p3epochs)

```
  
p3dynamic_grn now contains the dynamic neuroectoderm network.
  

##### Reconstruct the in vivo mesoderm network <a name="invivo"></a>
  
The data (limited to epiblast, primitive streak early, primitive streak late, and mesoderm presomitic)
```R
list12<-loadDataFromLoom("adGrosswendt_trimmed250_dpt_20210817.loom")
expDat<-list12[['expDat']]
sampTab<-list12[['sampTab']]

head(sampTab[,c("celltype","dpt_pseudotime")])

#                             celltype dpt_pseudotime
# WT65;SG4_AAACCTGCACCAGGCT-1 Epiblast    0.008824126
# WT65;SG4_AAACCTGTCAACACCA-1 Epiblast    0.081007242
# WT65;SG4_AAAGCAATCTTAGCCC-1 Epiblast    0.095422700
# WT65;SG4_AACACGTTCTACTATC-1 Endoderm    0.472867429
# WT65;SG4_AACCATGCACAACGTT-1 Epiblast    0.015073729
# WT65;SG4_AACTCAGTCGAATCCA-1 Epiblast    0.105524026

sampTab<-sampTab[sampTab$celltype %in% c("Epiblast","Primitive streak early","Primitive streak late","Mesoderm presomitic"),]

expDat<-expDat[,rownames(sampTab)]
expDat<-expDat[!(rowSums(expDat)==0),]

```

Network reconstruction
```R
xdyn <- findDynGenes(expDat, sampTab, group_column="celltype", pseudotime_column="dpt_pseudotime")
pThresh<-0.05
dgenes<-names(xdyn$genes)[xdyn$genes<pThresh]

grnDF <- reconstructGRN(expDat, mmTFs, dgenes, method="MI", zThresh=0)
grnDF <- crossweight(grnDF,expDat,xdyn,filter_thresh=0)

ccells <- xdyn$cells
expSmoothed <- grnKsmooth(expDat, ccells, BW=0.02)

# Plot a heatmap of the dynamic TFs
tfstoplot<-intersect(dgenes,mmTFs)
dynTFs<-xdyn
dynTFs$genes<-dynTFs$genes[names(dynTFs$genes) %in% tfstoplot]
hm_dyn(expSmoothed,dynTFs,topX=50)

genesToPlot<-sort(names(xdyn$genes[xdyn$genes<pThresh]))
dynGenes<-xdyn
dynGenes$genes<-dynGenes$genes[genesToPlot] 
hm_dyn(expSmoothed,dynGenes,topX=500, toScale=T)

# define epochs and extract dynamic GRN
xdyn<-define_epochs(xdyn,expDat,method="pseudotime",num_epochs=3,pseudotime_cuts=c(.2,.37))
epoch_assignments<-assign_epochs(expDat,xdyn)
dynamic_grn<-epochGRN(grnDF,epoch_assignments)

```














