####------------------------------------------------------------------------------------------------------
# this is used for reproducibility of fig4 and corresponding supplementary figs and tables
# scavenge analysis of Acute lymphoblastic leukemia trait using a hematopoiesis scATAC dataset
####------------------------------------------------------------------------------------------------------
library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(tidyverse)
library(SummarizedExperiment)
library(Matrix)
library(BuenColors)
library(cowplot)
library(diffloop)
library(igraph)
library(stringr)
setwd("/Users/fyu/Documents/GitHub/SCAVENGE-reproducibility/data")





load("Satpathy_2019_heme_ALL_gscore.rda")
load("Satpathy_2019_heme_knngraph.rda")
colnames(mutualknn30) <- rownames(mutualknn30) <- zscoreWeighted2$cell_name
zscoreWeighted2$cell_cluster2 <- zscoreWeighted2$cell_cluster

baso_df <- zscoreWeighted2
trait_name="ALL"
seedcutoff=0.05
workdir="/Users/fyu/Documents/GitHub/SCAVENGE-reproducibility/data/all"
umap_pointsize=0.1
umap_pointalpha= 0.6
umap_color="wolfgang_extra"

setwd(workdir)
original_df <- baso_df

if( any(is.na(baso_df$Zscore))){
    message("NA occurs in the zscore")
    message("How many? [", sum(is.na(baso_df$Zscore)), "]")
    message("set NA to 0")
    baso_df$Zscore[is.na(baso_df$Zscore)] <- 0
}

zscore<- baso_df$Zscore
seed_idx <- seedindex(baso_df$Zscore, 0.05)
scale_factor <- cal_scalefactor(z_score=baso_df$Zscore, 0.01)

### Network propagation  
np_score <- randomWalk_sparse(intM=mutualknn30, queryGenes=rownames(mutualknn30)[seed_idx], gamma=0.05)
# **Trait relevant score (TRS) with scaled and normalized**  
# A few cells are singletons are removed from further analysis
omit_idx <- np_score==0
sum(omit_idx)
mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
np_score <- np_score[!omit_idx]
TRS <- np_score %>% capOutlierQuantile(., 0.99) %>% max_min_scale
TRS <- TRS * scale_factor
zscoreWeighted3 <- data.frame(zscoreWeighted2[!omit_idx, ], TRS)

write.csv(zscoreWeighted3, "covid19b1_table.csv", row.names=F)
