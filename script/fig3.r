####------------------------------------------------------------------------------------------------------
# this is used for reproducibility of fig4 and corresponding supplementary figs and tables
# scavenge analysis of covid19-b1 traits using a PBMC scATAC dataset
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


load("covidb1_you_2021_covidpbmc_gscore.rda")
load("you_2021_covidpbmc_knngraph.rda")
colnames(mutualknn30) <- rownames(mutualknn30) <- zscoreWeighted2$cell_name
zscoreWeighted2$cell_cluster2 <- zscoreWeighted2$cell_cluster

baso_df <- zscoreWeighted2
trait_name="covidb1"
seedcutoff=0.05
workdir="/Users/fyu/Documents/GitHub/SCAVENGE-reproducibility/data/covid19b1"
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


baso_df$seed_p0.05_top <- 

baso_df$gchromVAR_pvalue <- pnorm(baso_df$Zscore, lower.tail = FALSE)

head(baso_df$gchromVAR_pvalue)
baso_df$seed_p0.05_top <- baso_df$gchromVAR_pvalue<0.05 # 16899
message("top cells")
# how many seed you have selected
message("how many seed you have selected: ", sum(baso_df$seed_p0.05_top))
percent5_cutoff <- floor(0.05*length(baso_df$seed_p0.05_top))
message("top 5% cell number: ", percent5_cutoff)
percent1_cutoff <- floor(0.01*length(baso_df$seed_p0.05_top))
zfactor_cutoff_num <- zfactor_cutoff

percent5_cutoff_value <- sort(baso_df$gchromVAR_pvalue)[percent5_cutoff]
percent1_cutoff_value <- sort(baso_df$gchromVAR_pvalue)[percent1_cutoff]
zfactor_cutoff_num_value <- sort(baso_df$gchromVAR_pvalue)[zfactor_cutoff_num]
message("does the seeds more than 5%?: ", sum(baso_df$seed_p0.05_top) > percent5_cutoff)
if (sum(baso_df$seed_p0.05_top) > percent5_cutoff){
        message("we will use top 5% cells as the seeds")
        baso_df$seed_p0.05_top <- FALSE
        baso_df$seed_p0.05_top[baso_df$gchromVAR_pvalue <= percent5_cutoff_value] <- TRUE
} else {
        message("we will use cells with p < 0.05 as the seeds")
}

# Z score factor 
if (zfactor_cutoff_num_value > percent1_cutoff_value){
        message("we will use top 1% cells to calculate z facotr score")
        zfactor_seed_p0.05_top <- mean(baso_df$Zscore[baso_df$gchromVAR_pvalue <= percent1_cutoff_value])
} else {
        message("we will use top ",  zfactor_cutoff," cells to calculate z facotr score")
        zfactor_seed_p0.05_top <- mean(baso_df$Zscore[baso_df$gchromVAR_pvalue <= zfactor_cutoff_num_value])
}

message("z factor calculated with top1% cells: ", zfactor_seed_p0.05_top)

if(zfactor_seed_p0.05_top<0){
        message("[warning] your z factor less than 0")
}
seed_p0.05_top <- baso_df$seed_p0.05_top
gchromVAR_pvalue <- baso_df$gchromVAR_pvalue

## scatterplot top seed=p0.05 gamma=1 ------------
original_seed_p0.05_top_gamma1 <- randomWalk_sparse(intM=mutualknn30, queryGenes=(1:nrow(mutualknn30))[baso_df$seed_p0.05_top], gamma=1)
p <- ggplot(data=baso_df, aes(UMAP1, UMAP2, color=original_seed_p0.05_top_gamma1)) + geom_point(size=0.2, na.rm = TRUE, alpha = 0.4) + 
scale_color_gradientn(colors = viridis) + scale_alpha()+
    pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "_scatterplot_original_seed_p0.05_top.pdf"), width = 5, height = 5)

## scatterplot top seed=p0.05 gamma=0.05  trs_cap_0.99------------
original_seed_p0.05_top_gamma005 <- randomWalk_sparse(intM=mutualknn30, queryGenes=(1:nrow(mutualknn30))[baso_df$seed_p0.05_top], gamma=0.05)
original_seed_p0.05_top_gamma005_ceiling099 <- original_seed_p0.05_top_gamma005 %>% capOutlierQuantile(., 0.99)
p <- ggplot(data=baso_df, aes(UMAP1, UMAP2, color=original_seed_p0.05_top_gamma005_ceiling099)) + geom_point(size=0.2, na.rm = TRUE, alpha = 0.4) + 
    scale_color_gradientn(colors = viridis) + scale_alpha()+
    pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "_scatterplot_original_seed_p0.05_top_walkgamma005_ceiling099_viridis.pdf"), width = 5, height = 5)

## violinplot top seed=p0.05 gamma=0.05  trs_cap_0.99------------
pp <- ggplot(data=baso_df,  aes(x=cell_cluster2, y=original_seed_p0.05_top_gamma005_ceiling099))  +
    geom_violin(adjust=1, aes(fill=cell_cluster2)) +scale_fill_manual(values=cols)+
    geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
    ylab("Trait relevant score (TRS)") + xlab("") + ggtitle(trait_name) +
    ylim(c(0,quantile(baso_df[, original_seed_p0.05_top_gamma005_ceiling099], 1, na.rm=TRUE))) +
    theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "_violinplot_original_seed_p0.05_top_gamma005_ceiling099.pdf"), width = 10, height = 4)


original_seed_p0.05_top_gamma005_ceiling099_maxmin <- original_seed_p0.05_top_gamma005_ceiling099 %>% max_min_scale
original_seed_p0.05_top_gamma005_ceiling099_maxmin_zfactor <- original_seed_p0.05_top_gamma005_ceiling099_maxmin * zfactor_seed_p0.05_top
baso_df$original_seed_p0.05_top_gamma005_ceiling099_maxmin <- original_seed_p0.05_top_gamma005_ceiling099_maxmin
baso_df$original_seed_p0.05_top_gamma005_ceiling099_maxmin_zfactor <- original_seed_p0.05_top_gamma005_ceiling099_maxmin_zfactor
## violinplot top seed=p0.05 gamma=0.05  trs_cap_0.99_maxmin------------
pp <- ggplot(data=baso_df,  aes(x=cell_cluster2, y=original_seed_p0.05_top_gamma005_ceiling099_maxmin))  +
    geom_violin(adjust=1, aes(fill=cell_cluster2)) +scale_fill_manual(values=cols)+
    geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
    ylab("scaled TRS") + xlab("") + ggtitle(trait_name) +
    ylim(c(0, 1)) +
    theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "_violinplot_original_seed_p0.05_top_gamma005_ceiling099_maxmin.pdf"), width = 10, height = 4)
## violinplot top seed=p0.05 gamma=0.05  trs_cap_0.99_maxmin_zscore------------
pp <- ggplot(data=baso_df,  aes(x=cell_cluster2, y=original_seed_p0.05_top_gamma005_ceiling099_maxmin_zfactor))  +
    geom_violin(adjust=1, aes(fill=cell_cluster2)) +scale_fill_manual(values=cols)+
    geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
    ylab("z-factor TRS") + xlab("") + ggtitle(trait_name) +
    theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "_violinplot_original_seed_p0.05_top_gamma005_ceiling099_maxmin_zfactor.pdf"), width = 10, height = 4)

message("bottom cells")
# --------------------- bottom cells
seed_p0.05_bottom <- rank(baso_df$Zscore) <= sum(baso_df$seed_p0.05)
baso_df$seed_p0.05_bottom <- seed_p0.05_bottom
zfactor_seed_p0.05_bottom <- mean(baso_df$Zscore[baso_df$seed_p0.05_bottom])
message("z factor: ", zfactor_seed_p0.05_bottom)
## scatterplot bottom seed=p0.05 gamma=1 ------------
original_seed_p0.05_bottom_gamma1 <- randomWalk_sparse(intM=mutualknn30, queryGenes=(1:nrow(mutualknn30))[baso_df$seed_p0.05_bottom], gamma=1)
p <- ggplot(data=baso_df, aes(UMAP1, UMAP2, color=original_seed_p0.05_bottom_gamma1)) + geom_point(size=0.2, na.rm = TRUE, alpha = 0.4) + 
scale_color_gradientn(colors = viridis) + scale_alpha()+
    pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "_scatterplot_original_seed_p0.05_bottom.pdf"), width = 5, height = 5)

## scatterplot bottom seed=p0.05 gamma=0.05  trs_cap_0.99------------
original_seed_p0.05_bottom_gamma005 <- randomWalk_sparse(intM=mutualknn30, queryGenes=(1:nrow(mutualknn30))[baso_df$seed_p0.05_bottom], gamma=0.05)
original_seed_p0.05_bottom_gamma005_ceiling099 <- original_seed_p0.05_bottom_gamma005 %>% capOutlierQuantile(., 0.99)
p <- ggplot(data=baso_df, aes(UMAP1, UMAP2, color=original_seed_p0.05_bottom_gamma005_ceiling099)) + geom_point(size=0.2, na.rm = TRUE, alpha = 0.4) + 
    scale_color_gradientn(colors = viridis) + scale_alpha()+
    pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "_scatterplot_original_seed_p0.05_bottom_walkgamma005_ceiling099_viridis.pdf"), width = 5, height = 5)

## violinplot bottom seed=p0.05 gamma=0.05  trs_cap_0.99------------
pp <- ggplot(data=baso_df,  aes(x=cell_cluster2, y=original_seed_p0.05_bottom_gamma005_ceiling099))  +
    geom_violin(adjust=1, aes(fill=cell_cluster2)) +scale_fill_manual(values=cols)+
    geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
    ylab("Trait relevant score (TRS)") + xlab("") + ggtitle(trait_name) +
    ylim(c(0,quantile(baso_df[, original_seed_p0.05_bottom_gamma005_ceiling099], 1, na.rm=TRUE))) +
    theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "_violinplot_original_seed_p0.05_bottom_gamma005_ceiling099.pdf"), width = 10, height = 4)

original_seed_p0.05_bottom_gamma005_ceiling099_maxmin <- original_seed_p0.05_bottom_gamma005_ceiling099 %>% max_min_scale
original_seed_p0.05_bottom_gamma005_ceiling099_maxmin_zfactor <- original_seed_p0.05_bottom_gamma005_ceiling099_maxmin * zfactor_seed_p0.05_bottom
baso_df$original_seed_p0.05_bottom_gamma005_ceiling099_maxmin <- original_seed_p0.05_bottom_gamma005_ceiling099_maxmin
baso_df$original_seed_p0.05_bottom_gamma005_ceiling099_maxmin_zfactor <- original_seed_p0.05_bottom_gamma005_ceiling099_maxmin_zfactor
## violinplot bottom seed=p0.05 gamma=0.05  trs_cap_0.99_maxmin------------
pp <- ggplot(data=baso_df,  aes(x=cell_cluster2, y=original_seed_p0.05_bottom_gamma005_ceiling099_maxmin))  +
    geom_violin(adjust=1, aes(fill=cell_cluster2)) +scale_fill_manual(values=cols)+
    geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
    ylab("scaled TRS") + xlab("") + ggtitle(trait_name) +
    ylim(c(0, 1)) +
    theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "_violinplot_original_seed_p0.05_bottom_gamma005_ceiling099_maxmin.pdf"), width = 10, height = 4)
## violinplot bottom seed=p0.05 gamma=0.05  trs_cap_0.99_maxmin_zscore------------
pp <- ggplot(data=baso_df,  aes(x=cell_cluster2, y=original_seed_p0.05_bottom_gamma005_ceiling099_maxmin_zfactor))  +
    geom_violin(adjust=1, aes(fill=cell_cluster2)) +scale_fill_manual(values=cols)+
    geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
    ylab("z-factor TRS") + xlab("") + ggtitle(trait_name) +
    theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "_violinplot_original_seed_p0.05_bottom_gamma005_ceiling099_maxmin_zfactor.pdf"), width = 10, height = 4)


temp_df <- data.frame(original_df, gchromVAR_pvalue, seed_p0.05_top, original_seed_p0.05_top_gamma005, original_seed_p0.05_top_gamma005_ceiling099, original_seed_p0.05_top_gamma005_ceiling099_maxmin, original_seed_p0.05_top_gamma005_ceiling099_maxmin_zfactor, seed_p0.05_bottom, original_seed_p0.05_bottom_gamma005, original_seed_p0.05_bottom_gamma005_ceiling099, original_seed_p0.05_bottom_gamma005_ceiling099_maxmin, original_seed_p0.05_bottom_gamma005_ceiling099_maxmin_zfactor)
return(temp_df)
}


covidb1_df <- zscoreWeighted2
covidb1_trs <- trait_netpropagation(trait_df=covidb1_df, trait_name="covidb1", seedcutoff=0.05, workdir="/Users/fyu/Documents/GitHub/SCAVENGE-reproducibility/data/covid19b1", umap_pointsize=0.1, umap_pointalpha= 0.6, umap_color="wolfgang_extra")
covidb1_trs <- trait_netpropagation(trait_df=covidb1_df, trait_name="covidb1", seedcutoff=0.05, workdir="/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/covidb1_you_2021_covidpbmc/covidb1-5", umap_pointsize=0.1, umap_pointalpha= 0.6, umap_color="blackbright")


identical(covidb1_trs$cell_cluster2, meta_df$cell_cluster)
# [1] TRUE

# p <- ggplot(data=covidb1_df, aes(UMAP1, UMAP2, color=cell_cluster2)) + geom_point(size=0.2, na.rm = TRUE, alpha = 0.4) + 
#     scale_color_discrete() +
#      pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
# ggsave(p, file = "/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/covidb1_you_2021_covidpbmc/cellcluter.png", width = 8, height = 5)

colnames(covidb1_trs)
#  [1] "cell_name"                   "cell_cluster"
#  [3] "UMAP1"                       "UMAP2"
#  [5] "trait"                       "Zscore"
#  [7] "Health_state"                "Sample_state"
#  [9] "cell_cluster2"               "gchromVAR_pvalue"
# [11] "seed_p0.05_top"              "seed_p0.05_random"
# [13] "seed_p0.05_bottom"           "Zscore_scaled"
# [15] "zscore_topseed"              "zscore_bottomseed"
# [17] "zscore_randomseed"           "zscore_topseed_np"
# [19] "zscore_bottomseed_np"        "zscore_randomseed_np"
# [21] "zscore_topseed_np_scaled"    "zscore_bottomseed_np_scaled"
# [23] "zscore_randomseed_np_scaled"
tapply(covidb1_trs$Zscore, covidb1_trs$Health_state, mean)
#           D           H
#  0.02196495 -0.08175314
tapply(covidb1_trs$zscore_topseed_np_scaled, covidb1_trs$Health_state, mean)
#        D        H
# 1.803430 1.554548
tapply(covidb1_trs$Zscore, covidb1_trs$Health_state, median)
#          D          H
# -0.2845169 -0.3236276
tapply(covidb1_trs$zscore_topseed_np_scaled, covidb1_trs$Health_state, median)
#        D        H
# 1.452532 1.183256

tapply(covidb1_trs$Zscore, covidb1_trs$Sample_state, mean)
#           H           M           S
# -0.08175314  0.01716032  0.02931954
tapply(covidb1_trs$zscore_topseed_np_scaled, covidb1_trs$Sample_state, mean)
#        H        M        S
# 1.554548 1.787685 1.827532
tapply(covidb1_trs$Zscore, covidb1_trs$Sample_state, median)
#          H          M          S
# -0.3236276 -0.2848795 -0.2839558
tapply(covidb1_trs$zscore_topseed_np_scaled, covidb1_trs$Sample_state, median)
#        H        M        S
# 1.183256 1.426905 1.491395

trait_name="COVID19-B1"
library(ggpubr)
# Visualize: Specify the comparisons you want
my_comparisons <- list( c("H", "M"), c("H", "S"), c("M", "S") )
p <- ggboxplot(covidb1_trs, x = "Sample_state", y = "zscore_topseed_np_scaled",
          fill = "Sample_state", palette = "npg", outlier.shape = NA,
          order = c("H", "M", "S")) + ylim(0, 10) +
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value
ggsave(p, file = "/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/covidb1_you_2021_covidpbmc/TRSboxplot_Sample_state2.pdf", width = 4, height = 4)

covid19hg_B1_df <- covidb1_trs
save(covid19hg_B1_df, file="/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/mono_you_2021_covidpbmc/covid_b1_trs_df.rda")

# ----------------
# ---------------- 11082021
# ----------------
# true cells
# some codes from 20210921-cleandata_largeheme_trait_truecell
library(stringr)
source("/broad/sankaranlab/fyu/script/general_script/general_utils.r")
load("/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata/you_2021_covidpbmc_knngraph.rda")
setwd("/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/mono_you_2021_covidpbmc")
load("covid_b1_trs_df.rda")
get_sigcell_simple <- function(knn_sparse_mat=mutualknn30, 
                        topidx=baso_df$seed_p0.05, 
                        topseed_trs=baso_df$original_seed_p0.05_gamma005_ceiling099, 
                        permutation_times=1000,
                        true_cell_threshold=950, 
                        out_rda="true_cell_df.rda", 
                        mycores=4){
        cell_mat <- data.frame(cell=1:nrow(knn_sparse_mat), degree=colSums(knn_sparse_mat))
        cell_table <- data.frame(table(cell_mat$degree))


        seed_mat_top  <- data.frame(seed=which(topidx), degree=colSums(knn_sparse_mat[, topidx]))
        summary(seed_mat_top $degree)
        seed_table_top <- data.frame(table(seed_mat_top$degree))
        xx_top <- tapply(cell_mat[, 1], cell_mat[, 2], list) 
        xx2_top <- xx_top[names(xx_top) %in% seed_table_top$Var1]
        # permutation_score_top <- data.frame(matrix(nrow=nrow(knn_sparse_mat), ncol=1000))
        permutation_score_top <- parallel::mclapply(1:permutation_times, mc.cores = mycores, function(i){
                sampled_cellid <- xx2_top %>% 
                            mapply(sample, ., seed_table_top$Freq) %>% 
                            unlist %>% 
                            sort
            xx <- randomWalk_sparse(intM=knn_sparse_mat, queryGenes=as.numeric(sampled_cellid), gamma=0.05)
            if (i %% 100 == 0) {print(i)}
            return(xx)
            }
        )
        names(permutation_score_top) <- paste0("permutation_", 1:permutation_times)
        permutation_score_top <- as.data.frame(permutation_score_top)

        permutation_df_top <- data.frame(matrix(nrow=nrow(knn_sparse_mat), ncol=permutation_times))

        permutation_df_top <- apply(permutation_score_top, 2, function(x) { temp <- x<=topseed_trs; return(temp) } )
        message("more than 999: ", (sum(rowSums(permutation_df_top)>=999)*100)/nrow(permutation_df_top))
        message("more than 950: ", (sum(rowSums(permutation_df_top)>=950)*100)/nrow(permutation_df_top))
        message("more than 990: ", (sum(rowSums(permutation_df_top)>=990)*100)/nrow(permutation_df_top))
        true_cell_top_idx <- rowSums(permutation_df_top)>=true_cell_threshold
        message("Fold of true cell over seed: ", sum(true_cell_top_idx)/sum(topidx)) # fold of true cell over seed
        # [1] 2.384638
        message("How many propertion of seed were true cells: ", sum(true_cell_top_idx & topidx)*100/sum(topidx)) # how many propertion of seed were true cells
        # [1] 67.50104
        message("How many propertion true cell over all cells: ", (sum(true_cell_top_idx)*100)/nrow(permutation_df_top)) # how many propertion true cell over all cells
        # [1] 18.74055
        message("How many propertion of seed over all cells: ", (sum(topidx)*100)/nrow(permutation_df_top)) # how many propertion of seed over all cells
        # [1] 7.858867

        true_cell_top_filter_idx <- 
        ture_cell_df <- data.frame(topidx, 
                            topseed_trs, 
                            true_cell_top_idx)
        save(ture_cell_df, permutation_score_top, file=out_rda)
        return(ture_cell_df)
}


covid19hg_B1_df_ture_cell_df <- get_sigcell_simple(knn_sparse_mat=mutualknn30,
                                                    topidx=covid19hg_B1_df$seed_p0.05_top, 
                                                    topseed_trs=covid19hg_B1_df$zscore_topseed_np, 
                                                    permutation_times=1000,
                                                    true_cell_threshold=950, 
                                                    out_rda="/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/mono_you_2021_covidpbmc/covid-b1_true_cell_df.rda", 
                                                    mycores=4
)
# more than 999: 2.05225900028143
# more than 950: 10.8728595241703
# more than 990: 3.97787543567207
# Fold of true cell over seed: 2.17471314137259
# How many propertion of seed were true cells: 86.3390344230353
# How many propertion true cell over all cells: 10.8728595241703
# How many propertion of seed over all cells: 4.99967527547464

##### deal with only top_idx cells and non
setwd("/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/mono_you_2021_covidpbmc")
load("/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/mono_you_2021_covidpbmc/covid_b1_trs_df.rda")
load("/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/mono_you_2021_covidpbmc/covid-b1_true_cell_df.rda")
baso_df2 <- data.frame(covid19hg_B1_df, true_cell_top_idx=ture_cell_df$true_cell_top_idx)
baso_df2$true_cellplot = "Depleted cells"

baso_df2$true_cellplot[baso_df2$true_cell_top_idx] = "Enriched cells"
baso_df2$true_cellplot <- factor(baso_df2$true_cellplot, levels=c("Enriched cells", "Depleted cells"))
table(baso_df2$true_cellplot)
# Enriched cells Depleted cells
#          10045          82341

# tapply(baso_df2$true_cell_top_idx, baso_df2$cell_cluster, sum)
baso_df3 <- data.frame(baso_df2$cell_cluster, baso_df2$true_cellplot)
colnames(baso_df3) <- c("cell_cluster", "true_cellplot")
library(reshape2)
baso_df3_counts <- melt(table(baso_df3))
names(baso_df3_counts) <- names(baso_df3)
colnames(baso_df3_counts)[ncol(baso_df3_counts)] <- "cell_count"
baso_df3_counts %>% head
baso_df3_counts$cell_num <- table(baso_df2$cell_cluster)[match(as.character(baso_df3_counts$cell_cluster), names(table(baso_df2$cell_cluster)))]
baso_df3_counts$percent <- baso_df3_counts$cell_count/baso_df3_counts$cell_num
#      cell_cluster  true_cellplot cell_count cell_num     percent
# 1  CD14+ Monocytes Enriched cells       4452    11385 0.391040843
# 2  CD16+ Monocytes Enriched cells        785     4716 0.166454623
# 3               DC Enriched cells        370     2127 0.173953926
# 4     Effector CD8 Enriched cells       1746    22423 0.077866476
# 5             MAIT Enriched cells        434     8431 0.051476693
# 6         Memory B Enriched cells         17     3422 0.004967855
temp_df <- baso_df3_counts[baso_df3_counts$true_cellplot=="Enriched cells", ]
baso_df3_counts$cell_cluster <- factor(as.character(baso_df3_counts$cell_cluster), levels=rev(temp_df$cell_cluster[order(temp_df$percent)]))
mean(baso_df3_counts[2:15,5])
[1] 0.05441591
mycolor=c("#2d9b77", "#6a67aa") # green; purple
pp <- ggplot(data=baso_df3_counts, aes(x=cell_cluster, y=cell_count, fill=true_cellplot)) +
geom_bar(stat="identity", color=NA, position=position_dodge())+ scale_fill_manual(values=mycolor)  +
  theme_minimal()+ pretty_plot(fontsize = 10) + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0("covid19hg_B1_truetopcell_barplot_count_noa.pdf"), width = 7, height = 4, useDingbats=FALSE)

mycolor=c("#f8f9fa", "#343a40") # light gray/ dark gray 
baso_df3_counts$true_cellplot <- factor(baso_df3_counts$true_cellplot, levels=c("Depleted cells", "Enriched cells"))
pp <- ggplot(data=baso_df3_counts, aes(x=cell_cluster, y=percent, fill=true_cellplot)) +
geom_bar(position="stack", stat="identity", color=NA)+scale_fill_manual(values=mycolor)  +
  theme_minimal()+ pretty_plot(fontsize = 10) + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0("covid19hg_B1_truetopcell_barplot_percent_noa.pdf"), width = 7, height = 4, useDingbats=FALSE)

##### look at cd14+ cells
baso_df2_mono14 = baso_df2[baso_df2$cell_cluster == "CD14+ Monocytes", ] # 11385    25
table(baso_df2_mono14[, c("Health_state", "true_cellplot")])
#             true_cellplot
# Health_state Enriched cells Depleted cells
#            D           4031           5750
#            H            421           1183
table(baso_df2_mono14[, c("Health_state", "seed_p0.05_top")])
#             seed_p0.05_top
# Health_state FALSE TRUE
#            D  8965  816
#            H  1491  113

table(baso_df2_mono14[, c("Sample_state", "true_cellplot")])
#             true_cellplot
# Sample_state Enriched cells Depleted cells
#            H            421           1183
#            M           2380           3230
#            S           1651           2520



baso_df2_mono14_samplepercent <- data.frame(sampleid=str_split(rownames(baso_df2_mono14), "#", simplify=T)[, 1], baso_df2_mono14$true_cellplot, baso_df2_mono14$zscore_topseed_np_scaled, baso_df2_mono14$cell_cluster)
baso_df2_mono14_samplepercent_enriched <- baso_df2_mono14_samplepercent[baso_df2_mono14_samplepercent$baso_df2_mono14.true_cellplot=="Enriched cells", ]
baso_df2_mono14_samplepercent_depleted <- baso_df2_mono14_samplepercent[baso_df2_mono14_samplepercent$baso_df2_mono14.true_cellplot=="Depleted cells", ]
table(baso_df2_mono14_samplepercent_enriched$sampleid)*100/nrow(baso_df2_mono14_samplepercent_enriched)
table(baso_df2_mono14_samplepercent_depleted$sampleid)*100/nrow(baso_df2_mono14_samplepercent_depleted)

library(ggpubr)
# Visualize: Specify the comparisons you want

p <- ggboxplot(baso_df2_mono14_samplepercent, x = "sampleid", y = "baso_df2_mono14.zscore_topseed_np_scaled",
          fill = "sampleid", palette = "npg", outlier.shape = NA) + ylim(0, 10)
#   stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#   stat_compare_means(label.y = 50)     # Add global p-value
ggsave(p, file = "/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/covidb1_you_2021_covidpbmc/TRSboxplot_mono14_Sampleid2.pdf", width = 9, height = 4)





421/1183
# [1] 0.3558749
2380/3230
# [1] 0.7368421
1651/2520
# [1] 0.6551587

table(baso_df2_mono14[, c("Sample_state", "true_cellplot")]) %>% chisq.test
# 	Pearson's Chi-squared test
# data:  .
# X-squared = 137.72, df = 2, p-value < 2.2e-16

table(baso_df2_mono14[, c("Health_state", "true_cellplot")]) %>% chisq.test
# 	Pearson's Chi-squared test with Yates' continuity correction
# data:  .
# X-squared = 128.98, df = 1, p-value < 2.2e-16

table(baso_df2_mono14[, c("Health_state", "seed_p0.05_top")]) %>% chisq.test
# 	Pearson's Chi-squared test with Yates' continuity correction
# data:  .
# X-squared = 2.9264, df = 1, p-value = 0.08714

odds.ratio(table(baso_df2_mono14[, c("Sample_state", "true_cellplot")]), conf.level = 0.95)

pdf("mosicplot_mono14_enrichornot_samplestate.pdf")
mosaicplot(table(baso_df2_mono14[, c("Sample_state", "true_cellplot")]), col=c("firebrick", "goldenrod1"), cex.axis = 1.2, sub = "Condition", dir = c("h","v"), ylab = "Relative frequency")
dev.off()
pdf("mosicplot_mono14_enrichornot_healthstate.pdf")
mosaicplot(table(baso_df2_mono14[, c("Health_state", "true_cellplot")]), col=c("firebrick", "goldenrod1"), cex.axis = 1.2, sub = "Condition", dir = c("h","v"), ylab = "Relative frequency")
dev.off()
pdf("mosicplot_mono14_seedornot_healthstate.pdf")
mosaicplot(table(baso_df2_mono14[, c("Health_state", "seed_p0.05_top")]), col=c("firebrick", "goldenrod1"), cex.axis = 1.2, sub = "Condition", dir = c("h","v"), ylab = "Relative frequency")
dev.off()

library(epitools)
health_state_tab <- table(baso_df2_mono14[, c("Health_state", "true_cellplot")])
oddsratio(health_state_tab, method = "wald")$measure[-1,]
# estimate    lower    upper
# 1.969916 1.750146 2.217282
sample_state_tab <- table(baso_df2_mono14[, c("Sample_state", "true_cellplot")])
oddsratio(sample_state_tab[c(2, 1),], method = "wald")$measure[-1,]
# estimate    lower    upper
# 2.070509 1.830532 2.341945
oddsratio(sample_state_tab[c(3, 1),], method = "wald")$measure[-1,]
# estimate    lower    upper
# 1.840980 1.620815 2.091053





baso_df2_mono14$true_cellplot = "Depleted cells"
baso_df2_mono14$true_cellplot[baso_df2_mono14$true_cell_top_idx] = "Enriched cells"
baso_df2_mono14$true_cellplot <- factor(baso_df2_mono14$true_cellplot, levels=c("Enriched cells", "Depleted cells"))
trait_name="covid19hg_B1"
mycolor=c("#343a40", "#dee2e6")
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=true_cellplot)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_manual(values=mycolor)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "_scatterplot_true_cell-mono14.pdf"), width = 6.5, height = 5, useDingbats=FALSE)





# ------------------------------------------------------------------------
# integrate with archr obj
# ------------------------------------------------------------------------
ish -l h_vmem=64G
# Granja_2021_largeheme
use Anaconda3
source activate py3
reuse R-4.0
# ~/R/x86_64-pc-linux-gnu-library/4.0
reuse Python-3.6
# pip install umap-learn --user
use .gsl-2.6
R

source("/broad/sankaranlab/fyu/script/general_script/general_utils.r")
library(ArchR)
addArchRThreads(threads = 1) 
addArchRGenome("hg19")
# BiocManager::install("JASPAR2018")
proj <- loadArchRProject(path = "/broad/sankaranlab/fyu/varSC/data/you_2021_covidpbmc/COVID19-fragment/covidpbmc_sample11")

proj <- addMotifAnnotations(proj, motifSet="JASPAR2018", version = 1, name = "Motif_JASPAR2018")
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(proj, force = TRUE, threads = 1)
saveArchRProject(ArchRProj = proj, outputDirectory = "/broad/sankaranlab/fyu/varSC/data/you_2021_covidpbmc/COVID19-fragment/covidpbmc_sample11", load = FALSE)
getAvailableMatrices(ArchRProj = proj)
proj_MotifMatrix <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "Motif_JASPAR2018Matrix",
)
# save(proj_MotifMatrix, file="/broad/sankaranlab/fyu/varSC/data/you_2021_covidpbmc/COVID19-fragment/covidpbmc_sample11/SE-scATAC-Motifs_jaspar2018.rda")
motif_data <- proj_MotifMatrix@assays@data$z

identical(rownames(baso_df2), colnames(motif_data))
# [1] FALSE
sample_idx <- match(rownames(baso_df2), colnames(motif_data)) 
motif_data2 <- motif_data[, sample_idx]
identical(rownames(baso_df2), colnames(motif_data2))
# [1] TRUE
save(motif_data2, file="/broad/sankaranlab/fyu/varSC/data/you_2021_covidpbmc/COVID19-fragment/covidpbmc_sample11/mat-scATAC-Motifs_jaspar2018_sorted.rda")



# # motif_data2@assays@data$z

# # rownames(proj@cellColData)

# sample_idx <- match(rownames(proj@cellColData), colnames(motif_data)) 
# motif_data2 <- motif_data[, sample_idx]
# head(colnames(motif_data2))
# # save(motif_data2, file="/broad/sankaranlab/fyu/general_data/largeHemeArrows/SE-scATAC-Large-Heme-Motifs2.rds")



# ------------------------------
# ------------------------------
idx_enrich <- colnames(motif_data2) %in% rownames(baso_df2_mono14)[baso_df2_mono14$true_cellplot == "Enriched cells"]
idx_deplete <- colnames(motif_data2) %in% rownames(baso_df2_mono14)[baso_df2_mono14$true_cellplot == "Depleted cells"]
# identical(colnames(motif_data2)[colnames(motif_data2) %in% rownames(baso_df2_mono14)], rownames(baso_df2_mono14))

motif_enrich_z <- motif_data2[, idx_enrich] # 4452
motif_deplete_z <- motif_data2[, idx_deplete] # 6933

xx <- t.test(motif_enrich_z[1, ] %>% unlist, motif_deplete_z[1, ] %>% unlist)

mycores=2
res_t.test_list <- parallel::mclapply(1:nrow(motif_enrich_z), mc.cores = mycores, function(i){
        xx <- t.test(motif_enrich_z[i, ] %>% unlist, motif_deplete_z[i, ] %>% unlist)
        res_t.test = c(xx$statistic, xx$p.value)
        if (i %% 100 == 0) {print(i)}
        return(res_t.test)
    }
)
res_t.test_df <- as.data.frame(res_t.test_list) %>% t %>% as.data.frame

rownames(res_t.test_df) <- rownames(motif_enrich_z)

setwd("/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/mono_you_2021_covidpbmc")
res_t.test_df$fdr <- p.adjust(res_t.test_df[, 2], method = "bonferroni")
res_t.test_df$log10fdr <- -log10(res_t.test_df$fdr)
res_t.test_df$log10fdr[res_t.test_df[, 1]<0] <- -res_t.test_df$log10fdr[res_t.test_df[, 1]<0]
res_t.test_df$posorneg <- res_t.test_df[, 1]>0
res_t.test_df <- res_t.test_df[order(res_t.test_df$log10fdr), ]
res_t.test_df$order <- 1:nrow(res_t.test_df)
trait_name="covid19hg_B1"
save(res_t.test_df, )
eithertop10 <- rbind.data.frame(head(res_t.test_df, 10), tail(res_t.test_df, 10))
eithertop10$tf_name <- rownames(eithertop10)

eithertop50 <- rbind.data.frame(head(res_t.test_df, 50), tail(res_t.test_df, 50))
save(res_t.test_df, baso_df2_mono14, file="diff_compare_mono14.rda")
[1] "EOMES_274"       "TBR1_276"        "TBX2_151"        "TBX21_153"
  [5] "TBX5_281"        "RUNX3_144"       "TBX4_280"        "RUNX2_273"
  [9] "TCF7L2_58"       "TBX20_152"       "TBX1_279"        "LEF1_237"
 [13] "TCF3_306"        "TCF4_307"        "ZEB1_412"        "MGA_275"
 [17] "TBX15_277"       "FOXI1_100"       "FOXL1_101"       "FOXD2_323"
 [21] "FOXO4_324"       "FOXO6_325"       "FOXP3_326"       "FOXK2_375"
 [25] "FOXK1_374"       "FOXP1_376"       "SPDEF_146"       "CTCFL_373"
 [29] "ID4_300"         "SNAI2_209"       "ERG_227"         "GATA6_378"
 [33] "FOXO3_102"       "ETS1_228"        "SCRT1_207"       "FLI1_99"
 [37] "ERF_226"         "FOXD1_3"         "FOXG1_79"        "FOXP2_72"
 [41] "FIGLA_296"       "ETV2_230"        "SRY_17"          "GATA5_235"
 [45] "CTCF_24"         "FOXF2_2"         "ELK4_62"         "TCF7L1_452"
 [49] "ZBTB7C_159"      "ZBTB7A_411"      "FOS..JUND_434"   "HIC2_204"
 [53] "FOSL1_40"        "RXRB_327"        "ESR1_189"        "PLAG1_32"
 [57] "RARA..RXRA_30"   "TEAD4_283"       "TEAD1_154"       "NR2C2_52"
 [61] "RARA..RXRG_442"  "NR1H4_393"       "MEF2C_49"        "TEF_318"
 [65] "TEAD2_407"       "MEF2A_241"       "TEAD3_282"       "BACH2_371"
 [69] "NFIC..TLX1_22"   "ESRRB_190"       "PPARA..RXRA_441" "NR2F2_394"
 [73] "MEF2D_242"       "RARA.var.2_195"  "ESR2_63"         "NR4A2_31"
 [77] "JUN_44"          "MAF..NFE2_50"    "MEF2B_122"       "HNF4G_43"
 [81] "RARA_194"        "MAFG_121"        "MAFK_384"        "MAFF_383"
 [85] "NRL_317"         "HLF_108"         "ELF5_94"         "MZF1_5"
 [89] "NR2F1_191"       "ZNF263_60"       "EWSR1.FLI1_25"   "ATF4_308"
 [93] "CEBPE_313"       "CEBPB_311"       "CEBPD_312"       "CEBPG_314"
 [97] "CEBPA_61"        "SPIB_16"         "SPI1_147"        "SPIC_148"

marker tf:
ELK4_62, ltf
CEBPD_312, CEBPE_313, JUN_44, FOS..JUND_434, SPI1_147

marker genes:

CD83

S100A8, S100A12, SELL, PLAC8

jaspar 
https://jaspar.genereg.net/matrix/MA0125.1/
https://jaspar.genereg.net/matrix/MA0043.1/

manual_selected_top <- rbind.data.frame(
                                rbind.data.frame(head(res_t.test_df, 50), tail(res_t.test_df, 50))[c(1, 4, 6, 9, 12, 14, 15, 50,    51, 53, 77, 78, 86, 94, 95, 98, 99), ], 
                                res_t.test_df[c(300, 322), ]
                                )
manual_selected_top$tf_name <- str_split(rownames(manual_selected_top), "_", simplify=T)[, 1]

grep("IRF", rownames(res_t.test_df))
# [1]  70  79  85  87 104 118 300 322
grep("IRF", rownames(res_t.test_df), value=T)
# [1] "IRF5_451" "IRF7_240" "IRF9_115" "IRF4_450" "IRF8_114" "IRF2_4"   "IRF3_449"
# [8] "IRF1_65"

# c("#f37010", "#3893c8") #yellow and blue
p <- ggplot(res_t.test_df, aes(x = order, y = log10fdr, fill = posorneg)) +
            geom_col(position = "identity") +
            scale_fill_manual(values = c("#f37010", "#3893c8"), guide = FALSE) + 
            ggrepel::geom_text_repel(data=eithertop10, aes(order, log10fdr, label = tf_name), color = 'black',
                        size = 2, max.overlaps=300) + xlab("TF motif variable enrichment between trait- relevant and non-relevant cells") +ylab("-log10(Adjusted P value)")+ pretty_plot(fontsize = 10)
ggsave(p, file = paste0(trait_name, "_heterogeneity_motif-mono14.pdf"), width = 7.5, height = 4, useDingbats=FALSE)

manual_selected_top <- manual_selected_top[-c(18, 19), ]
p <- ggplot(res_t.test_df, aes(x = order, y = log10fdr, fill = posorneg)) +
            geom_col(position = "identity") +
            scale_fill_manual(values = c("#647eba", "#f37010"), guide = FALSE) + 
            ggrepel::geom_label_repel(data=manual_selected_top, aes(order, log10fdr, label = tf_name), color = 'black', fill="white",
                        size = 1.5, 
                        point.padding = 0, # additional padding around each point
    min.segment.length = 0, # draw all line segments
    max.time = 1, max.iter = 1e5, # stop after 1 second, or after 100,000 iterations
    box.padding = 0.3 # additional padding around each text label
    ) + xlab("TF motif variable enrichment between trait- relevant and non-relevant cells") +ylab("-log10(Adjusted P value)")+ pretty_plot(fontsize = 10)
ggsave(p, file = paste0(trait_name, "_heterogeneity_motif-mono14_manual.pdf"), width = 6, height = 3.5, useDingbats=FALSE)


#  [1] "EOMES_274"     "TBX21_153"     "RUNX3_144"     "TCF7L2_58"
#  [5] "LEF1_237"      "TCF4_307"      "ZEB1_412"      "ZBTB7A_411"
#  [9] "FOS..JUND_434" "FOSL1_40"      "JUN_44"        "MAF..NFE2_50"
# [13] "HLF_108"       "CEBPB_311"     "CEBPD_312"     "SPIB_16"
# [17] "SPI1_147"

# --------------------------------------------------------
trait_name="covid-b1"
tf_name="CEBPB_311"
tf_idx <- which(rownames(motif_enrich_z)==tf_name)
motif_enrich_z_tf <- data.frame("Enriched cells", unlist(motif_enrich_z[tf_idx, ]))
motif_deplete_z_tf <- data.frame("Depleted cells", unlist(motif_deplete_z[tf_idx, ]))
colnames(motif_enrich_z_tf) <- colnames(motif_deplete_z_tf) <- c("cell_group", "score")
motif_box_plot <- rbind.data.frame(motif_enrich_z_tf, motif_deplete_z_tf)
pp <- ggplot(data=motif_box_plot,  aes(x=cell_group, y=score))  +
        geom_violin(adjust=1, aes(fill=cell_group, color=cell_group)) + scale_color_manual(values=c("#3893c8", "#f37010")) +scale_fill_manual(values=c("#3893c8", "#f37010"))+
        geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
        ylab("Motif enrichment") + xlab("") + ggtitle(tf_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "-", tf_name,"-violinplot_heterogeneity_motif-mono14.pdf"), width = 2, height = 3, useDingbats=FALSE)
# scatterplot
tf_name="CEBPB_311"
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$CEBPB_311 <- motif_box_plot$score[match_idx]
baso_df2_mono14$CEBPB_311[baso_df2_mono14$CEBPB_311>quantile(baso_df2_mono14$CEBPB_311, 0.99)] <- quantile(baso_df2_mono14$CEBPB_311, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=CEBPB_311)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=viridis)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "-", tf_name,"_scatterplot_tfenrichment-mono14.pdf"), width = 6.5, height = 5, useDingbats=FALSE)
# --------------------------------------------------------
trait_name="covid-b1"
tf_name="SPI1_147"
tf_idx <- which(rownames(motif_enrich_z)==tf_name)
motif_enrich_z_tf <- data.frame("Enriched cells", unlist(motif_enrich_z[tf_idx, ]))
motif_deplete_z_tf <- data.frame("Depleted cells", unlist(motif_deplete_z[tf_idx, ]))
colnames(motif_enrich_z_tf) <- colnames(motif_deplete_z_tf) <- c("cell_group", "score")
motif_box_plot <- rbind.data.frame(motif_enrich_z_tf, motif_deplete_z_tf)
pp <- ggplot(data=motif_box_plot,  aes(x=cell_group, y=score))  +
        geom_violin(adjust=1, aes(fill=cell_group, color=cell_group)) + scale_color_manual(values=c("#3893c8", "#f37010")) +scale_fill_manual(values=c("#3893c8", "#f37010"))+
        geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
        ylab("Motif enrichment") + xlab("") + ggtitle(tf_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "-", tf_name,"-violinplot_heterogeneity_motif-mono14.pdf"), width = 2, height = 3, useDingbats=FALSE)
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$SPI1_147 <- motif_box_plot$score[match_idx]
baso_df2_mono14$SPI1_147[baso_df2_mono14$SPI1_147>quantile(baso_df2_mono14$SPI1_147, 0.99)] <- quantile(baso_df2_mono14$SPI1_147, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=SPI1_147)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=viridis)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "-", tf_name,"_scatterplot_tfenrichment-mono14.pdf"), width = 6.5, height = 5, useDingbats=FALSE)
# --------------------------------------------------------
trait_name="covid-b1"
tf_name="SPIB_16"
tf_idx <- which(rownames(motif_enrich_z)==tf_name)
motif_enrich_z_tf <- data.frame("Enriched cells", unlist(motif_enrich_z[tf_idx, ]))
motif_deplete_z_tf <- data.frame("Depleted cells", unlist(motif_deplete_z[tf_idx, ]))
colnames(motif_enrich_z_tf) <- colnames(motif_deplete_z_tf) <- c("cell_group", "score")
motif_box_plot <- rbind.data.frame(motif_enrich_z_tf, motif_deplete_z_tf)
pp <- ggplot(data=motif_box_plot,  aes(x=cell_group, y=score))  +
        geom_violin(adjust=1, aes(fill=cell_group, color=cell_group)) + scale_color_manual(values=c("#3893c8", "#f37010")) +scale_fill_manual(values=c("#3893c8", "#f37010"))+
        geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
        ylab("Motif enrichment") + xlab("") + ggtitle(tf_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "-", tf_name,"-violinplot_heterogeneity_motif-mono14.pdf"), width = 2, height = 3, useDingbats=FALSE)
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$SPIB_16 <- motif_box_plot$score[match_idx]
baso_df2_mono14$SPIB_16[baso_df2_mono14$SPIB_16>quantile(baso_df2_mono14$SPIB_16, 0.99)] <- quantile(baso_df2_mono14$SPIB_16, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=SPIB_16)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=viridis)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "-", tf_name,"_scatterplot_tfenrichment-mono14.pdf"), width = 6.5, height = 5, useDingbats=FALSE)
# --------------------------------------------------------
trait_name="covid-b1"
tf_name="JUN_44"
tf_idx <- which(rownames(motif_enrich_z)==tf_name)
motif_enrich_z_tf <- data.frame("Enriched cells", unlist(motif_enrich_z[tf_idx, ]))
motif_deplete_z_tf <- data.frame("Depleted cells", unlist(motif_deplete_z[tf_idx, ]))
colnames(motif_enrich_z_tf) <- colnames(motif_deplete_z_tf) <- c("cell_group", "score")
motif_box_plot <- rbind.data.frame(motif_enrich_z_tf, motif_deplete_z_tf)
pp <- ggplot(data=motif_box_plot,  aes(x=cell_group, y=score))  +
        geom_violin(adjust=1, aes(fill=cell_group, color=cell_group)) + scale_color_manual(values=c("#3893c8", "#f37010")) +scale_fill_manual(values=c("#3893c8", "#f37010"))+
        geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
        ylab("Motif enrichment") + xlab("") + ggtitle(tf_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "-", tf_name,"-violinplot_heterogeneity_motif-mono14.pdf"), width = 2, height = 3, useDingbats=FALSE)
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$JUN_44 <- motif_box_plot$score[match_idx]
baso_df2_mono14$JUN_44[baso_df2_mono14$JUN_44>quantile(baso_df2_mono14$JUN_44, 0.99)] <- quantile(baso_df2_mono14$JUN_44, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=JUN_44)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=viridis)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "-", tf_name,"_scatterplot_tfenrichment-mono14.pdf"), width = 6.5, height = 5, useDingbats=FALSE)

# --------------------------------------------------------
trait_name="covid-b1"
tf_name="EOMES_274"
tf_idx <- which(rownames(motif_enrich_z)==tf_name)
motif_enrich_z_tf <- data.frame("Enriched cells", unlist(motif_enrich_z[tf_idx, ]))
motif_deplete_z_tf <- data.frame("Depleted cells", unlist(motif_deplete_z[tf_idx, ]))
colnames(motif_enrich_z_tf) <- colnames(motif_deplete_z_tf) <- c("cell_group", "score")
motif_box_plot <- rbind.data.frame(motif_enrich_z_tf, motif_deplete_z_tf)
pp <- ggplot(data=motif_box_plot,  aes(x=cell_group, y=score))  +
        geom_violin(adjust=1, aes(fill=cell_group, color=cell_group)) + scale_color_manual(values=c("#3893c8", "#f37010")) +scale_fill_manual(values=c("#3893c8", "#f37010"))+
        geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
        ylab("Motif enrichment") + xlab("") + ggtitle(tf_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "-", tf_name,"-violinplot_heterogeneity_motif-mono14.pdf"), width = 2, height = 3, useDingbats=FALSE)
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$EOMES_274 <- motif_box_plot$score[match_idx]
baso_df2_mono14$EOMES_274[baso_df2_mono14$EOMES_274>quantile(baso_df2_mono14$EOMES_274, 0.99)] <- quantile(baso_df2_mono14$EOMES_274, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=EOMES_274)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=blueYellow)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "-", tf_name,"_scatterplot_tfenrichment-mono14.pdf"), width = 6.5, height = 5, useDingbats=FALSE)
# --------------------------------------------------------
trait_name="covid-b1"
tf_name="TBX21_153"
tf_idx <- which(rownames(motif_enrich_z)==tf_name)
motif_enrich_z_tf <- data.frame("Enriched cells", unlist(motif_enrich_z[tf_idx, ]))
motif_deplete_z_tf <- data.frame("Depleted cells", unlist(motif_deplete_z[tf_idx, ]))
colnames(motif_enrich_z_tf) <- colnames(motif_deplete_z_tf) <- c("cell_group", "score")
motif_box_plot <- rbind.data.frame(motif_enrich_z_tf, motif_deplete_z_tf)
pp <- ggplot(data=motif_box_plot,  aes(x=cell_group, y=score))  +
        geom_violin(adjust=1, aes(fill=cell_group, color=cell_group)) + scale_color_manual(values=c("#3893c8", "#f37010")) +scale_fill_manual(values=c("#3893c8", "#f37010"))+
        geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
        ylab("Motif enrichment") + xlab("") + ggtitle(tf_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "-", tf_name,"-violinplot_heterogeneity_motif-mono14.pdf"), width = 2, height = 3, useDingbats=FALSE)
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$TBX21_153 <- motif_box_plot$score[match_idx]
baso_df2_mono14$TBX21_153[baso_df2_mono14$TBX21_153>quantile(baso_df2_mono14$TBX21_153, 0.99)] <- quantile(baso_df2_mono14$TBX21_153, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=TBX21_153)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=blueYellow)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "-", tf_name,"_scatterplot_tfenrichment-mono14.pdf"), width = 6.5, height = 5, useDingbats=FALSE)
# --------------------------------------------------------
trait_name="covid-b1"
tf_name="RUNX3_144"
tf_idx <- which(rownames(motif_enrich_z)==tf_name)
motif_enrich_z_tf <- data.frame("Enriched cells", unlist(motif_enrich_z[tf_idx, ]))
motif_deplete_z_tf <- data.frame("Depleted cells", unlist(motif_deplete_z[tf_idx, ]))
colnames(motif_enrich_z_tf) <- colnames(motif_deplete_z_tf) <- c("cell_group", "score")
motif_box_plot <- rbind.data.frame(motif_enrich_z_tf, motif_deplete_z_tf)
pp <- ggplot(data=motif_box_plot,  aes(x=cell_group, y=score))  +
        geom_violin(adjust=1, aes(fill=cell_group, color=cell_group)) + scale_color_manual(values=c("#3893c8", "#f37010")) +scale_fill_manual(values=c("#3893c8", "#f37010"))+
        geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
        ylab("Motif enrichment") + xlab("") + ggtitle(tf_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "-", tf_name,"-violinplot_heterogeneity_motif-mono14.pdf"), width = 2, height = 3, useDingbats=FALSE)
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$RUNX3_144 <- motif_box_plot$score[match_idx]
baso_df2_mono14$RUNX3_144[baso_df2_mono14$RUNX3_144>quantile(baso_df2_mono14$RUNX3_144, 0.99)] <- quantile(baso_df2_mono14$RUNX3_144, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=RUNX3_144)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=blueYellow)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "-", tf_name,"_scatterplot_tfenrichment-mono14.pdf"), width = 6.5, height = 5, useDingbats=FALSE)
# --------------------------------------------------------
trait_name="covid-b1"
tf_name="TCF7L2_58"
tf_idx <- which(rownames(motif_enrich_z)==tf_name)
motif_enrich_z_tf <- data.frame("Enriched cells", unlist(motif_enrich_z[tf_idx, ]))
motif_deplete_z_tf <- data.frame("Depleted cells", unlist(motif_deplete_z[tf_idx, ]))
colnames(motif_enrich_z_tf) <- colnames(motif_deplete_z_tf) <- c("cell_group", "score")
motif_box_plot <- rbind.data.frame(motif_enrich_z_tf, motif_deplete_z_tf)
pp <- ggplot(data=motif_box_plot,  aes(x=cell_group, y=score))  +
        geom_violin(adjust=1, aes(fill=cell_group, color=cell_group)) + scale_color_manual(values=c("#3893c8", "#f37010")) +scale_fill_manual(values=c("#3893c8", "#f37010"))+
        geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
        ylab("Motif enrichment") + xlab("") + ggtitle(tf_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "-", tf_name,"-violinplot_heterogeneity_motif-mono14.pdf"), width = 2, height = 3, useDingbats=FALSE)
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$TCF7L2_58 <- motif_box_plot$score[match_idx]
baso_df2_mono14$TCF7L2_58[baso_df2_mono14$TCF7L2_58>quantile(baso_df2_mono14$TCF7L2_58, 0.99)] <- quantile(baso_df2_mono14$TCF7L2_58, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=TCF7L2_58)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=blueYellow)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "-", tf_name,"_scatterplot_tfenrichment-mono14.pdf"), width = 6.5, height = 5, useDingbats=FALSE)
# --------------------------------------------------------
trait_name="covid-b1"
tf_name="TCF4_307"
tf_idx <- which(rownames(motif_enrich_z)==tf_name)
motif_enrich_z_tf <- data.frame("Enriched cells", unlist(motif_enrich_z[tf_idx, ]))
motif_deplete_z_tf <- data.frame("Depleted cells", unlist(motif_deplete_z[tf_idx, ]))
colnames(motif_enrich_z_tf) <- colnames(motif_deplete_z_tf) <- c("cell_group", "score")
motif_box_plot <- rbind.data.frame(motif_enrich_z_tf, motif_deplete_z_tf)
pp <- ggplot(data=motif_box_plot,  aes(x=cell_group, y=score))  +
        geom_violin(adjust=1, aes(fill=cell_group, color=cell_group)) + scale_color_manual(values=c("#3893c8", "#f37010")) +scale_fill_manual(values=c("#3893c8", "#f37010"))+
        geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
        ylab("Motif enrichment") + xlab("") + ggtitle(tf_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "-", tf_name,"-violinplot_heterogeneity_motif-mono14.pdf"), width = 2, height = 3, useDingbats=FALSE)
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$TCF4_307 <- motif_box_plot$score[match_idx]
baso_df2_mono14$TCF4_307[baso_df2_mono14$TCF4_307>quantile(baso_df2_mono14$TCF4_307, 0.99)] <- quantile(baso_df2_mono14$TCF4_307, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=TCF4_307)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=blueYellow)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "-", tf_name,"_scatterplot_tfenrichment-mono14.pdf"), width = 6.5, height = 5, useDingbats=FALSE)





# --------------------------------------------------------
# boxplot & scatterplot
Jasminium_polyanthum_A <- c("#f7fcfe", "#fef5f9", "#fad5e5", "#f8bdd7", "#d35c95", "#af2963", "#821434", "#45101b", "#230c11")

tf_name="CEBPB_202"
tf_idx <- which(rownames(motif_enrich_z)==tf_name)
motif_enrich_z_tf <- data.frame("Enriched cells", unlist(motif_enrich_z[tf_idx, ]))
motif_deplete_z_tf <- data.frame("Depleted cells", unlist(motif_deplete_z[tf_idx, ]))
colnames(motif_enrich_z_tf) <- colnames(motif_deplete_z_tf) <- c("cell_group", "score")
motif_box_plot <- rbind.data.frame(motif_enrich_z_tf, motif_deplete_z_tf)
pp <- ggplot(data=motif_box_plot,  aes(x=cell_group, y=score))  +
        geom_violin(adjust=1, aes(fill=cell_group, color=cell_group)) + scale_color_manual(values=c("#f37010", "#3893c8")) +scale_fill_manual(values=c("#f37010", "#3893c8"))+
        geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
        ylab("Motif enrichment") + xlab("") + ggtitle(tf_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "-", tf_name,"-violinplot_heterogeneity_motif-Late.Ery.pdf"), width = 2.5, height = 5, useDingbats=FALSE)
tf_name="CEBPB_202"
match_idx <- match(rownames(baso_df2_Late.Ery), rownames(motif_box_plot))
baso_df2_Late.Ery$CEBPB_202 <- motif_box_plot$score[match_idx]
baso_df2_Late.Ery$CEBPB_202[baso_df2_Late.Ery$CEBPB_202>quantile(baso_df2_Late.Ery$CEBPB_202, 0.99)] <- quantile(baso_df2_Late.Ery$CEBPB_202, 0.99)
# Late.Ery
p <- ggplot(data=baso_df2_Late.Ery, aes(UMAP1, UMAP2, color=CEBPB_202)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=Jasminium_polyanthum_A)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "-", tf_name,"_scatterplot_tfenrichment-Late.Ery.pdf"), width = 6.5, height = 5, useDingbats=FALSE)



tf_name="IRF8_964"
tf_idx <- which(rownames(motif_enrich_z)==tf_name)
motif_enrich_z_tf <- data.frame("Enriched cells", unlist(motif_enrich_z[tf_idx, ]))
motif_deplete_z_tf <- data.frame("Depleted cells", unlist(motif_deplete_z[tf_idx, ]))
colnames(motif_enrich_z_tf) <- colnames(motif_deplete_z_tf) <- c("cell_group", "score")
motif_box_plot <- rbind.data.frame(motif_enrich_z_tf, motif_deplete_z_tf)
pp <- ggplot(data=motif_box_plot,  aes(x=cell_group, y=score))  +
        geom_violin(adjust=1, aes(fill=cell_group, color=cell_group)) + scale_color_manual(values=c("#f37010", "#3893c8")) +scale_fill_manual(values=c("#f37010", "#3893c8"))+
        geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
        ylab("Motif enrichment") + xlab("") + ggtitle(tf_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "-", tf_name,"-violinplot_heterogeneity_motif-Late.Ery.pdf"), width = 2.5, height = 5, useDingbats=FALSE)
# --------------------------------------------------------
# scatterplot
tf_name="IRF8_964"
match_idx <- match(rownames(baso_df2_Late.Ery), rownames(motif_box_plot))
baso_df2_Late.Ery$IRF8_964 <- motif_box_plot$score[match_idx]
baso_df2_Late.Ery$IRF8_964[baso_df2_Late.Ery$IRF8_964>quantile(baso_df2_Late.Ery$IRF8_964, 0.99)] <- quantile(baso_df2_Late.Ery$IRF8_964, 0.99)
# Late.Ery
p <- ggplot(data=baso_df2_Late.Ery, aes(UMAP1, UMAP2, color=IRF8_964)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=Jasminium_polyanthum_A)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "-", tf_name,"_scatterplot_tfenrichment-Late.Ery.pdf"), width = 6.5, height = 5, useDingbats=FALSE)



tf_name="GATA2_600"
tf_idx <- which(rownames(motif_enrich_z)==tf_name)
motif_enrich_z_tf <- data.frame("Enriched cells", unlist(motif_enrich_z[tf_idx, ]))
motif_deplete_z_tf <- data.frame("Depleted cells", unlist(motif_deplete_z[tf_idx, ]))
colnames(motif_enrich_z_tf) <- colnames(motif_deplete_z_tf) <- c("cell_group", "score")
motif_box_plot <- rbind.data.frame(motif_enrich_z_tf, motif_deplete_z_tf)
pp <- ggplot(data=motif_box_plot,  aes(x=cell_group, y=score))  +
        geom_violin(adjust=1, aes(fill=cell_group, color=cell_group)) + scale_color_manual(values=c("#f37010", "#3893c8")) +scale_fill_manual(values=c("#f37010", "#3893c8"))+
        geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
        ylab("Motif enrichment") + xlab("") + ggtitle(tf_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "-", tf_name,"-violinplot_heterogeneity_motif-Late.Ery.pdf"), width = 2.5, height = 5, useDingbats=FALSE)
# --------------------------------------------------------
# scatterplot
tf_name="GATA2_600"
match_idx <- match(rownames(baso_df2_Late.Ery), rownames(motif_box_plot))
baso_df2_Late.Ery$GATA2_600 <- motif_box_plot$score[match_idx]
baso_df2_Late.Ery$GATA2_600[baso_df2_Late.Ery$GATA2_600>quantile(baso_df2_Late.Ery$GATA2_600, 0.99)] <- quantile(baso_df2_Late.Ery$GATA2_600, 0.99)
# Late.Ery
p <- ggplot(data=baso_df2_Late.Ery, aes(UMAP1, UMAP2, color=GATA2_600)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors = jdb_palette("brewer_marine"))  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "-", tf_name,"_scatterplot_tfenrichment-Late.Ery.pdf"), width = 6.5, height = 5, useDingbats=FALSE)


tf_name="GATA1_588"
tf_idx <- which(rownames(motif_enrich_z)==tf_name)
motif_enrich_z_tf <- data.frame("Enriched cells", unlist(motif_enrich_z[tf_idx, ]))
motif_deplete_z_tf <- data.frame("Depleted cells", unlist(motif_deplete_z[tf_idx, ]))
colnames(motif_enrich_z_tf) <- colnames(motif_deplete_z_tf) <- c("cell_group", "score")
motif_box_plot <- rbind.data.frame(motif_enrich_z_tf, motif_deplete_z_tf)
pp <- ggplot(data=motif_box_plot,  aes(x=cell_group, y=score))  +
        geom_violin(adjust=1, aes(fill=cell_group, color=cell_group)) + scale_color_manual(values=c("#f37010", "#3893c8")) +scale_fill_manual(values=c("#f37010", "#3893c8"))+
        geom_boxplot(aes(fill=NULL), width=0.075, outlier.shape=NA) +
        ylab("Motif enrichment") + xlab("") + ggtitle(tf_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4))
ggsave(pp, file = paste0(trait_name, "-", tf_name,"-violinplot_heterogeneity_motif-Late.Ery.pdf"), width = 2.5, height = 5, useDingbats=FALSE)
# --------------------------------------------------------
# scatterplot
tf_name="GATA1_588"
match_idx <- match(rownames(baso_df2_Late.Ery), rownames(motif_box_plot))
baso_df2_Late.Ery$GATA1_588 <- motif_box_plot$score[match_idx]
baso_df2_Late.Ery$GATA1_588[baso_df2_Late.Ery$GATA1_588>quantile(baso_df2_Late.Ery$GATA1_588, 0.99)] <- quantile(baso_df2_Late.Ery$GATA1_588, 0.99)
# Late.Ery
p <- ggplot(data=baso_df2_Late.Ery, aes(UMAP1, UMAP2, color=GATA1_588)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors = jdb_palette("brewer_marine")) +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = paste0(trait_name, "-", tf_name,"_scatterplot_tfenrichment-Late.Ery.pdf"), width = 6.5, height = 5, useDingbats=FALSE)









# ---------------------------------------------
# gene subset of Mono
# ---------------------------------------------
# Mono 14
source("/broad/sankaranlab/fyu/script/general_script/general_utils.r")
library(ArchR)
addArchRThreads(threads = 2) 
addArchRGenome("hg19")
# BiocManager::install("JASPAR2018")
proj <- loadArchRProject(path = "/broad/sankaranlab/fyu/varSC/data/you_2021_covidpbmc/COVID19-fragment/covidpbmc_sample11")

idx_Mono <- match(rownames(baso_df2_mono14), getCellNames(proj))
length(idx_Mono)
# [1] 11385
identical(getCellNames(proj)[idx_Mono], rownames(baso_df2_mono14))
# ----
# do not create the folder before run this command
proj_mono14 <- subsetArchRProject(
        ArchRProj = proj,
        cells = getCellNames(proj)[idx_Mono],
        outputDirectory = "/broad/sankaranlab/fyu/varSC/data/you_2021_covidpbmc/COVID19-fragment/covidpbmc_sample11-mono14",
        dropCells = TRUE,
        logFile = NULL,
        threads = getArchRThreads(),
        force = T
        )
identical(getCellNames(proj_mono14), rownames(baso_df2_mono14))
# [1] TRUE
table(proj_mono14@cellColData$BioCluster)
# CD14+ Monocytes
#           11385
proj_mono14@cellColData  <- DataFrame(getCellColData(proj_mono14), baso_df2_mono14[, c(3, 4, 6, 21:25)])
proj_mono14@cellColData$true_cellplot <- as.character(proj_mono14@cellColData$true_cellplot)
markersGS <- getMarkerFeatures(ArchRProj = proj_mono14,
                               useMatrix = "GeneScoreMatrix", 
                               groupBy = "true_cellplot",
                               bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon", 
)

markerList <- getMarkers(markersGS,
                         cutOff = "FDR <= 0.05 & Log2FC >= 0.25",
                         returnGR = FALSE)

markerList@listData$`Enriched cells`
markerList@listData$`Depleted cells`

proj_mono14 <- getGroupBW(proj_mono14, groupBy = "true_cellplot") # export the bigwig files
markersGS <- getMarkerFeatures(ArchRProj = proj_Mono,
                               useMatrix = "GeneScoreMatrix", 
                               groupBy = "enrichornot",
                               bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon", 
)
quantile(assay(markersGS, "Pval")$Enriched, na.rm=T)
#           0%          25%          50%          75%         100%
# 3.539117e-05 2.504780e-01 5.013906e-01 7.567904e-01 1.000000e+00
quantile(assay(markersGS, "Log2FC")$Enriched, na.rm=T)
#            0%           25%           50%           75%          100%
# -3.9433142707 -0.1592069982 -0.0008539071  0.1476454350  5.5166474414
markerList <- getMarkers(markersGS,
                         cutOff = "Pval <= 0.05 & Log2FC >= 0.25",
                         returnGR = FALSE)
markerList$Enriched

proj_mono14 <- getGroupBW(proj_mono14, groupBy = "true_cellplot", tileSize = 20,
  maxCells = 5000,
  ceiling = 2) # export the bigwig files

proj_mono14 <- getGroupBW(proj_mono14, groupBy = "true_cellplot", tileSize = 20,
  maxCells = 5000,
  ceiling = 2) # export the bigwig files


proj_mono14 <- loadArchRProject(path = "/broad/sankaranlab/fyu/varSC/data/you_2021_covidpbmc/COVID19-fragment/covidpbmc_sample11-mono14")
proj_mono14@cellColData  <- DataFrame(getCellColData(proj_mono14), baso_df2_mono14[, c(3, 4, 6, 21:25)])
proj_mono14@cellColData$true_cellplot <- as.character(proj_mono14@cellColData$true_cellplot)
head(proj_mono14@cellColData)
proj_mono14@cellColData$true_cellplot_1k <- "Ambiguous"
proj_mono14@cellColData$true_cellplot_1k[rank(proj_mono14@cellColData$zscore_topseed_np_scaled) <= 1000] <- "Depleted"
proj_mono14@cellColData$true_cellplot_1k[rank(proj_mono14@cellColData$zscore_topseed_np_scaled) >= (length(proj_mono14@cellColData$true_cellplot_1k) -999)] <- "Enriched"
proj_mono14@cellColData$true_cellplot_1k <- as.character(proj_mono14@cellColData$true_cellplot_1k)
proj_mono14 <- getGroupBW(proj_mono14, groupBy = "true_cellplot_1k", tileSize = 50,
  maxCells = 1000,
  ceiling = 4) # export the bigwig files

#----------- chr5 30k
plot_region <- c("chr3:168361070-168363071", "chr3:168382423-168384424", "chr5:163124305-163126306", "chr7:134666121-134668122", "chr7:134703170-134705171", "chr9:23546543-23548544", "chr9:23552800-23554801", "chr13:46864283-46866284", "chr13:46861803-46863804")
plot_region <- data.frame(stringr::str_split(plot_region, ":", simplify=T)[, 1], stringr::str_split(stringr::str_split(plot_region, ":", simplify=T)[, 2], "-", simplify=T))
colnames(plot_region) <- c("chr", "start", "end")
plot_region_g <- plot_region %>% makeGRangesFromDataFrame

start(plot_region_g) <- start(plot_region_g) - 15000
end(plot_region_g) <- end(plot_region_g) + 15000

p <- plotBrowserTrack(proj_mono14, region= plot_region_g[3], groupBy = "true_cellplot_1k", useGroups=c("Depleted", "Enriched"), log2Norm=T, tileSize=20, plotSummary=c("bulkTrack", "scTrack"), scTileSize=1, scCellsMax = 1000)
plotPDF(p, name="test_gb_chr5_true_cellplot_30k.pdf", ArchRProj = proj_mono14)

#----------- chr5 500k
plot_region <- c("chr3:168361070-168363071", "chr3:168382423-168384424", "chr5:163124305-163126306", "chr7:134666121-134668122", "chr7:134703170-134705171", "chr9:23546543-23548544", "chr9:23552800-23554801", "chr13:46864283-46866284", "chr13:46861803-46863804")
plot_region <- data.frame(stringr::str_split(plot_region, ":", simplify=T)[, 1], stringr::str_split(stringr::str_split(plot_region, ":", simplify=T)[, 2], "-", simplify=T))
colnames(plot_region) <- c("chr", "start", "end")
plot_region_g <- plot_region %>% makeGRangesFromDataFrame

start(plot_region_g) <- start(plot_region_g) - 250000
end(plot_region_g) <- end(plot_region_g) + 250000

p <- plotBrowserTrack(proj_mono14, region= plot_region_g[3], groupBy = "true_cellplot_1k", useGroups=c("Depleted", "Enriched"), 
                        log2Norm=T, tileSize=1000, useMatrix="GeneScoreMatrix", plotSummary=c("bulkTrack", "scTrack", "geneTrack", "featureTrack"), scTileSize=1, scCellsMax = 1000)
plotPDF(p, name="test_gb_chr5_true_cellplot_500k.pdf", ArchRProj = proj_mono14)

p <- plotBrowserTrack(proj_mono14, groupBy = "true_cellplot_1k", useGroups=c("Depleted", "Enriched"), 
                        geneSymbol=c("HMMR", "NUDCD2", "CCNG1", "MAT2B"), log2Norm=T, tileSize=1000, useMatrix="GeneScoreMatrix", plotSummary=c("bulkTrack", "scTrack", "geneTrack", "featureTrack"), scTileSize=1, scCellsMax = 1000)
plotPDF(p, name="test_gb_chr5_true_cellplot_500k_gene.pdf", ArchRProj = proj_mono14)

CTGTATT(C)AGATC
CTGTATTAGATC
CTGTATCAGATC
chr5    163125300   163125310

p <- plotBrowserTrack(proj_mono14, groupBy = "true_cellplot_1k", useGroups=c("Depleted", "Enriched"), 
                        geneSymbol=c("HMMR", "NUDCD2", "CCNG1", "MAT2B"), log2Norm=T, tileSize=1000, useMatrix="GeneScoreMatrix", plotSummary=c("bulkTrack", "scTrack", "geneTrack", "featureTrack"), scTileSize=1, scCellsMax = 1000)
plotPDF(p, name="test_gb_chr5_true_cellplot_500k_gene.pdf", ArchRProj = proj_mono14)

#----------- motif presence 
peak_motif_present <- readRDS("/broad/sankaranlab/fyu/varSC/data/you_2021_covidpbmc/COVID19-fragment/covidpbmc_sample11-mono14/Annotations/Motif_JASPAR2018-Matches-In-Peaks.rds")
trait_file <- read.table("/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata2/finemappedtraits_hg19/covid19hg_B1.p1e-6.PP001_filter.txt")
colnames(trait_file)[1:3] <- c("chr", "start", "end")
trait_file_g <- makeGRangesFromDataFrame(trait_file)
ol <- findOverlaps(rowRanges(peak_motif_present), trait_file_g) #note some genes in the peak_obj overlapping
ol
# Hits object with 28 hits and 0 metadata columns:
#        queryHits subjectHits
#        <integer>   <integer>
#    [1]     76281          22
#    [2]     76319         155
#    [3]     76326         175
#    [4]     76326         172
#    [5]     76332         128
#    ...       ...         ...
#   [24]     76380          47
#   [25]     76388         114
#   [26]     76391         126
#   [27]     76392          24
#   [28]    274663           3
#   -------
#   queryLength: 396642 / subjectLength: 265
trait_file[sort(subjectHits(ol)), ]
queryHits(ol) %>% length
# [1] 28
queryHits(ol) %>% unique %>% length
# [1] 23
subjectHits(ol) %>% length
# [1] 28
subjectHits(ol) %>% unique %>% length
# [1] 28
peak_motif_present_overlap <- assay(peak_motif_present)[queryHits(ol) %>% unique, ]
colSums(peak_motif_present_overlap) %>% sort %>% tail(., 100)
    #   ZNF713_364              WT1_378              SP1_379
    #            5                    5                    5
    #    BRCA1_461             SPI1_466             FLI1_493
    #            5                    5                    5
    #     ETV4_505           POU2F1_919           POU3F3_938
    #            5                    5                    5
    #     IRF3_957             IRF4_961           PPARG_1053
    #            5                    5                    5
    #    RXRA_1105           NFKB1_1147             CIC_1207
    #            5                    5                    5
    #    SOX4_1217           STAT3_1280           TEAD1_1319
    #            5                    5                    5
    #    HLTF_1328           CPEB1_1344          BCL11B_1365
    #            5                    5                    5
    #   FOXN3_1372           FOXD4_1388         FOXD4L1_1417
    #            5                    5                    5
    # FOXD4L3_1447         FOXD4L5_1476         FOXD4L6_1505
    #            5                    5                    5
    #    MYF5_1555            KLF2_1566           FOXJ1_1591
    #            5                    5                    5
    #   FOXN2_1614           FOXL2_1625           FOXL2_1627
    #            5                    5                    5
    # FOXD4L4_1658         FOXD4L2_1687           HOXA4_1733
    #            5                    5                    5
    #      TCF3_40            PTF1A_104            BACH1_183
    #            6                    6                    6
    #    BACH1_185           ZBTB7B_311            FOXC2_563
    #            6                    6                    6
    #     IRF5_960           HNF4A_1022           NR4A2_1074
    #            6                    6                    6
    #   TBPL2_1315            BPTF_1341           FOXD4_1406
    #            6                    6                    6
    # FOXD4L1_1436         FOXD4L3_1465         FOXD4L5_1494
    #            6                    6                    6
    # FOXD4L6_1523           FOXJ1_1609           FOXL2_1647
    #            6                    6                    6
    # FOXD4L4_1676         FOXD4L2_1705            SOX6_1760
    #            6                    6                    6
    #   ZNF263_234             CTCF_265              SP3_352
    #            7                    7                    7
    #     SPI1_465              EHF_484             SPIB_489
    #            7                    7                    7
    #      FEV_502            FOXO1_540            FOXJ3_581
    #            7                    7                    7
    #      HDX_781             IRF1_954             IRF8_965
    #            7                    7                    7
    #   RUNX2_1178           RUNX1_1182 ENSG00000250096_1183
    #            7                    7                    7
    #    NFIC_1192           SOX10_1214           FOXN3_1379
    #            7                    7                    7
    #   FOXD4_1383           FOXD4_1409         FOXD4L1_1412
    #            7                    7                    7
    # FOXD4L1_1439         FOXD4L3_1443         FOXD4L3_1468
    #            7                    7                    7
    # FOXD4L5_1472         FOXD4L5_1497         FOXD4L6_1501
    #            7                    7                    7
    # FOXD4L6_1526           FOXJ1_1586           FOXJ1_1612
    #            7                    7                    7
    #   FOXN2_1621           FOXL2_1650         FOXD4L4_1654
    #            7                    7                    7
    # FOXD4L4_1679         FOXD4L2_1683         FOXD4L2_1708
    #            7                    7                    7
    #   ZNF263_235           ZNF148_317            FOXP1_524
    #            8                    8                    8
    #   STAT2_1281             IRF4_963             SRY_1257
    #            8                   10                   10
    #      MAZ_266
    #           11


odd_factor <- nrow(peak_motif_present)/(queryHits(ol) %>% unique %>% length)
oddratio_overlapmotif <- odd_factor * (colSums(peak_motif_present_overlap)/colSums(assay(peak_motif_present)))

tail(sort(oddratio_overlapmotif), 20)
 HOMEZ_853   CGBP_423  CXXC1_428  BRCA1_461  CEBPG_179   IRF5_959  BATF3_167
  3.951720   3.958977   3.958977   3.964802   4.006343   4.022698   4.035405
  OSR1_304   PAX6_889 NR4A2_1074 PROX1_1138   LHX8_761  GATA1_589  RXRG_1059
  4.105047   4.317082   4.387374   4.599361   4.776211   4.778416   4.973651
  OTX1_675  NAIF1_974   LHX6_654   TET1_427  HSFY2_948  NPAS4_112
  5.093120   6.239256   6.877032   7.135004   9.271669   9.596719

# ---------------------------------------------
# gene activity score of cd14 and cd16
# ---------------------------------------------

proj_GeneScoreMatrix <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix"
)

proj_GeneScoreMatrix

proj <- addGeneScoreMatrix(proj)








plot_region <- c("chr3:168361070-168363071", "chr3:168382423-168384424", "chr5:163124305-163126306", "chr7:134666121-134668122", "chr7:134703170-134705171", "chr9:23546543-23548544", "chr9:23552800-23554801", "chr13:46864283-46866284", "chr13:46861803-46863804")
plot_region <- data.frame(stringr::str_split(plot_region, ":", simplify=T)[, 1], stringr::str_split(stringr::str_split(plot_region, ":", simplify=T)[, 2], "-", simplify=T))
colnames(plot_region) <- c("chr", "start", "end")
plot_region_g <- plot_region %>% makeGRangesFromDataFrame

start(plot_region_g) <- start(plot_region_g) - 50000
end(plot_region_g) <- end(plot_region_g) + 50000

p <- plotBrowserTrack(proj_mono14, region= plot_region_g[5], groupBy = "true_cellplot", log2Norm=F, tileSize=50, plotSummary=c("bulkTrack", "scTrack"), scTileSize=4, scCellsMax = 5000)
plotPDF(p, name="test_gb_highppregions.pdf", ArchRProj = proj_mono14)


p <- plotBrowserTrack(proj_mono14, groupBy = "true_cellplot", geneSymbol = "CD247", upstream = 100000, downstream = 100000)
plotPDF(p, name="test_gb_cd247.pdf", ArchRProj = proj_mono14)
p <- plotBrowserTrack(proj_mono14, groupBy = "true_cellplot", geneSymbol = "LRRC63", upstream = 500000, downstream = 500000)
plotPDF(p, name="test_gb_LRRC63.pdf", ArchRProj = proj_mono14)
p <- plotBrowserTrack(proj_mono14, groupBy = "true_cellplot", geneSymbol = "LINC01198", upstream = 500000, downstream = 500000)
plotPDF(p, name="test_gb_LINC01198.pdf", ArchRProj = proj_mono14)
p <- plotBrowserTrack(proj_mono14, groupBy = "true_cellplot", geneSymbol = "EGFEM1P", upstream = 500000, downstream = 500000)
plotPDF(p, name="test_gb_EGFEM1P.pdf", ArchRProj = proj_mono14)
p <- plotBrowserTrack(proj_mono14, groupBy = "true_cellplot", geneSymbol = "AGBL3", upstream = 500000, downstream = 500000)
plotPDF(p, name="test_gb_AGBL3.pdf", ArchRProj = proj_mono14)
p <- plotBrowserTrack(proj_mono14, groupBy = "true_cellplot", geneSymbol = "CALD1", upstream = 500000, downstream = 500000)
plotPDF(p, name="test_gb_CALD1.pdf", ArchRProj = proj_mono14)
p <- plotBrowserTrack(proj_mono14, groupBy = "true_cellplot", geneSymbol = "NR_121602.1", upstream = 500000, downstream = 500000)
plotPDF(p, name="test_gb_NR_121602.1.pdf", ArchRProj = proj_mono14)
p <- plotBrowserTrack(proj_mono14, groupBy = "true_cellplot", geneSymbol = "CCNG1", upstream = 500000, downstream = 500000)
plotPDF(p, name="test_gb_CCNG1.pdf", ArchRProj = proj_mono14)
p <- plotBrowserTrack(proj_mono14, groupBy = "true_cellplot", geneSymbol = "MAT2B", upstream = 500000, downstream = 500000)
plotPDF(p, name="test_gb_MAT2B.pdf", ArchRProj = proj_mono14)

name = "Plot",
  width = 6,
  height = 6,
  ArchRProj = NULL
# ------------------------------------------------------
# you_2021_covidpbmc (hg19)
# mono
source("/broad/sankaranlab/fyu/script/general_script/general_utils.r")
# GC-content corrected enrichment score (gchromVAR)
load("/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata/you_2021_covidpbmc_SE_gvar.rda") # 
load("/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata/you_2021_covidpbmc_cellmetainfo.rda") # 

trait_files="/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata2/finemappedtraits_hg19/mono.PP001.bed"
bcx <- importBedScore(rowRanges(SE_gvar), trait_files, colidx=5)
# head(bcx@assays@data@listData$weights[rowSums(bcx@assays@data@listData$weights) !=0,], 10)
ukbb_wDEV <- computeWeightedDeviations(SE_gvar, bcx, background_peaks = SE_gvar_bg)

# Reformat results
cell_name <- rownames(t(assays(ukbb_wDEV)[["z"]]))
cell_cluster <- meta_df$cell_cluster
UMAP1 <- meta_df$x
UMAP2 <- meta_df$y
# temp_df2 <- data.frame(cell_name, cell_cluster, UMAP1, UMAP2, t(assays(ukbb_wDEV)[["z"]]))
# zscoreWeighted2 <- melt(temp_df2, c(1, 2, 3, 4))
zscoreWeighted2 <- data.frame(cell_name, cell_cluster, UMAP1, UMAP2, "mono", t(assays(ukbb_wDEV)[["z"]]), meta_df$Health_state, meta_df$Sample_state)

tapply(zscoreWeighted2$mono.PP001, zscoreWeighted2$cell_cluster, mean)


colnames(zscoreWeighted2) <- c("cell_name", "cell_cluster", "UMAP1", "UMAP2", "trait", "Zscore", "Health_state", "Sample_state") # 92386 6
table(zscoreWeighted2$trait)
# mono
#      92386

save(ukbb_wDEV, zscoreWeighted2, meta_df, file="/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata3/mono_you_2021_covidpbmc_gscore.rda")
# load("/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata3/covidb1_you_2021_covidpbmc_gscore.rda")




load("/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata/you_2021_covidpbmc_knngraph.rda")
zscoreWeighted2$cell_cluster2 <- zscoreWeighted2$cell_cluster
# celltypeanno <- c("Isocoritical excitatory", "Striatal inhibitory (major)", "Hippocampal excitatory 1", "Hippocampal excitatory 2", "Nigral neurons (unclassified)", "Nigral cells (unclassified)", "Neurons (unclassified)", "OPCs1", "OPCs2", "Nigral OPCs", "Isocortical inhibitory", "Striatal inhibitory (minor)", "Astrocytes (unclassified)", "Nigral astrocytes", "Isocortical astrocytes", "Striatal astrocytes", "Astrocytes2 (unclassified)", "Potential doublets", "Oligodendrocytes1", "Oligodendrocytes2", "Oligodendrocytes3", "Oligodendrocytes4", "Oligodendrocytes5", "Microglia")
# zscoreWeighted2$cell_cluster2 <- factor(zscoreWeighted2$cell_cluster2, levels=celltypeanno)

mono_df <- zscoreWeighted2
mono_trs <- trait_netpropagation(trait_df=mono_df, seedcutoff=0.005, mutual_knn_graph_mat=mutualknn30, trait_name="mono", workdir="/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/mono_you_2021_covidpbmc/covidb1", umap_pointsize=0.05, umap_pointalpha= 0.55, umap_color="Photinia_fraseri_B")



p <- ggplot(data=covidb1_df, aes(UMAP1, UMAP2, color=cell_cluster2)) + geom_point(size=0.2, na.rm = TRUE, alpha = 0.4) + 
    scale_color_discrete() +
     pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = "/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/covidb1_you_2021_covidpbmc/covidb1/test.png", width = 15, height = 5)


# > tapply(zscoreWeighted2$Zscore, zscoreWeighted2$Health_state, mean)
#            D            H
#  0.006445632 -0.093942972
# > tapply(zscoreWeighted2$Zscore, zscoreWeighted2$Sample_state, mean)
#            H            M            S
# -0.093942972 -0.006353707  0.026037964

zscoreWeighted3 <-  zscoreWeighted2[order(zscoreWeighted2$Zscore), ]
# > table(tail(zscoreWeighted3, 500)$cell_cluster)

# CD14+ Monocytes CD16+ Monocytes              DC    Effector CD8            MAIT
#             342              38              32              13              36
#        Memory B      Memory CD4              NK         Naive B         Naive T
#               5               1               1              10               2
#         Other B         Other T          Plasma       Undefined
#               5               2               5               8
# > table(tail(zscoreWeighted3, 1000)$cell_cluster)

# CD14+ Monocytes CD16+ Monocytes              DC    Effector CD8            MAIT
#             645              83              50              44              70
#        Memory B      Memory CD4              NK         Naive B         Naive T
#              14               9               2              21              11
#         Other B         Other T          Plasma       Undefined
#              10               7               9              25
# > table(tail(zscoreWeighted3, 300)$cell_cluster)

# CD14+ Monocytes CD16+ Monocytes              DC    Effector CD8            MAIT
#             206              21              22               8              23
#        Memory B      Memory CD4              NK         Naive B         Naive T
#               2               1               1               4               1
#         Other B         Other T          Plasma       Undefined
#               3               2               1               5

zscoreWeighted3$gchromVAR_pvalue <- pnorm(zscoreWeighted3$Zscore, lower.tail = FALSE)
    # head(trait_df$gchromVAR_pvalue)
sum(zscoreWeighted3$gchromVAR_pvalue<0.05)
[1] 7946

zscoreWeighted3$seedornot = 0
zscoreWeighted3$seedornot[(length(zscoreWeighted3$seedornot)-461 ):length(zscoreWeighted3$seedornot)] = 1
zscoreWeighted3$seedornot2 <- sample(zscoreWeighted3$seedornot)
setwd("/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/mono_you_2021_covidpbmc")
p <- ggplot(data=zscoreWeighted3, aes(UMAP1, UMAP2, color=seedornot)) + geom_point(size=0.05) +
    pretty_plot() + xlab("UMAP 1") + ylab("UMAP 2") + theme_void()  + theme(legend.position = "none")
ggsave(p, file = "test_seed461.png", width = 5, height = 5)
zscoreWeighted3 <- zscoreWeighted3[order(zscoreWeighted3$seedornot2), ]
p <- ggplot(data=zscoreWeighted3, aes(UMAP1, UMAP2, color=seedornot2)) + geom_point(size=0.05) +
    pretty_plot() + xlab("UMAP 1") + ylab("UMAP 2") + theme_void()  + theme(legend.position = "none")
ggsave(p, file = "test_seed461_2.png", width = 5, height = 5)

p <- ggplot(data=zscoreWeighted3, aes(UMAP1, UMAP2, color=cell_cluster2)) + geom_point(size=0.05) +
    pretty_plot() + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = "test_cellcluster_3.png", width = 5, height = 5)

p <- ggplot(data=zscoreWeighted2, aes(UMAP1, UMAP2, color=cell_cluster2)) + geom_point(size=0.05) +
    pretty_plot() + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = "test_cellcluster_4.png", width = 5, height = 5)

df <- getEmbedding(projncov, embedding = "UMAPHarmony", returnDF = TRUE)
colnames(df) <- c("real_umap1", "real_umap2")
# rownames(df)
df <- df[match(rownames(xtest), rownames(df)), ]
xtest2 <- data.frame(xtest, df)
p <- ggplot(data=xtest2, aes(real_umap1, real_umap2, color=BioCluster)) + geom_point(size=0.05) +
    pretty_plot() + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = "test_cellcluster_4.png", width = 5, height = 5)


df <- getEmbedding(projncov, embedding = "UMAPHarmony", returnDF = TRUE)
colnames(df) <- c("real_umap1", "real_umap2")
# rownames(df)
df <- df[match(rownames(zscoreWeighted3), rownames(df)), ]
zscoreWeighted3 <- data.frame(zscoreWeighted3, df)
p <- ggplot(data=zscoreWeighted3, aes(real_umap1, real_umap2, color=cell_cluster2)) + geom_point(size=0.05) +
    pretty_plot() + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p, file = "test_cellcluster_5.png", width = 5, height = 5)








