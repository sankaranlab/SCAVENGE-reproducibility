####------------------------------------------------------------------------------------------------------
# this is used for reproducibility of fig5 and corresponding supplementary figs and tables
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

zscore<- baso_df$Zscore
seed_idx <- seedindex(baso_df$Zscore, 0.05)
scale_factor <- cal_scalefactor(z_score=baso_df$Zscore, 0.01)

### Network propagation  
np_score <- randomWalk_sparse(intM=mutualknn30, queryCells=rownames(mutualknn30)[seed_idx], gamma=0.05)
# **Trait relevant score (TRS) with scaled and normalized**  
# A few cells are singletons are removed from further analysis
omit_idx <- np_score==0
sum(omit_idx)
meta_df <- meta_df[!omit_idx, ]
seed_idx <- seed_idx[!omit_idx]
mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
np_score <- np_score[!omit_idx]
TRS <- np_score %>% capOutlierQuantile(., 0.99) %>% max_min_scale
TRS <- TRS * scale_factor
zscoreWeighted3 <- data.frame(zscoreWeighted2[!omit_idx, ], TRS, np_score, seed_idx)
write.csv(zscoreWeighted3, "ALL_table.csv", row.names=F)

# compare the difference between before and after scavenge
## scatterplot 
blueYellow = c("#352A86", "#352A86", "#352A86", "#343DAE", "#343DAE", "#0262E0", "#0262E0", "#1389D2", "#1389D2", "#2DB7A3", "#2DB7A3", "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")
p <- ggplot(data=zscoreWeighted3, aes(UMAP1, UMAP2, color=Zscore)) + geom_point(size=0.2, na.rm = TRUE, alpha = 0.4) + 
    scale_color_gradientn(colors = blueYellow) + scale_alpha()+
    pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
p
p <- ggplot(data=zscoreWeighted3, aes(UMAP1, UMAP2, color=TRS)) + geom_point(size=0.2, na.rm = TRUE, alpha = 0.4) + 
    scale_color_gradientn(colors = blueYellow) + scale_alpha()+
    pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
p

violin_color=c("#7295e3", "#62b93c", "#b05cd3", "#47b65d", "#cd4dad", "#b8b939", "#636ad9", "#df9529", "#814ea3", "#618929", "#da417f", "#52c297", "#dc404f", "#42c0c7", "#d8592c", "#50a4d3", "#ae392f", "#2f856a", "#c58cd5", "#8eba6e", "#9b4577", "#3d7f46", "#de81ac", "#9d8428", "#5767a8", "#e18f5c", "#6b6e2e", "#e37f7d", "#c9ab69", "#a64e5a", "#9b5d2d")
pp <- ggplot(data=zscoreWeighted3,  aes(x=cell_cluster2, y=Zscore))  +
        geom_boxplot(aes(fill=cell_cluster2, color=cell_cluster2), outlier.shape=NA, width=0.75) + scale_fill_manual(values=violin_color) + scale_color_manual(values=violin_color)+
        ylab("Trait relevant score (TRS)") + xlab("") + ggtitle(trait_name) + ylim(c(-2, 4))+
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4)) +
        stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + theme(legend.position = "none")
pp
pp <- ggplot(data=zscoreWeighted3,  aes(x=cell_cluster2, y=TRS))  +
        geom_boxplot(aes(fill=cell_cluster2, color=cell_cluster2), outlier.shape=NA, width=0.75) + scale_fill_manual(values=violin_color) + scale_color_manual(values=violin_color)+
        ylab("Trait relevant score (TRS)") + xlab("") + ggtitle(trait_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4)) +
        stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + theme(legend.position = "none")
pp



trajectory <- c("HSC/MPP", "LMPP", "CLP", "Pro-B", "Pre-B", "Naive B", "Memory B", "Plasma cell")
ALL_trs3 <- zscoreWeighted3[zscoreWeighted3$cell_cluster %in% trajectory, ]
ALL_trs3 <- ALL_trs3[ALL_trs3$UMAP1<2, ]  # remove 3 cells
p <- ggplot(data=ALL_trs3, aes(UMAP1, UMAP2, color=TRS)) + geom_point(size=0.8, na.rm = TRUE, alpha = 0.65) + 
    scale_color_gradientn(colors = horizon) + scale_alpha()+
    pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
p



out <- readRDS("Aligned-Trajectory.rds")
trajectory <- c("HSC/MPP", "LMPP", "CLP", "Pro-B", "Pre-B", "Naive B", "Memory B", "Plasma B")

pseudo_outdf <- out[[1]][out[[1]]$Group %in% trajectory, ] # 14808
# range(pseudo_outdf$x)
# # [1] -12.436569   9.259022
# range(pseudo_outdf$y)
# # [1] -3.745459  9.509473

# dim(ALL_trs4)
# # [1] 14297    21
# range(ALL_trs3$UMAP1)
# # [1] -12.4365691  -0.2068717
# range(ALL_trs3$UMAP2)
# # [1] -2.547125  9.509473
pseudo_outdf <- pseudo_outdf[pseudo_outdf$x<2, ] # 14804
pseudo_outdf <- pseudo_outdf[match(ALL_trs3$cell_name, rownames(pseudo_outdf)), ]
ALL_trs3 <- data.frame(ALL_trs3, pseudotime=pseudo_outdf$pseudotime)
ggplot(ALL_trs3, aes(UMAP1, UMAP2, color=pseudotime)) + geom_point(size=0.8, na.rm = TRUE) +
    pretty_plot() + scale_color_gradientn(colors = jdb_palette("blue_cyan")) + scale_alpha() +
    geom_path(data=data.frame(out[[2]]), aes(x, y, color=NULL), 
        arrow = arrow(type = "open", angle = 20, length = unit(0.1, "inches")), size = 3)

# --- boxplot
trajectory <- c("HSC/MPP", "LMPP", "CLP", "Pro-B", "Pre-B", "Naive B", "Memory B", "Plasma cell")
ALL_trs3$cell_cluster2 <- factor(ALL_trs3$cell_cluster2, levels=rev(trajectory))

colfunc <- colorRampPalette(jdb_palette("blue_cyan"))
# colfunc(8)
pp <- ggplot(data=ALL_trs3,  aes(x=cell_cluster2, y=TRS))  +
    geom_boxplot(aes(fill=cell_cluster2, color=cell_cluster2), outlier.shape=NA, width=0.75, lwd=.5) + 
    scale_fill_manual(values=rev(colfunc(8))) + 
    scale_color_manual(values=rev(colfunc(8))) +
    ylab("Trait relevant score (TRS)") + xlab("") + ggtitle(trait_name) +
    theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4)) +
    coord_flip() +
    stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + theme(legend.position = "none")
pp

# scatter plot
pp <- ggplot(data=ALL_trs3,  aes(x=pseudotime, y=TRS, color=pseudotime))  +
    geom_point(size=0.8, na.rm = TRUE) + geom_smooth(colour = "black") + 
    pretty_plot() + scale_color_gradientn(colors=jdb_palette("blue_cyan")) +
    ylab("Trait relevant score (TRS)") + xlab("Pseudo time") + ggtitle(trait_name) +
    theme(axis.title.y=element_blank()) 
pp


# --------------------------------------------
# --------------------------------------------
# --------------------------------------------
# motif 
dev <- readRDS("chromVAR_Heme_All_SummarizedExperiment.final.rds")
motif_data <- dev@assays$data$z
motif_data <- motif_data[, colnames(motif_data) %in% ALL_trs3$cell_name ]  # 1764 14804
identical(ALL_trs3$cell_name, colnames(motif_data))
[1] FALSE
sample_idx <- match(ALL_trs3$cell_name, colnames(motif_data)) 
motif_data2 <- motif_data[, sample_idx]
identical(ALL_trs3$cell_name, colnames(motif_data2))

mycores=1
res_cor.test_list <- parallel::mclapply(1:nrow(motif_data2), mc.cores = mycores, function(i){
        xx <- cor.test(motif_data2[i, ] %>% unlist, ALL_trs3$zscore_topseed_np_scaled)
        res_cor.test = c(xx$estimate, xx$p.value)
        if (i %% 100 == 0) {print(i)}
        return(res_cor.test)
    }
)
res_cor.test_df <- as.data.frame(res_cor.test_list) %>% t %>% as.data.frame
rownames(res_cor.test_df) <- rownames(motif_data2)
colnames(res_cor.test_df) <- c("rho", "p.value")

res_cor.test_df$fdr <- p.adjust(res_cor.test_df[, 2], method = "bonferroni")
res_cor.test_df$log10fdr <- -log10(res_cor.test_df$fdr)

# res_cor.test_df$log10fdr[res_cor.test_df[, 1]<0] <- -res_cor.test_df$log10fdr[res_cor.test_df[, 1]<0]
res_cor.test_df$posorneg <- res_cor.test_df[, 1]>0
res_cor.test_df <- res_cor.test_df[order(0-res_cor.test_df$rho), ]
res_cor.test_df$order <- 1:nrow(res_cor.test_df)


load("cor-trs-Motifs_jaspar2018_sorted_Btrajectory.rda")
pos_tf <- c("PAX5_1131", "ID4_110", "EBF1_101", "TCF3_38", "ZEB1_231", "IRF4_963", "MEF2C_980", "BCL11A_284")
neg_tf <- c("NFE2_165", "GATA3_593", "CEBPE_146", "JUN_206", "RUNX1_1182", "HOXA9_626", "ERG_496", "MYB_995")

pos_tf_idx <- str_split(pos_tf, "_", simplify=T)[, 2] %>% as.numeric
neg_tf_idx <- str_split(neg_tf, "_", simplify=T)[, 2] %>% as.numeric
pos_tf_name <- str_split(pos_tf, "_", simplify=T)[, 1]
neg_tf_name <- str_split(neg_tf, "_", simplify=T)[, 1]
tf <- c(pos_tf, neg_tf)
tf_name <- c(pos_tf_name, neg_tf_name)

manual_selected_top <- res_cor.test_df[tf, ]
manual_selected_top <- data.frame(manual_selected_top, tf_name=tf_name)
my_col <- c("#007AB7", "#007AB7", "#208DC3", "#9ACDE7", "#E3E3E3", "#E3E3E3", "#FCEDAA", "#F8D32F", "#EDAD0B", "#EDAD0B")

p <- ggplot(res_cor.test_df, aes(x = order, y = rho, fill = rho)) + 
            geom_col(position = "identity") + pretty_plot() + scale_fill_gradientn(colors = my_col)+
            ggrepel::geom_label_repel(data=manual_selected_top, aes(order, rho, label = tf_name), color = 'black', fill="white",
                        size = 1.5, 
                        point.padding = 0, # additional padding around each point
    min.segment.length = 0, # draw all line segments
    max.time = 1, max.iter = 1e5, # stop after 1 second, or after 100,000 iterations
    box.padding = 0.3 # additional padding around each text label
    ) + xlab("Correlation between TF motif activaties and ALL-TRS across B-trajectory cells") +ylab("Spearman correlation")+ pretty_plot(fontsize = 10)
p




