####------------------------------------------------------------------------------------------------------
# this is used for reproducibility of fig3 and corresponding supplementary figs and tables
# scavenge analysis of 22 blood cell traits using two independent hematopoiesis scATAC datasets
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


# -----------
# # blood22
# -----------
# -------------------------------------------------------
# Granja_2019_heme (hg19)
# ukbb_wDEV, zscoreWeighted2, meta_df
load("blood22_Granja_2019_heme_gscore.rda")
load("Granja_2019_heme_knngraph.rda")
zscoreWeighted2$cell_cluster2 <- zscoreWeighted2$cell_cluster
celltypeanno <- c("HSC", "Early.Eryth", "Late.Eryth", "Early.Baso", "CMP.LMPP", "CLP.1", "GMP", "GMP.Neut", "pDC", "cDC", "CD14.Mono.1", "CD14.Mono.2", "CLP.2", "Pre.B","B", "Plasma", "CD8.N", "CD4.N1", "CD4.N2", "CD4.M", "CD8.EM", "CD8.CM","NK")
zscoreWeighted2$cell_cluster2 <- factor(zscoreWeighted2$cell_cluster2, levels=celltypeanno)

# baso_df, baso_trs , eo_df , eo_trs, hct_df, hct_trs, hgb_df, hgb_trs, hlr_df, hlr_trs, irf_df, irf_trs, lymph_df , lymph_trs, mch_df, mch_trs, mchc_df, mchc_trs , mcv_df, mcv_trs, mono_df, mono_trs , mpv_df, mpv_trs, mrv_df, mrv_trs, mscv_df, mscv_trs, neut_df, neut_trs, pct_df, pct_trs, pdw_df, pdw_trs, plt_df, plt_trs, rbc_df, rbc_trs, rdw_cv_df, rdw_cv_trs, ret_df, ret_trs, wbc_df, wbc_trs
load("blood22_Granja_2019_heme_nptrs.rda")

trait_name22 <- c("baso", "eo", "hct", "hgb", "hlr", "irf", "lymph", "mch", "mchc", "mcv", "mono", "mpv", "mrv", "mscv", "neut", "pct", "pdw", "plt", "rbc", "rdw_cv", "ret", "wbc")
temp_mat <- matrix(0, nrow=23, ncol=length(trait_name22))
rownames(temp_mat) <- celltypeanno
colnames(temp_mat) <- trait_name22
blood22_trs_mean_mat <- blood22_trs_median_mat <- temp_mat
trait_list <- list(baso_trs, eo_trs, hct_trs, hgb_trs, hlr_trs, irf_trs, lymph_trs, mch_trs, mchc_trs, mcv_trs, mono_trs, mpv_trs, mrv_trs, mscv_trs, neut_trs, pct_trs, pdw_trs, plt_trs, rbc_trs, rdw_cv_trs, ret_trs, wbc_trs)
for (i in 1:22){
    
    blood22_trs_mean_mat[, i] <- with(trait_list[[i]], tapply(X = zscore_topseed_np_scaled, 
                         INDEX = cell_cluster2, 
                         FUN = mean 
                         ))
    blood22_trs_median_mat[, i] <- with(trait_list[[i]], tapply(X = zscore_topseed_np_scaled, 
                         INDEX = cell_cluster2, 
                         FUN = median 
                         ))
}
xx <- data.frame(matrix(NA, nrow=nrow(baso_trs), ncol=22))
for (i in 1:22){
    xx[, i] <- trait_list[[i]]$`zscore_topseed_np_scaled`
}
colnames(xx) <- trait_name22
blood22_trs_table <- data.frame(baso_trs[, c(1, 3, 4, 7)], xx)
write.csv(blood22_trs_table, "blood22_trs_table.csv", row.names=F)

# local R
load("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/blood22_Granja_2019_heme/blood22_Granja_2019_heme_dataheatmap.rda")
library("stringr")
library("circlize")
library(ComplexHeatmap)

cell_types <- c("Platelet", 
                "Red cell",                
                "Myeloid/Compound white cell",               
                "Lymphoid white cell"
                )
cell_relevance_anno <- cell_types[c(3, 3, 2, 2, 2, 2, 4, 2, 2, 2, 3, 1, 2, 2, 3, 1, 1, 1, 2, 2, 2, 3)]
# cbind(colnames(blood22_trs_median_mat), cell_relevance_anno)

pdf(file="heatmap22bloodtraits.pdf", width=10.8, height=6)
ha = HeatmapAnnotation("anno" = cell_relevance_anno,
                       which="row",
                       col = list("anno" = c("Platelet" = "#4b994a", 
                                             "Red cell" = "#2e82b3", 
                                             "Myeloid/Compound white cell" = "#Ea9e01", 
                                             "Lymphoid white cell"= "#Bd3046")),
                       #gp = gpar(col = "white", lwd = .1),
                       border = F)
Heatmap(scale(blood22_trs_median_mat) %>% t, rect_gp = gpar(col = "white", lwd = .2), clustering_distance_rows = "pearson", column_names_side = "top",
        col=jdb_palette("brewer_violet"), row_dend_side = "left", column_dend_side = "bottom", row_km = 4, 
        right_annotation = ha, column_names_rot = 45, border = F, 
        column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 9)
)
dev.off()


# -------------------------------------------------------
# Satpathy_2019_heme (hg19)
load("blood22_Satpathy_2019_heme_gscore.rda")
load("Satpathy_2019_heme_knngraph.rda")
load("blood22_Satpathy_2019_heme_nptrs.rda")

trait_name22 <- c("baso", "eo", "hct", "hgb", "hlr", "irf", "lymph", "mch", "mchc", "mcv", "mono", "mpv", "mrv", "mscv", "neut", "pct", "pdw", "plt", "rbc", "rdw_cv", "ret", "wbc")
temp_mat <- matrix(0, nrow=31, ncol=length(trait_name22))
rownames(temp_mat) <- celltypeanno
colnames(temp_mat) <- trait_name22
blood22_trs_mean_mat <- blood22_trs_median_mat <- temp_mat
trait_list <- list(baso_trs, eo_trs, hct_trs, hgb_trs, hlr_trs, irf_trs, lymph_trs, mch_trs, mchc_trs, mcv_trs, mono_trs, mpv_trs, mrv_trs, mscv_trs, neut_trs, pct_trs, pdw_trs, plt_trs, rbc_trs, rdw_cv_trs, ret_trs, wbc_trs)
for (i in 1:22){
    
    blood22_trs_mean_mat[, i] <- with(trait_list[[i]], tapply(X = zscore_topseed_np_scaled, 
                         INDEX = cell_cluster2, 
                         FUN = mean 
                         ))
    blood22_trs_median_mat[, i] <- with(trait_list[[i]], tapply(X = zscore_topseed_np_scaled, 
                         INDEX = cell_cluster2, 
                         FUN = median 
                         ))
}
xx <- data.frame(matrix(NA, nrow=nrow(baso_trs), ncol=22))
for (i in 1:22){
    xx[, i] <- trait_list[[i]]$`zscore_topseed_np_scaled`
}
colnames(xx) <- trait_name22
blood22_trs_table <- data.frame(baso_trs[, c(1, 3, 4, 7)], xx)
write.csv(blood22_trs_table, "blood22_trs_table2.csv", row.names=F)

library("stringr")
library("circlize")
library(ComplexHeatmap)
cell_types <- c("Platelet", 
                "Red cell",                
                "Myeloid/Compound white cell",               
                "Lymphoid white cell"
                )
cell_relevance_anno <- cell_types[c(3, 3, 2, 2, 2, 2, 4, 2, 2, 2, 3, 1, 2, 2, 3, 1, 1, 1, 2, 2, 2, 3)]
# cbind(colnames(blood22_trs_median_mat), cell_relevance_anno)
pdf(file="heatmap22bloodtraits.pdf", width=10.8, height=6)
ha = HeatmapAnnotation("anno" = cell_relevance_anno,
                       which="row",
                       col = list("anno" = c("Platelet" = "#4b994a", 
                                             "Red cell" = "#2e82b3", 
                                             "Myeloid/Compound white cell" = "#Ea9e01", 
                                             "Lymphoid white cell"= "#Bd3046")),
                       #gp = gpar(col = "white", lwd = .1),
                       border = F)
Heatmap(scale(blood22_trs_median_mat) %>% t, rect_gp = gpar(col = "white", lwd = .2), clustering_distance_rows = "pearson", column_names_side = "top",
        col=jdb_palette("brewer_blue"), row_dend_side = "left", column_dend_side = "bottom", row_km = 4, 
        right_annotation = ha, column_names_rot = 45, border = F, 
        column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 9)
)
dev.off()








