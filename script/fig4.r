####------------------------------------------------------------------------------------------------------
# this is used for reproducibility of fig4 and corresponding supplementary figs and tables
# scavenge analysis of covid19-b1 traits using a covid PBMC scATAC dataset
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
    zscoreWeighted2$Zscore[is.na(baso_df$Zscore)] <- 0      
    baso_df$Zscore[is.na(baso_df$Zscore)] <- 0
}

seed_idx <- seedindex(baso_df$Zscore, 0.05)
scale_factor <- cal_scalefactor(z_score=baso_df$Zscore, 0.01)

### Network propagation  
np_score <- randomWalk_sparse(intM=mutualknn30, rownames(mutualknn30)[seed_idx], gamma=0.05)
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
# write.csv(zscoreWeighted3, "covid19b1_table.csv", row.names=F)

# compare the difference between before and after scavenge
## scatterplot 
Photinia_fraseri_B <- c("#57612f", "#b5b335", "#cecc2a", "#f8f39e", "#fdfadd", "#fddedb", "#f7a9a8", "#e54e4d", "#a32122")
p <- ggplot(data=zscoreWeighted3, aes(UMAP1, UMAP2, color=Zscore)) + geom_point(size=0.2, na.rm = TRUE, alpha = 0.4) + 
    scale_color_gradientn(colors = Photinia_fraseri_B) + scale_alpha()+
    pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")

p <- ggplot(data=zscoreWeighted3, aes(UMAP1, UMAP2, color=TRS)) + geom_point(size=0.2, na.rm = TRUE, alpha = 0.4) + 
    scale_color_gradientn(colors = Photinia_fraseri_B) + scale_alpha()+
    pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")


violin_color=c("#d8b7ec", "#349615", "#6b30b7", "#adcc25", "#2347ba", "#2ecf6b", "#db1d9d", "#85d6a8", "#966ef7", "#dc960f", "#307ce2", "#e5bc54", "#2d6da9", "#f33d45", "#267e68", "#b30c2f", "#b2d187", "#f87cd7", "#b45c0c", "#dc90e9", "#c89b74", "#a61b62", "#f67c71", "#653c69", "#925158")
pp <- ggplot(data=zscoreWeighted3,  aes(x=cell_cluster2, y=Zscore))  +
        geom_boxplot(aes(fill=cell_cluster2, color=cell_cluster2), outlier.shape=NA, width=0.75) + scale_fill_manual(values=violin_color) + scale_color_manual(values=violin_color)+
        ylab("Trait relevant score (TRS)") + xlab("") + ggtitle(trait_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4)) +
        stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + theme(legend.position = "none")
pp
pp <- ggplot(data=zscoreWeighted3,  aes(x=cell_cluster2, y=TRS))  +
        geom_boxplot(aes(fill=cell_cluster2, color=cell_cluster2), outlier.shape=NA, width=0.75) + scale_fill_manual(values=violin_color) + scale_color_manual(values=violin_color)+
        ylab("Trait relevant score (TRS)") + xlab("") + ggtitle(trait_name) +
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) + pretty_plot(fontsize = 10) + theme(axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=.4)) +
        stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + theme(legend.position = "none")
pp





covidb1_trs <- zscoreWeighted3
trait_name="COVID19-B1"
library(ggpubr)
# Visualize: Specify the comparisons you want
my_comparisons <- list( c("H", "M"), c("H", "S"), c("M", "S") )
p <- ggboxplot(covidb1_trs, x = "Sample_state", y = "TRS",
          fill = "Sample_state", palette = "npg", outlier.shape = NA,
          order = c("H", "M", "S")) + ylim(0, 10) +
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

# covid19hg_B1_df <- covidb1_trs
# save(covid19hg_B1_df, file="/broad/sankaranlab/fyu/varSC/data/20210919-scatac_cleandata4/mono_you_2021_covidpbmc/covid_b1_trs_df.rda")


covidb1_permu <- get_sigcell_simple(knn_sparse_mat=mutualknn30, seed_idx=covidb1_trs$seed_idx, topseed_npscore=covidb1_trs$np_score, permutation_times=1000, true_cell_significance=0.05, rda_output=F, mycores=8, rw_gamma=0.05)
covidb1_mat2 <- data.frame(covidb1_trs, true_cell_top_idx=covidb1_permu$true_cell_top_idx)

baso_df2 <- covidb1_mat2
baso_df2$true_cellplot = "Depleted cells"

baso_df2$true_cellplot[baso_df2$true_cell_top_idx] = "Enriched cells"
baso_df2$true_cellplot <- factor(baso_df2$true_cellplot, levels=c("Enriched cells", "Depleted cells"))
table(baso_df2$true_cellplot)

baso_df3 <- data.frame(baso_df2$cell_cluster, baso_df2$true_cellplot)
colnames(baso_df3) <- c("cell_cluster", "true_cellplot")
library(reshape2)
baso_df3_counts <- melt(table(baso_df3))
names(baso_df3_counts) <- names(baso_df3)
colnames(baso_df3_counts)[ncol(baso_df3_counts)] <- "cell_count"
baso_df3_counts %>% head
baso_df3_counts$cell_num <- table(baso_df2$cell_cluster)[match(as.character(baso_df3_counts$cell_cluster), names(table(baso_df2$cell_cluster)))]
baso_df3_counts$percent <- baso_df3_counts$cell_count/baso_df3_counts$cell_num
temp_df <- baso_df3_counts[baso_df3_counts$true_cellplot=="Enriched cells", ]
baso_df3_counts$cell_cluster <- factor(as.character(baso_df3_counts$cell_cluster), levels=rev(temp_df$cell_cluster[order(temp_df$percent)]))
mean(baso_df3_counts[2:15,5])
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
pp


##### look at cd14+ cells
baso_df2_mono14 = baso_df2[baso_df2$cell_cluster == "CD14+ Monocytes", ] 

table(baso_df2_mono14[, c("Sample_state", "true_cellplot")]) %>% chisq.test


odds.ratio(table(baso_df2_mono14[, c("Sample_state", "true_cellplot")]), conf.level = 0.95)

mosaicplot(table(baso_df2_mono14[, c("Health_state", "true_cellplot")]), col=c("firebrick", "goldenrod1"), cex.axis = 1.2, sub = "Condition", dir = c("h","v"), ylab = "Relative frequency")
mosaicplot(table(baso_df2_mono14[, c("Health_state", "seed_idx")]), col=c("firebrick", "goldenrod1"), cex.axis = 1.2, sub = "Condition", dir = c("h","v"), ylab = "Relative frequency")


library(epitools)
health_state_tab <- table(baso_df2_mono14[, c("Health_state", "true_cellplot")])
oddsratio(health_state_tab, method = "wald")$measure[-1,]
sample_state_tab <- table(baso_df2_mono14[, c("Sample_state", "true_cellplot")])
oddsratio(sample_state_tab[c(2, 1),], method = "wald")$measure[-1,]





baso_df2_mono14$true_cellplot = "Depleted cells"
baso_df2_mono14$true_cellplot[baso_df2_mono14$true_cell_top_idx] = "Enriched cells"
baso_df2_mono14$true_cellplot <- factor(baso_df2_mono14$true_cellplot, levels=c("Enriched cells", "Depleted cells"))
trait_name="covid19hg_B1"
mycolor=c("#343a40", "#dee2e6")
# mono14
ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=true_cellplot)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_manual(values=mycolor)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")





# ------------------------------------------------------------------------
# TF motif

library(ArchR)
addArchRThreads(threads = 2) 
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
motif_data <- proj_MotifMatrix@assays@data$z

identical(rownames(baso_df2), colnames(motif_data))
# [1] FALSE
sample_idx <- match(rownames(baso_df2), colnames(motif_data)) 
motif_data2 <- motif_data[, sample_idx]
identical(rownames(baso_df2), colnames(motif_data2))
# [1] TRUE

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

manual_selected_top <- rbind.data.frame(
                                rbind.data.frame(head(res_t.test_df, 50), tail(res_t.test_df, 50))[c(1, 4, 6, 9, 12, 14, 15, 50,    51, 53, 77, 78, 86, 94, 95, 98, 99), ], 
                                res_t.test_df[c(300, 322), ]
                                )
manual_selected_top$tf_name <- str_split(rownames(manual_selected_top), "_", simplify=T)[, 1]
p <- ggplot(res_t.test_df, aes(x = order, y = log10fdr, fill = posorneg)) +
            geom_col(position = "identity") +
            scale_fill_manual(values = c("#f37010", "#3893c8"), guide = FALSE) + 
            ggrepel::geom_text_repel(data=eithertop10, aes(order, log10fdr, label = tf_name), color = 'black',
                        size = 2, max.overlaps=300) + xlab("TF motif variable enrichment between trait- relevant and non-relevant cells") +ylab("-log10(Adjusted P value)")+ pretty_plot(fontsize = 10)

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
# scatterplot
tf_name="CEBPB_311"
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$CEBPB_311 <- motif_box_plot$score[match_idx]
baso_df2_mono14$CEBPB_311[baso_df2_mono14$CEBPB_311>quantile(baso_df2_mono14$CEBPB_311, 0.99)] <- quantile(baso_df2_mono14$CEBPB_311, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=CEBPB_311)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=viridis)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
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
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$SPI1_147 <- motif_box_plot$score[match_idx]
baso_df2_mono14$SPI1_147[baso_df2_mono14$SPI1_147>quantile(baso_df2_mono14$SPI1_147, 0.99)] <- quantile(baso_df2_mono14$SPI1_147, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=SPI1_147)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=viridis)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
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
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$SPIB_16 <- motif_box_plot$score[match_idx]
baso_df2_mono14$SPIB_16[baso_df2_mono14$SPIB_16>quantile(baso_df2_mono14$SPIB_16, 0.99)] <- quantile(baso_df2_mono14$SPIB_16, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=SPIB_16)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=viridis)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
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
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$JUN_44 <- motif_box_plot$score[match_idx]
baso_df2_mono14$JUN_44[baso_df2_mono14$JUN_44>quantile(baso_df2_mono14$JUN_44, 0.99)] <- quantile(baso_df2_mono14$JUN_44, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=JUN_44)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=viridis)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")

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
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$EOMES_274 <- motif_box_plot$score[match_idx]
baso_df2_mono14$EOMES_274[baso_df2_mono14$EOMES_274>quantile(baso_df2_mono14$EOMES_274, 0.99)] <- quantile(baso_df2_mono14$EOMES_274, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=EOMES_274)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=blueYellow)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
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
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$TBX21_153 <- motif_box_plot$score[match_idx]
baso_df2_mono14$TBX21_153[baso_df2_mono14$TBX21_153>quantile(baso_df2_mono14$TBX21_153, 0.99)] <- quantile(baso_df2_mono14$TBX21_153, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=TBX21_153)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=blueYellow)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
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
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$RUNX3_144 <- motif_box_plot$score[match_idx]
baso_df2_mono14$RUNX3_144[baso_df2_mono14$RUNX3_144>quantile(baso_df2_mono14$RUNX3_144, 0.99)] <- quantile(baso_df2_mono14$RUNX3_144, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=RUNX3_144)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=blueYellow)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
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
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$TCF7L2_58 <- motif_box_plot$score[match_idx]
baso_df2_mono14$TCF7L2_58[baso_df2_mono14$TCF7L2_58>quantile(baso_df2_mono14$TCF7L2_58, 0.99)] <- quantile(baso_df2_mono14$TCF7L2_58, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=TCF7L2_58)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=blueYellow)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
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
# scatterplot
match_idx <- match(rownames(baso_df2_mono14), rownames(motif_box_plot))
baso_df2_mono14$TCF4_307 <- motif_box_plot$score[match_idx]
baso_df2_mono14$TCF4_307[baso_df2_mono14$TCF4_307>quantile(baso_df2_mono14$TCF4_307, 0.99)] <- quantile(baso_df2_mono14$TCF4_307, 0.99)
# mono14
p <- ggplot(data=baso_df2_mono14, aes(UMAP1, UMAP2, color=TCF4_307)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=blueYellow)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")


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
tf_name="CEBPB_202"
match_idx <- match(rownames(baso_df2_Late.Ery), rownames(motif_box_plot))
baso_df2_Late.Ery$CEBPB_202 <- motif_box_plot$score[match_idx]
baso_df2_Late.Ery$CEBPB_202[baso_df2_Late.Ery$CEBPB_202>quantile(baso_df2_Late.Ery$CEBPB_202, 0.99)] <- quantile(baso_df2_Late.Ery$CEBPB_202, 0.99)
# Late.Ery
p <- ggplot(data=baso_df2_Late.Ery, aes(UMAP1, UMAP2, color=CEBPB_202)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=Jasminium_polyanthum_A)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")

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
# --------------------------------------------------------
# scatterplot
tf_name="IRF8_964"
match_idx <- match(rownames(baso_df2_Late.Ery), rownames(motif_box_plot))
baso_df2_Late.Ery$IRF8_964 <- motif_box_plot$score[match_idx]
baso_df2_Late.Ery$IRF8_964[baso_df2_Late.Ery$IRF8_964>quantile(baso_df2_Late.Ery$IRF8_964, 0.99)] <- quantile(baso_df2_Late.Ery$IRF8_964, 0.99)
# Late.Ery
p <- ggplot(data=baso_df2_Late.Ery, aes(UMAP1, UMAP2, color=IRF8_964)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors=Jasminium_polyanthum_A)  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")

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
# --------------------------------------------------------
# scatterplot
tf_name="GATA2_600"
match_idx <- match(rownames(baso_df2_Late.Ery), rownames(motif_box_plot))
baso_df2_Late.Ery$GATA2_600 <- motif_box_plot$score[match_idx]
baso_df2_Late.Ery$GATA2_600[baso_df2_Late.Ery$GATA2_600>quantile(baso_df2_Late.Ery$GATA2_600, 0.99)] <- quantile(baso_df2_Late.Ery$GATA2_600, 0.99)
# Late.Ery
p <- ggplot(data=baso_df2_Late.Ery, aes(UMAP1, UMAP2, color=GATA2_600)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors = jdb_palette("brewer_marine"))  +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")


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
# --------------------------------------------------------
# scatterplot
tf_name="GATA1_588"
match_idx <- match(rownames(baso_df2_Late.Ery), rownames(motif_box_plot))
baso_df2_Late.Ery$GATA1_588 <- motif_box_plot$score[match_idx]
baso_df2_Late.Ery$GATA1_588[baso_df2_Late.Ery$GATA1_588>quantile(baso_df2_Late.Ery$GATA1_588, 0.99)] <- quantile(baso_df2_Late.Ery$GATA1_588, 0.99)
# Late.Ery
p <- ggplot(data=baso_df2_Late.Ery, aes(UMAP1, UMAP2, color=GATA1_588)) + geom_point(size=1, na.rm = TRUE, alpha = 0.9) + scale_color_gradientn(colors = jdb_palette("brewer_marine")) +
        pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")



# ---------------------------------------------
# gene subset of Mono
# ---------------------------------------------
# Mono 14
library(ArchR)
addArchRThreads(threads = 2) 
addArchRGenome("hg19")
proj <- loadArchRProject(path = "/broad/sankaranlab/fyu/varSC/data/you_2021_covidpbmc/COVID19-fragment/covidpbmc_sample11")

idx_Mono <- match(rownames(baso_df2_mono14), getCellNames(proj))
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
proj_mono14@cellColData$true_cellplot_1k[rank(proj_mono14@cellColData$TRS) <= 1000] <- "Depleted"
proj_mono14@cellColData$true_cellplot_1k[rank(proj_mono14@cellColData$TRS) >= (length(proj_mono14@cellColData$true_cellplot_1k) -999)] <- "Enriched"
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


p <- plotBrowserTrack(proj_mono14, groupBy = "true_cellplot_1k", useGroups=c("Depleted", "Enriched"), 
                        geneSymbol=c("HMMR", "NUDCD2", "CCNG1", "MAT2B"), log2Norm=T, tileSize=1000, useMatrix="GeneScoreMatrix", plotSummary=c("bulkTrack", "scTrack", "geneTrack", "featureTrack"), scTileSize=1, scCellsMax = 1000)
plotPDF(p, name="test_gb_chr5_true_cellplot_500k_gene.pdf", ArchRProj = proj_mono14)

p <- plotBrowserTrack(proj_mono14, groupBy = "true_cellplot", geneSymbol = "CCNG1", upstream = 500000, downstream = 500000)
plotPDF(p, name="test_gb_CCNG1.pdf", ArchRProj = proj_mono14)
p <- plotBrowserTrack(proj_mono14, groupBy = "true_cellplot", geneSymbol = "MAT2B", upstream = 500000, downstream = 500000)
plotPDF(p, name="test_gb_MAT2B.pdf", ArchRProj = proj_mono14)

