####------------------------------------------------------------------------------------------------------
# this is used for reproducibility of fig2 and corresponding supplementary figs and tables
# scavenge analysis with simulated data and 10X pbmc scATAC datasets
####------------------------------------------------------------------------------------------------------

####----------------------------------
#  simulations with bulk atac
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# use real bulk count matrix from mpn study
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
library(Metrics)
setwd("/Users/fyu/Documents/GitHub/SCAVENGE-reproducibility/data")


simulate_scatac <- function(n_cells, which_celltypes, n_frags_per_cell = 1000, 
                            rate_noise = 0, seed = 100, shuffle = FALSE, bulk=xx, peaks=peaks){
  # Reproducibility
  set.seed(seed)
#   which_celltypes <- sort(which_celltypes)
  stopifnot(rate_noise < 1) 
  stopifnot(n_frags_per_cell > 100)
  n_peaks <- dim(bulk)[1]
  #--
  # Set up cell labels
  #--
  message("here0")
  if(length(n_cells) > 1){
    stopifnot(length(which_celltypes) == length(n_cells))
    
    # Generate cell labels
    cell_labels <- sapply(1:length(which_celltypes), function(i){
      rep(which_celltypes[i], n_cells[i])
    }) %>% unlist() %>% sort()
    
  } else {
#     n_groups <- length(which_celltypes)
#     cell_labels <- sort(rep(which_celltypes, n_cells*n_groups))
      cell_labels <- sort(rep(which_celltypes, n_cells))
  }
  final_names <- paste0(cell_labels, "_", as.character(1:length(cell_labels)))
  
  #-------------------
  # Simulate true data
  #-------------------
  message("here1")
  # Generate cell-type specific peaks
  lapply(which_celltypes, function(celltype){
    
    # Apply different rates per cell depending on group label for generating cell-type specific peaks
    n_cells_this_celltype <- sum(cell_labels == celltype)
    counts_celltype <- bulk[, celltype]
    
    # Define probabilities
    #                        Prob observting frag                Total number of fragments epxpected; the 0.5s are for two alleles that will be simulated/added later
    prob_per_peaks <- counts_celltype/sum(counts_celltype) * (n_frags_per_cell*0.5 * (1-rate_noise)) + ((rate_noise*n_frags_per_cell)/n_peaks*0.5) 
    
    # Cap probabilities at something sensible
    prob_per_peaks <- ifelse(prob_per_peaks > 0.9, 0.9, prob_per_peaks)
    
    # Represent the two haplotypes as two random draws
    mat1 <- (matrix(rbinom(n_peaks*n_cells_this_celltype, size = 1, prob = prob_per_peaks),
                    ncol = n_cells_this_celltype, byrow = FALSE) )
    mat2 <- (matrix(rbinom(n_peaks*n_cells_this_celltype, size = 1, prob = prob_per_peaks),
                    ncol = n_cells_this_celltype, byrow = FALSE) )
    
    mat <- mat1 + mat2
    Matrix(mat)
  }) %>% do.call(what = "cbind") -> sparse_matrix
  return(sparse_matrix)
  message("here2")
  colnames(sparse_matrix) <- final_names
  message("here2")
  peaknames = paste(peaks$V1,peaks$V2,peaks$V3,sep = "_")
  rownames(sparse_matrix) <- peaknames
  message("here2")
  sparse_matrix
}
strong_connect <- function(mknn){
    library(igraph)
    mknn_graph <- graph_from_adjacency_matrix(mknn, mode = "undirected", diag = F)
    clu <- components(mknn_graph)
    set.seed(9527)
    if (clu$no!=1){
        print(clu$no)
        link_chain <- tapply(1:nrow(mknn), clu[[1]], function(x) {sample(x, 1)})
        mknn[link_chain[1:(length(link_chain)-1)] %>% c, link_chain[2:length(link_chain)] %>% c]
        for (i in 1:(length(link_chain)-1)){
            mknn[(link_chain[1:(length(link_chain)-1)] %>% c)[i], (link_chain[2:length(link_chain)] %>% c)[i]] <- 1
            mknn[(link_chain[2:length(link_chain)] %>% c)[i], (link_chain[1:(length(link_chain)-1)] %>% c)[i]] <- 1
        }
    }
    return(mknn)
}



# Import and run bulk ATAC
peaksdf <- fread("26August2017_EJCsamples_allReads_250bp.bed") 
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("26August2017_EJCsamples_allReads_250bp.counts.txt")) 
mono_trait <- "finemappedtraits_hg19/mono.PP001.bed"
trait <- "mono"

# simulate two cell types with 500 each

test_q40_n200_c9 <- simulate_scatac(n_cells=500, which_celltypes=c("Mono", "NK"), n_frags_per_cell = 10000, 
                                    rate_noise = 0.3, seed = 9527, shuffle = FALSE, bulk=counts)
columns <- paste(rep(c("Mono", "NK"), each=500), rep(1:500, times=2), sep="_")
print(length(columns))
dimnames(test_q40_n200_c9)[[2]] = columns
idx_no0 <- rowSums(test_q40_n200_c9)!=0; sum(!idx_no0)
        # Create objects for g-chromVAR
SE_2 <- SummarizedExperiment(assays = list(counts = test_q40_n200_c9[idx_no0, ]),
                        rowData = peaks[idx_no0 ], 
                        colData = DataFrame(names = colnames(test_q40_n200_c9[idx_no0, ])))
SE_2 <- addGCBias(SE_2, genome = BSgenome.Hsapiens.UCSC.hg19)
mpns <- importBedScore(rowRanges(SE_2), mono_trait, colidx = 5)
        # Run g-chromVAR
bg <- getBackgroundPeaks(SE_2,niterations=200)
dev <- computeWeightedDeviations(SE_2, mpns, background_peaks = bg)
        # Reformat results
zscoreWeighted <- reshape2::melt(t(assays(dev)[["z"]]))
colnames(zscoreWeighted) <- c("Celltype","Trait","zscore")
zscoreWeighted$Celltype2 <- stringr::str_split(zscoreWeighted$Celltype, "_", simplify=T)[, 1]
zscoreWeighted$rank <- rank(-zscoreWeighted$zscore)
plot(1:length(zscoreWeighted$zscore), zscoreWeighted$zscore)

matSVD <- function(mat=test_q40_n200_c9[idx_no0, ], nComponents = 50, binarize = TRUE){
  
  #TF IDF LSI adapted from flyATAC
  cs <- Matrix::colSums(mat)
  if(binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1 
  }
  #Calc TF IDF
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/Matrix::colSums(mat))
  idf   <- as(log(1 + ncol(mat) / Matrix::rowSums(mat)), "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svd <- irlba::irlba(tfidf, nComponents, nComponents)
  svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  return(matSVD)
}

lsimat50=matSVD(mat=test_q40_n200_c9[idx_no0, ], nComponents = 30, binarize = TRUE) # 1000    50
# save(lsimat50, file="test_q40_n200_c9_lsimat50.rda")
mutualknn30 <- getmutualknn(lsimat50, 30)
mutualknn30 <- strong_connect(mutualknn30)
mutualknn30_graph <- graph_from_adjacency_matrix(mutualknn30, mode = "undirected", diag = F)

mutualknn5 <- getmutualknn(lsimat50, 5)
mutualknn5 <- strong_connect(mutualknn5)
mutualknn5_graph <- graph_from_adjacency_matrix(mutualknn5, mode = "undirected", diag = F)

mutualknn10 <- getmutualknn(lsimat50, 10)
mutualknn10 <- strong_connect(mutualknn10)
mutualknn10_graph <- graph_from_adjacency_matrix(mutualknn10, mode = "undirected", diag = F)

mutualknn20 <- getmutualknn(lsimat50, 20)
mutualknn20 <- strong_connect(mutualknn20)
mutualknn20_graph <- graph_from_adjacency_matrix(mutualknn20, mode = "undirected", diag = F)

mutualknn50 <- getmutualknn(lsimat50, 50)
mutualknn50 <- strong_connect(mutualknn50)
mutualknn50_graph <- graph_from_adjacency_matrix(mutualknn50, mode = "undirected", diag = F)




qq <- 0.05
seed_p0.05_top <- seedindex(zscoreWeighted$zscore, qq)
zscore_topseed_np <- randomWalk_sparse(intM=mutualknn30, queryCells=rownames(mutualknn30)[seed_p0.05_top], gamma=0.05)
# only top seed times z factor
zscore_topseed_np_scaled <- zscore_topseed_np %>% capOutlierQuantile(., 0.99) %>% max_min_scale
zscoreWeighted2 <- zscoreWeighted
zscoreWeighted2$zscore_topseed_np_scaled <- zscore_topseed_np_scaled
sum(zscoreWeighted2$zscore_topseed_np_scaled==0) # remove singletons isolated in the graph 8
# zscoreWeighted2 <- zscoreWeighted2[zscoreWeighted2$zscore_topseed_np_scaled!=0, ]
zscoreWeighted2$rank_np2 <- BiocGenerics::rank(-zscoreWeighted2$zscore_topseed_np_scaled, ties.method="first")
zscoreWeighted3 <- zscoreWeighted2[zscoreWeighted2$zscore_topseed_np_scaled !=0, ]
# comparsion of ranks 
# TRS 
pp <- ggplot(zscoreWeighted2, aes(x = rank_np2, y = 0.3, fill = Celltype2)) +
  geom_col() + pretty_plot() + 
  scale_fill_manual(values = c("#ec4a9a", "#dddddd")) +
  xlab("rank")
pp
ggsave(pp, file = paste0("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-1.pdf"), width = 5, height = 3)

# z score
pp <- ggplot(zscoreWeighted2, aes(x = rank, y = 0.3, fill = Celltype2)) +
  geom_col() + pretty_plot() + 
  scale_fill_manual(values = c("#ec4a9a", "#dddddd")) +
  xlab("rank")
pp
ggsave(pp, file = paste0("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-2.pdf"), width = 5, height = 3)

# look at different quarters
# donutplots
auc(zscoreWeighted2$Celltype2=="Mono", zscoreWeighted2$zscore_topseed_np_scaled)
# [1] 0.993
auc(zscoreWeighted2$Celltype2=="Mono", zscoreWeighted2$zscore)
# [1] 0.7993467
(quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank_np2, 4), table))


$`(0.002,250]`

Mono   NK # 0.996
 249    1 

$`(250,500]`

Mono   NK # 0.996
 237   13 

$`(500,750]`

Mono   NK # 0.948
  14  236 

$`(750,1e+03]`

 NK # 0
250

tpr= 97.2; fpr (4+24)/((4+24)+(234+228)) = 0.057; accuracy 97.2
# Hole size
hsize <- 2
(df1 <- data.frame(value=c(249, 1), group=c("Mono", "NK"), hsize))
df2 <- data.frame(value=c(237, 13), group=c("Mono", "NK"), hsize)
df3 <- data.frame(value=c(14, 236), group=c("Mono", "NK"), hsize)
df4 <- data.frame(value=c(0, 250), group=c("Mono", "NK"), hsize)
pp <- ggplot(df1, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave(pp, file = paste0("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-3.pdf"), width = 5, height = 3)

pp <- ggplot(df2, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave(pp, file = paste0("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-4.pdf"), width = 5, height = 3)
pp <- ggplot(df3, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave(pp, file = paste0("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-5.pdf"), width = 5, height = 3)
pp <- ggplot(df4, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave(pp, file = paste0("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-6.pdf"), width = 5, height = 3)


Mono   NK # 0.816
 204   46 

$`(251,500]`

Mono   NK # 0.628
 157   93 

$`(500,750]`

Mono   NK # 0.392
  98  152 

$`(750,1e+03]`

Mono   NK # 0.164
  41  209 

tpr, 72.2 ; auc 0.799; accuracy 72.2
(df1 <- data.frame(value=c(204, 46), group=c("Mono", "NK"), hsize))
df2 <- data.frame(value=c(157, 93), group=c("Mono", "NK"), hsize)
df3 <- data.frame(value=c(98, 152), group=c("Mono", "NK"), hsize)
df4 <- data.frame(value=c(41, 209), group=c("Mono", "NK"), hsize)
pp <- ggplot(df1, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave(pp, file = paste0("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-7.pdf"), width = 5, height = 3)
pp <- ggplot(df2, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave(pp, file = paste0("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-8.pdf"), width = 5, height = 3)
pp <- ggplot(df3, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave(pp, file = paste0("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-9.pdf"), width = 5, height = 3)
pp <- ggplot(df4, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave(pp, file = paste0("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-10.pdf"), width = 5, height = 3)

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# investigation of effects unbalanced cell number
# simulate two cell types with 500 each
xx <- list(c(100, 900), c(200, 800), c(300, 700), c(400, 600), c(500, 500), c(600, 400), c(700, 300), c(800, 200), c(900, 100))
auc_z <- auc_v <- rep(0, 9)
for (i in 1:9){
    message("----------", i, "----------")
    test_q40_n200_c9 <- simulate_scatac(n_cells=xx[[i]], which_celltypes=c("Mono", "NK"), n_frags_per_cell = 10000, 
                                    rate_noise = 0.3, seed = 9527, shuffle = FALSE, bulk=counts)
    columns <- paste(rep(c("Mono", "NK"), xx[[i]]), 1:1000, sep="_")
    dimnames(test_q40_n200_c9)[[2]] = columns
    idx_no0 <- rowSums(test_q40_n200_c9)!=0; sum(!idx_no0)
            # Create objects for g-chromVAR
    SE_2 <- SummarizedExperiment(assays = list(counts = test_q40_n200_c9[idx_no0, ]),
                            rowData = peaks[idx_no0 ], 
                            colData = DataFrame(names = colnames(test_q40_n200_c9[idx_no0, ])))
    SE_2 <- addGCBias(SE_2, genome = BSgenome.Hsapiens.UCSC.hg19)
    mpns <- importBedScore(rowRanges(SE_2), mono_trait, colidx = 5)
            # Run g-chromVAR
    bg <- getBackgroundPeaks(SE_2,niterations=200)
    dev <- computeWeightedDeviations(SE_2, mpns, background_peaks = bg)
            # Reformat results
    zscoreWeighted <- reshape2::melt(t(assays(dev)[["z"]]))
    colnames(zscoreWeighted) <- c("Celltype","Trait","zscore")
    zscoreWeighted$Celltype2 <- stringr::str_split(zscoreWeighted$Celltype, "_", simplify=T)[, 1]
    zscoreWeighted$rank <- rank(-zscoreWeighted$zscore)
    plot(1:length(zscoreWeighted$zscore), zscoreWeighted$zscore)
    lsimat50=matSVD(mat=test_q40_n200_c9[idx_no0, ], nComponents = 30, binarize = TRUE) # 1000    50

    mutualknn50 <- getmutualknn(lsimat50, 30)
    mutualknn50 <- strong_connect(mutualknn50)
    mutualknn50_graph <- graph_from_adjacency_matrix(mutualknn50, mode = "undirected", diag = F)


    qq <- 0.05
    seed_p0.05_top <- seedindex(zscoreWeighted$zscore, qq)
    zscore_topseed_np <- randomWalk_sparse(intM=mutualknn50, queryCells=rownames(mutualknn50)[seed_p0.05_top], gamma=0.05)
    # only top seed times z factor
    zscore_topseed_np_scaled <- zscore_topseed_np %>% capOutlierQuantile(., 0.99) %>% max_min_scale
    zscoreWeighted2 <- zscoreWeighted
    zscoreWeighted2$zscore_topseed_np_scaled <- zscore_topseed_np_scaled
    sum(zscoreWeighted2$zscore_topseed_np_scaled==0) # remove singletons isolated in the graph 8
    # zscoreWeighted2 <- zscoreWeighted2[zscoreWeighted2$zscore_topseed_np_scaled!=0, ]
    zscoreWeighted2$rank_np2 <- BiocGenerics::rank(-zscoreWeighted2$zscore_topseed_np_scaled, ties.method="first")
    zscoreWeighted3 <- zscoreWeighted2[zscoreWeighted2$zscore_topseed_np_scaled !=0, ]
    # (quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank_np2, 10), table))
    # print(quantile4tab)
    (quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank, 10), table))
    print(quantile4tab)
    auc_v[i] <- auc(zscoreWeighted2$Celltype2=="Mono", zscoreWeighted2$zscore_topseed_np_scaled)
    # [1] 0.9772018
    auc_z[i] <- auc(zscoreWeighted2$Celltype2=="Mono", zscoreWeighted2$zscore)
    # pp <- ggplot(zscoreWeighted2, aes(x = rank_np2, y = 0.3, fill = Celltype2)) +
    #     geom_col() + pretty_plot() + 
    #     scale_fill_manual(values = c("#ec4a9a", "#dddddd")) +
    #     xlab("rank")
    # ggsave(pp, file = paste0("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/robustness-unbalancelabel3-np-", i, ".pdf"), width = 5, height = 3)
    # pp <- ggplot(zscoreWeighted2, aes(x = rank, y = 0.3, fill = Celltype2)) +
    #     geom_col() + pretty_plot() + 
    #     scale_fill_manual(values = c("#ec4a9a", "#dddddd")) +
    #     xlab("rank")
    # ggsave(pp, file = paste0("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/robustness-unbalancelabel3-z-", i, ".pdf"), width = 5, height = 3)

}
> auc_z
[1] 0.7531556 0.7947062 0.7791810 0.7964000 0.7994320
[6] 0.7941792 0.7932762 0.8131438 0.8133111
> auc_v
[1] 0.7538833 0.9606594 0.9510048 0.9569417 0.9936080
[6] 0.9999958 0.8855143 0.9999938 0.9936222
c(0.7531556, 0.7947062, 0.7791810, 0.7964000, 0.7994320, 0.7941792, 0.7932762, 0.8131438, 0.8133111) %>% mean 
[1] 0.7929761
c(0.7538833, 0.9606594,  0.9510048,  0.9569417,  0.9936080,  0.9999958,  0.8855143,  0.9999938,  0.9936222) %>% mean 
[1] 0.9439137

tpr_v <- c(0.49, .905, .787, .933, .972, .998, .851, .999, .988)
fpr_v <- c(0.057, 0.024, 0.091, .045, .028, .0025, 0.347, .005, .11)

tpr_z <- c(0.34, .485, .563, .658, .722, .753, .819, .868, 0.93)
fpr_z <- c(0.073, 0.129, 0.189, 0.228, .278, .37, 0.423, .48, .62)

# c(2500, 2500, 5000, 5000, 7500, 7500, 10000, 10000, 25000, 25000, 50000, 50000)

TPR <- (matrix(c(tpr_z, tpr_v), ncol=2) %>% t %>% c)*100
FPR <- (matrix(c(fpr_z, fpr_v), ncol=2) %>% t %>% c)*100
df3 <- data.frame(reads=rep(seq(10, 90, 10), each=2), 
                  auROC=(matrix(c(auc_z, auc_v), ncol=2) %>% t %>% c)*100,
                  TPR=TPR,
                  FPR=FPR,
                  type=rep(c("Z score", "SCAVENGE"), 9))
df3$type <- factor(df3$type, levels=c("Z score", "SCAVENGE") )

p <- ggplot(df3, aes(x=reads, y=auROC, fill=type)) + 
   geom_bar(stat="identity", position=position_dodge())+ coord_cartesian(ylim=c(50, 100)) + 
   scale_fill_brewer(palette="Purples") + pretty_plot()
ggsave(p, file = "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-11.pdf", width = 4, height = 3)

p <- ggplot(df3, aes(x=reads, y=TPR, fill=type)) + 
   geom_bar(stat="identity", position=position_dodge())+ coord_cartesian(ylim=c(50, 100)) + 
   scale_fill_brewer(palette="PuRd") + pretty_plot()
ggsave(p, file = "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-11-2.pdf", width = 4, height = 3)

p <- ggplot(df3, aes(x=reads, y=FPR, fill=type)) + 
   geom_bar(stat="identity", position=position_dodge())+ coord_cartesian(ylim=c(0, 65)) + 
   scale_fill_brewer(palette="Greys") + pretty_plot()
ggsave(p, file = "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-11-3.pdf", width = 4, height = 3)



# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# investigation of different proportion of seed cells
test_q40_n200_c9 <- simulate_scatac(n_cells=500, which_celltypes=c("Mono", "NK"), n_frags_per_cell = 10000, 
                                    rate_noise = 0.3, seed = 9527, shuffle = FALSE, bulk=counts)
columns <- paste(rep(c("Mono", "NK"), each=500), rep(1:500, times=2), sep="_")
print(length(columns))
dimnames(test_q40_n200_c9)[[2]] = columns
idx_no0 <- rowSums(test_q40_n200_c9)!=0; sum(!idx_no0)

lsimat50=matSVD(mat=test_q40_n200_c9[idx_no0, ], nComponents = 30, binarize = TRUE) # 1000    50

mutualknn50 <- getmutualknn(lsimat50, 30)
mutualknn50 <- strong_connect(mutualknn50)
mutualknn50_graph <- graph_from_adjacency_matrix(mutualknn50, mode = "undirected", diag = F)

qq <- (1:50)/100
mm <- rep(0, 50) # tpr # tp/(tp+fn)
mm2 <- rep(0, 50) # fpr # fp/(fp+tn) tpr=fpr as cell number is balance in two cell types
auc_v <-  rep(0, 50) # auc
nn <- rep(0, 50)
aa <- list(); length(aa) <- 50
for (i in 1:length(qq)){
    seed_p0.05_top <- zscoreWeighted$rank<=(1000*qq[i])
    zscore_topseed_np <- randomWalk_sparse(intM=mutualknn50, queryCells=rownames(mutualknn50)[seed_p0.05_top], gamma=0.05)
    # only top seed times z factor
    zscore_topseed_np_scaled <- zscore_topseed_np %>% capOutlierQuantile(., 0.99) %>% max_min_scale
    zscoreWeighted2 <- zscoreWeighted
    zscoreWeighted2$zscore_topseed_np_scaled <- zscore_topseed_np_scaled
    sum(zscoreWeighted2$zscore_topseed_np_scaled==0) %>% print
    # zscoreWeighted2 <- zscoreWeighted2[zscoreWeighted2$zscore_topseed_np_scaled!=0, ]
    zscoreWeighted2$rank_np2 <- BiocGenerics::rank(-zscoreWeighted2$zscore_topseed_np_scaled, ties.method="first")

    quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank_np2, 2), table)
    mm[i] <- sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2))[1]
    quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank, 2), table)
    nn[i] <- sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2))[1]
    aa[[i]] <- table(zscoreWeighted$Celltype2[seed_p0.05_top])
    auc_v[i] <- auc(zscoreWeighted2$Celltype2=="Mono", zscoreWeighted2$zscore_topseed_np_scaled)
}
auc_v
[1] 0.999996 0.998352 0.998720 0.998876 0.993632 0.985316
 [7] 0.994288 0.965752 0.977104 0.948548 0.945716 0.905984
[13] 0.917532 0.928736 0.930376 0.912276 0.884284 0.865368
[19] 0.868364 0.858196 0.835808 0.826488 0.813260 0.808484
[25] 0.800428 0.806868 0.801884 0.800804 0.791468 0.796248
[31] 0.787364 0.788424 0.781360 0.775852 0.769056 0.771572
[37] 0.770072 0.767252 0.754256 0.743172 0.732644 0.735424
[43] 0.731288 0.729520 0.716512 0.707368 0.704404 0.706364
[49] 0.697484 0.691816

terminal_x <- 20
seed_mono <- sapply(aa, "[", 1) # number of mono in seed 
seed_nk <- sapply(aa, "[", 2) # number of nk in seed 
seed_nk[is.na(seed_nk)] <- 0

loop_category <- mapply(rep, 1:terminal_x, each=(seed_mono[1:terminal_x]+seed_nk[1:terminal_x])) %>% unlist
loop_num <- rep(1, sum((seed_mono[1:terminal_x]+seed_nk[1:terminal_x])))
loop_type <- list(); length(loop_type) <- terminal_x
for(i in 1:terminal_x){
    loop_type[[i]] <- c(rep("Mono", seed_mono[1:terminal_x][i]), rep("NK", seed_nk[1:terminal_x][i]))
}
loop_type <- unlist(loop_type) %>% factor(., levels=c("NK", "Mono"))

mydata <- data.frame(loop_category, loop_num, loop_type)
pp <- ggplot(mydata, aes(x = loop_category, y = loop_num, fill = loop_type))  + 
    geom_bar(position="fill", stat="identity") +scale_y_continuous(labels = scales::percent) +
    # scale_fill_viridis(discrete = T)
  scale_fill_manual(values = c("#94c954", "#168a53")) + pretty_plot()  + coord_cartesian(ylim=c(0.5, 1))
ggsave(pp, file = "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-12.pdf", width = 5, height = 3)

# ---- for tpr
terminal_x <- 20
df2 <- data.frame(accuracy=mm, percent=1:50)
df2 <- df2[1:terminal_x, ]
pp <- ggplot(data=df2, aes(x=percent, y=accuracy/100)) +
        geom_bar(stat="identity", fill="#2ea196") + pretty_plot() + coord_cartesian(ylim=c(0.5, 1))
pp <- pp + geom_hline(yintercept=.72, linetype="dashed", 
                color = "#ec4a9a", size=0.8)
ggsave(pp, file = "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-13.pdf", width = 5, height = 3)

# ---- for auc
terminal_x <- 20
df2 <- data.frame(auc=auc_v, percent=1:50)
df2 <- df2[1:terminal_x, ]
pp <- ggplot(data=df2, aes(x=percent, y=auc)) +
        geom_bar(stat="identity", fill="#2ea196") + pretty_plot() + coord_cartesian(ylim=c(0.5, 1))
pp <- pp + geom_hline(yintercept=0.799, linetype="dashed", 
                color = "#ec4a9a", size=0.8)
ggsave(pp, file = "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-14.pdf", width = 5, height = 3)



# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# test for k used for knn-graph construction
mutualknn30 <- getmutualknn(lsimat50, 30) %>% strong_connect
mutualknn30_graph <- graph_from_adjacency_matrix(mutualknn30, mode = "undirected", diag = F)

mutualknn5 <- getmutualknn(lsimat50, 5) %>% strong_connect
mutualknn5_graph <- graph_from_adjacency_matrix(mutualknn5, mode = "undirected", diag = F)

mutualknn10 <- getmutualknn(lsimat50, 10) %>% strong_connect
mutualknn10_graph <- graph_from_adjacency_matrix(mutualknn10, mode = "undirected", diag = F)

mutualknn20 <- getmutualknn(lsimat50, 20) %>% strong_connect
mutualknn20_graph <- graph_from_adjacency_matrix(mutualknn20, mode = "undirected", diag = F)

mutualknn50 <- getmutualknn(lsimat50, 50) %>% strong_connect
mutualknn50_graph <- graph_from_adjacency_matrix(mutualknn50, mode = "undirected", diag = F)

qq <- 0.05
seed_p0.05_top <- zscoreWeighted$rank<=1000*qq
zscore_topseed_np <- randomWalk_sparse(intM=mutualknn50, queryCells=(rownames(mutualknn50)[seed_p0.05_top], gamma=0.05)
zscore_topseed_np_scaled <- zscore_topseed_np %>% capOutlierQuantile(., 0.99) %>% max_min_scale

graph_list <- list(mutualknn5, mutualknn10, mutualknn20, mutualknn30, mutualknn50)
mm <- rep(0, 50)
nn <- rep(0, 50)
auc_v <-  rep(0, 5) # auc
auc_z <-  rep(0, 5)
aa <- list(); length(aa) <- 50
for (i in 1:length(graph_list)){
    seed_p0.05_top <- zscoreWeighted$rank<=1000*0.05
    zscore_topseed_np <- randomWalk_sparse(intM=graph_list[[i]], queryCells=rownames(graph_list[[i]])[seed_p0.05_top], gamma=0.05)
    # only top seed times z factor
    zscore_topseed_np_scaled <- zscore_topseed_np %>% capOutlierQuantile(., 0.99) %>% max_min_scale
    zscoreWeighted2 <- zscoreWeighted
    zscoreWeighted2$zscore_topseed_np_scaled <- zscore_topseed_np_scaled
    sum(zscoreWeighted2$zscore_topseed_np_scaled==0) %>% print
    zscoreWeighted2 <- zscoreWeighted2[zscoreWeighted2$zscore_topseed_np_scaled!=0, ]
    zscoreWeighted2$rank_np2 <- BiocGenerics::rank(-zscoreWeighted2$zscore_topseed_np_scaled, ties.method="first")
    message("np score")
    print(quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank_np2, 2), table))
    # mm[i] <- sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2))[1]
    message("raw score")
    print(quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank, 2), table))
    auc_v[i] <- auc(zscoreWeighted2$Celltype2=="Mono", zscoreWeighted2$zscore_topseed_np_scaled)
    auc_z[i] <- auc(zscoreWeighted2$Celltype2=="Mono", zscoreWeighted2$zscore)

    # nn[i] <- sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2))[1]
    # aa[[i]] <- table(zscoreWeighted$Celltype2[seed_p0.05_top])
}
# k=5 (474/500)=0.948; 72.2
np score
$`(0.002,500]`

Mono   NK 
 474   26 

$`(500,1e+03]`

Mono   NK 
  26  473 

raw score
$`(0.001,500]`

Mono   NK 
 361  138 

$`(500,1e+03]`

Mono   NK 
 139  361 

# k=10 0.992; (327/(327+129))=0.717
np score
np score
$`(0.002,500]`

Mono   NK 
 494    6 

$`(500,1e+03]`

Mono   NK 
   6  493 

raw score
$`(0.001,500]`

Mono   NK 
 361  139 

$`(500,1e+03]`

Mono   NK 
 139  360 

# k=20 0.984
np score
$`(0.002,500]`

Mono   NK 
 492    8 

$`(500,1e+03]`

Mono   NK 
   8  491 

raw score
$`(0.001,500]`

Mono   NK 
 361  139 

$`(500,1e+03]`

Mono   NK 
 139  360 

# k=30 0.97
np score
$`(0.002,500]`

Mono   NK 
 485   15 

$`(500,1e+03]`

Mono   NK 
  15  484 

raw score
$`(0.001,500]`

Mono   NK 
 361  139 

$`(500,1e+03]`

Mono   NK 
 139  360 

# k=50 =0.952
np score
$`(0.002,500]`

Mono   NK 
 476   24 

$`(500,1e+03]`

Mono   NK 
  24  475 

raw score
$`(0.001,500]`

Mono   NK 
 361  139 

$`(500,1e+03]`

Mono   NK 
 139  360 


library(ggplot2)
df <- data.frame(Numberofk=1:5, TPR=c(0.948, 0.992, 0.984, 0.97, 0.952)*100)
df2 <- data.frame(Numberofk=1:5, auc=c(auc_v)*100)
p<-ggplot(data=df, aes(x=Numberofk, y=TPR)) +
  geom_bar(stat="identity", fill="steelblue") + pretty_plot() + geom_text(aes(label=TPR), vjust=1.6, color="white", size=3.5)+ coord_flip() + coord_cartesian(ylim=c(50, 100))
ggsave(p, file = "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-15.pdf", width = 4, height = 3)
p<-ggplot(data=df2, aes(x=Numberofk, y=auc)) +
  geom_bar(stat="identity", fill="steelblue") + pretty_plot() + geom_text(aes(label=auc), vjust=1.6, color="white", size=3.5)+ coord_flip() + coord_cartesian(ylim=c(50, 100))
ggsave(p, file = "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-16.pdf", width = 4, height = 3)


setwd("/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20211227-bulksimulation/data")
# Import and run
# ATAC
peaksdf <- fread("/Users/fyu/Documents/GitHub/mpn-gwas/data/atac/26August2017_EJCsamples_allReads_250bp.bed") # unified peaks; 556270      3
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("/Users/fyu/Documents/GitHub/mpn-gwas/data/atac/26August2017_EJCsamples_allReads_250bp.counts.txt")) # count table for peaks across 18 cell types; 556270     18
mono_trait <- "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/finemappedtraits_hg19/mono.PP001.bed"
trait <- "mono"

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# test for n to various values including 500, 1000, 5000, 7500, 10000, and 20000 to test the effects of sequencing depth. 
read_number <- c(1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000, 100000)
auc_z <- auc_v <-  rep(0, 9) # auc
for (i in 1:length(read_number)){
    message(read_number[i])
    test_q40_n200_c9 <- simulate_scatac(n_cells=500, which_celltypes=c("Mono", "NK"), n_frags_per_cell = read_number[i], 
                                      rate_noise = 0.3, seed = 9527, shuffle = FALSE, bulk=counts)
    columns <- paste(rep(c("Mono", "NK"), each=500), rep(1:500, times=2), sep="_")
    print(length(columns))
    dimnames(test_q40_n200_c9)[[2]] = columns
    idx_no0 <- rowSums(test_q40_n200_c9)!=0; sum(!idx_no0)
            # Create objects for g-chromVAR
    SE_2 <- SummarizedExperiment(assays = list(counts = test_q40_n200_c9[idx_no0, ]),
                            rowData = peaks[idx_no0 ], 
                            colData = DataFrame(names = colnames(test_q40_n200_c9[idx_no0, ])))
    SE_2 <- addGCBias(SE_2, genome = BSgenome.Hsapiens.UCSC.hg19)
    mytrait <- importBedScore(rowRanges(SE_2), mono_trait, colidx = 5)
            # Run g-chromVAR
    bg <- getBackgroundPeaks(SE_2,niterations=200)
    dev <- computeWeightedDeviations(SE_2, mytrait, background_peaks = bg)
            # Reformat results
    zscoreWeighted <- reshape2::melt(t(assays(dev)[["z"]]))
    colnames(zscoreWeighted) <- c("Celltype","Trait","zscore")
    zscoreWeighted$Celltype2 <- stringr::str_split(zscoreWeighted$Celltype, "_", simplify=T)[, 1]
    zscoreWeighted$rank <- rank(-zscoreWeighted$zscore)
    plot(1:length(zscoreWeighted$zscore), zscoreWeighted$zscore)

    lsimat30=matSVD(mat=test_q40_n200_c9[idx_no0, ], nComponents = 30, binarize = TRUE) # 1000    50
    mutualknn50 <- getmutualknn(lsimat30, 30) %>% strong_connect
    mutualknn50_graph <- graph_from_adjacency_matrix(mutualknn50, mode = "undirected", diag = F)

    seed_p0.05_top <- zscoreWeighted$rank<=1000*0.05
    zscore_topseed_np <- randomWalk_sparse(intM=mutualknn50, queryCells=rownames(mutualknn50)[seed_p0.05_top], gamma=0.05)
    # only top seed times z factor
    zscore_topseed_np_scaled <- zscore_topseed_np %>% capOutlierQuantile(., 0.99) %>% max_min_scale
    zscoreWeighted2 <- zscoreWeighted
    zscoreWeighted2$zscore_topseed_np_scaled <- zscore_topseed_np_scaled
    sum(zscoreWeighted2$zscore_topseed_np_scaled==0) %>% print
    # zscoreWeighted2 <- zscoreWeighted2[zscoreWeighted2$zscore_topseed_np_scaled!=0, ]
    zscoreWeighted2$rank_np2 <- BiocGenerics::rank(-zscoreWeighted2$zscore_topseed_np_scaled, ties.method="first")
    message("np score")
    print(quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank_np2, 2), table))
    # mm[i] <- sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2))[1]
    message("raw score")
    print(quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank, 2), table))
    # nn[i] <- sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2))[1]
    # aa[[i]] <- table(zscoreWeighted$Celltype2[seed_p0.05_top])
    auc_v[i] <- auc(zscoreWeighted2$Celltype2=="Mono", zscoreWeighted2$zscore_topseed_np_scaled) 
    auc_z[i] <- auc(zscoreWeighted$Celltype2=="Mono", zscoreWeighted$zscore) 
}

> auc_v
[1] 0.626386 0.722396 0.791452 0.908560 0.993608 0.999600
[7] 0.999072 0.999924 0.999948
> auc_z
[1] 0.627852 0.666536 0.712436 0.768184 0.799432 0.905876
[7] 0.974160 0.994868 0.998824
auc_v <- c(0.59, 0.70, 0.78, 0.969, 0.977, 0.978, 0.998, 0.996, 1)
auc_z <- c(0.63, 0.67, 0.71, 0.791, 0.801, 0.906, 0.974, 0.994, 0.998)
1000
here0
here1
[1] 1000
Binarizing matrix...
Computing Term Frequency IDF...
Computing SVD using irlba...
[info] fast knn
[1] 37
Stationary step: 177
Stationary Delta: 9.56576952981579e-06
[1] 1
np score
$`(0.001,500]`

Mono   NK 
 274  226 

$`(500,1e+03]`

Mono   NK 
 226  274 

raw score
$`(0.001,500]`

Mono   NK 
 299  201 

$`(500,1e+03]`

Mono   NK 
 201  299 

2500
here0
here1
[1] 1000
Binarizing matrix...
Computing Term Frequency IDF...
Computing SVD using irlba...
[info] fast knn
[1] 24
Stationary step: 166
Stationary Delta: 9.98494098712555e-06
[1] 1
np score
$`(0.001,500]`

Mono   NK 
 303  197 

$`(500,1e+03]`

Mono   NK 
 197  303 

raw score
$`(0.001,500]`

Mono   NK 
 304  196 

$`(500,1e+03]`

Mono   NK 
 196  304 

5000
here0
here1
[1] 1000
Binarizing matrix...
Computing Term Frequency IDF...
Computing SVD using irlba...
[info] fast knn
[1] 15
Stationary step: 157
Stationary Delta: 9.5012592103299e-06
[1] 1
np score
$`(0.001,500]`

Mono   NK 
 328  172 

$`(500,1e+03]`

Mono   NK 
 172  328 

raw score
$`(0.001,500]`

Mono   NK 
 329  171 

$`(500,1e+03]`

Mono   NK 
 171  329 

7500
here0
here1
[1] 1000
Binarizing matrix...
Computing Term Frequency IDF...
Computing SVD using irlba...
[info] fast knn
[1] 15
Stationary step: 151
Stationary Delta: 9.84549352645718e-06
[1] 2
np score
$`(0.001,500]`

Mono   NK 
 381  119 

$`(500,1e+03]`

Mono   NK 
 119  381 

raw score
$`(0.001,500]`

Mono   NK 
 344  156 

$`(500,1e+03]`

Mono   NK 
 156  344 

10000
here0
here1
[1] 1000
Binarizing matrix...
Computing Term Frequency IDF...
Computing SVD using irlba...
[info] fast knn
[1] 9
Stationary step: 127
Stationary Delta: 9.47141269831945e-06
[1] 1
np score
$`(0.001,500]`

Mono   NK 
 486   14 

$`(500,1e+03]`

Mono   NK 
  14  486 

raw score
$`(0.001,500]`

Mono   NK 
 361  139 

$`(500,1e+03]`

Mono   NK 
 139  361 

25000
here0
here1
[1] 1000
Binarizing matrix...
Computing Term Frequency IDF...
Computing SVD using irlba...
[info] fast knn
[1] 9
Stationary step: 79
Stationary Delta: 9.40132994156217e-06
[1] 1
np score
$`(0.001,500]`

Mono   NK 
 499    1 

$`(500,1e+03]`

Mono   NK 
   1  499 

raw score
$`(0.001,500]`

Mono   NK 
 409   91 

$`(500,1e+03]`

Mono   NK 
  91  409 

50000
here0
here1
[1] 1000
Binarizing matrix...
Computing Term Frequency IDF...
Computing SVD using irlba...
[info] fast knn
[1] 4
Stationary step: 49
Stationary Delta: 9.91747503269439e-06
[1] 1
np score
$`(0.001,500]`

Mono   NK 
 497    3 

$`(500,1e+03]`

Mono   NK 
   3  497 

raw score
$`(0.001,500]`

Mono   NK 
 463   37 

$`(500,1e+03]`

Mono   NK 
  37  463 

75000
here0
here1
[1] 1000
Binarizing matrix...
Computing Term Frequency IDF...
Computing SVD using irlba...
[info] fast knn
[1] 2
Stationary step: 66
Stationary Delta: 9.54365057893913e-06
[1] 1
np score
$`(0.001,500]`

Mono   NK 
 499    1 

$`(500,1e+03]`

Mono   NK 
   1  499 

raw score
$`(0.001,500]`

Mono   NK 
 480   20 

$`(500,1e+03]`

Mono   NK 
  20  480 

1e+05
here0
here1
[1] 1000
Binarizing matrix...
Computing Term Frequency IDF...
Computing SVD using irlba...
[info] fast knn
[1] 3
Stationary step: 64
Stationary Delta: 9.54369822855601e-06
[1] 1
np score
$`(0.001,500]`

Mono   NK 
 499    1 

$`(500,1e+03]`

Mono   NK 
   1  499 

raw score
$`(0.001,500]`

Mono   NK 
 493    7 

$`(500,1e+03]`

Mono   NK 
   7  493 

# c(2500, 2500, 5000, 5000, 7500, 7500, 10000, 10000, 25000, 25000, 50000, 50000)
> auc_v
[1] 0.626386 0.722396 0.791452 0.908560 0.993608 0.999600
[7] 0.999072 0.999924 0.999948
> auc_z
[1] 0.627852 0.666536 0.712436 0.768184 0.799432 0.905876
[7] 0.974160 0.994868 0.998824
auc_v <- c(0., 0.722, 0.78, 0.969, 0.977, 0.978, 0.998, 0.996, 1)
auc_z <- c(0.63, 0.67, 0.71, 0.791, 0.801, 0.906, 0.974, 0.994, 0.998)
TPR=(matrix(c(c(0.658, 0.688, 0.722, 0.818, 0.926), c(0.656, 0.762, 0.972, 0.998, 0.994)), ncol=2) %>% t %>% c)
auROC=(matrix(c(c(0.712436, 0.768184, 0.799432, 0.905876, 0.974160), c(0.791452, 0.908560, 0.993608, 0.999600, 0.999072)), ncol=2) %>% t %>% c)
df3 <- data.frame(reads=c("5000", "5000", "7500", "7500", "10000", "10000", "25000", "25000", "50000", "50000"), 
                  TPR=TPR*100,
                  auROC=auROC*100,
                  type=c("Z score", "SCAVENGE", "Z score", "SCAVENGE", "Z score", "SCAVENGE", "Z score", "SCAVENGE", "Z score", "SCAVENGE"))
df3$type <- factor(df3$type, levels=c("Z score", "SCAVENGE") )
df3$reads <- factor(df3$reads, levels=c("5000", "7500", "10000", "25000", "50000") )

p <- ggplot(df3, aes(x=reads, y=TPR, fill=type)) + 
   geom_bar(stat="identity", position=position_dodge())
  
p1 <- p + scale_fill_brewer(palette="Paired") + pretty_plot()
p2 <- p + scale_fill_brewer(palette="Reds") + pretty_plot()
p3 <- p + scale_fill_brewer(palette="Greens") + pretty_plot()
p <- ggplot(df3, aes(x=reads, y=TPR, fill=type)) + 
   geom_bar(stat="identity", position=position_dodge())+ coord_cartesian(ylim=c(50, 100)) + 
   scale_fill_brewer(palette="Greens") + pretty_plot()
ggsave(p, file = "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-17.pdf", width = 4, height = 3)
p <- ggplot(df3, aes(x=reads, y=auROC, fill=type)) + 
   geom_bar(stat="identity", position=position_dodge())+ coord_cartesian(ylim=c(50, 100)) + 
   scale_fill_brewer(palette="Greens") + pretty_plot()
ggsave(p, file = "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-18.pdf", width = 4, height = 3)


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# We set q to various values including 0.1, 0.3, 0.5 to test the robustness to noise
auc_z <- auc_v <-  rep(0, 4) # auc
q_number <- c(0.25, 0.30, 0.35, 0.4)
for (i in 1:length(q_number)){
    message(q_number[i])
    test_q40_n200_c9 <- simulate_scatac(n_cells=500, which_celltypes=c("Mono", "NK"), n_frags_per_cell = 10000, 
                                      rate_noise = q_number[i], seed = 9527, shuffle = FALSE, bulk=counts)
    columns <- paste(rep(c("Mono", "NK"), each=500), rep(1:500, times=2), sep="_")
    print(length(columns))
    dimnames(test_q40_n200_c9)[[2]] = columns
    idx_no0 <- rowSums(test_q40_n200_c9)!=0; sum(!idx_no0)
            # Create objects for g-chromVAR
    SE_2 <- SummarizedExperiment(assays = list(counts = test_q40_n200_c9[idx_no0, ]),
                            rowData = peaks[idx_no0 ], 
                            colData = DataFrame(names = colnames(test_q40_n200_c9[idx_no0, ])))
    SE_2 <- addGCBias(SE_2, genome = BSgenome.Hsapiens.UCSC.hg19)
    mytrait <- importBedScore(rowRanges(SE_2), mono_trait, colidx = 5)
            # Run g-chromVAR
    bg <- getBackgroundPeaks(SE_2,niterations=200)
    dev <- computeWeightedDeviations(SE_2, mytrait, background_peaks = bg)
            # Reformat results
    zscoreWeighted <- reshape2::melt(t(assays(dev)[["z"]]))
    colnames(zscoreWeighted) <- c("Celltype","Trait","zscore")
    zscoreWeighted$Celltype2 <- stringr::str_split(zscoreWeighted$Celltype, "_", simplify=T)[, 1]
    zscoreWeighted$rank <- rank(-zscoreWeighted$zscore)
    plot(1:length(zscoreWeighted$zscore), zscoreWeighted$zscore)

    lsimat30=matSVD(mat=test_q40_n200_c9[idx_no0, ], nComponents = 30, binarize = TRUE) # 1000    50
    mutualknn50 <- getmutualknn(lsimat30, 30) %>% strong_connect
    mutualknn50_graph <- graph_from_adjacency_matrix(mutualknn50, mode = "undirected", diag = F)

    seed_p0.05_top <- zscoreWeighted$rank<=1000*0.05
    zscore_topseed_np <- randomWalk_sparse(intM=mutualknn50, queryCells=rownames(mutualknn50)[seed_p0.05_top], gamma=0.05)
    # only top seed times z factor
    zscore_topseed_np_scaled <- zscore_topseed_np %>% capOutlierQuantile(., 0.99) %>% max_min_scale
    zscoreWeighted2 <- zscoreWeighted
    zscoreWeighted2$zscore_topseed_np_scaled <- zscore_topseed_np_scaled
    sum(zscoreWeighted2$zscore_topseed_np_scaled==0) %>% print
    # zscoreWeighted2 <- zscoreWeighted2[zscoreWeighted2$zscore_topseed_np_scaled!=0, ]
    zscoreWeighted2$rank_np2 <- BiocGenerics::rank(-zscoreWeighted2$zscore_topseed_np_scaled, ties.method="first")
    message("np score")
    print(quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank_np2, 2), table))
    # mm[i] <- sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2))[1]
    message("raw score")
    print(quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank, 2), table))
    # nn[i] <- sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2))[1]
    # aa[[i]] <- table(zscoreWeighted$Celltype2[seed_p0.05_top])
    auc_v[i] <- auc(zscoreWeighted2$Celltype2=="Mono", zscoreWeighted2$zscore_topseed_np_scaled) 
    auc_z[i] <- auc(zscoreWeighted2$Celltype2=="Mono", zscoreWeighted2$zscore) 
}
0.25
here0
here1
[1] 1000
Binarizing matrix...
Computing Term Frequency IDF...
Computing SVD using irlba...
[info] fast knn
[1] 18
Stationary step: 139
Stationary Delta: 9.81179470304452e-06
[1] 1
np score
$`(0.001,500]`

Mono   NK 
 499    1 

$`(500,1e+03]`

Mono   NK 
   1  499 

raw score
$`(0.001,500]`

Mono   NK 
 363  137 

$`(500,1e+03]`

Mono   NK 
 137  363 

0.3
here0
here1
[1] 1000
Binarizing matrix...
Computing Term Frequency IDF...
Computing SVD using irlba...
[info] fast knn
[1] 9
Stationary step: 127
Stationary Delta: 9.47141269831945e-06
[1] 1
np score
$`(0.001,500]`

Mono   NK 
 486   14 

$`(500,1e+03]`

Mono   NK 
  14  486 

raw score
$`(0.001,500]`

Mono   NK 
 361  139 

$`(500,1e+03]`

Mono   NK 
 139  361 

0.35
here0
here1
[1] 1000
Binarizing matrix...
Computing Term Frequency IDF...
Computing SVD using irlba...
[info] fast knn
[1] 7
Stationary step: 84
Stationary Delta: 9.98085703499836e-06
[1] 1
np score
$`(0.001,500]`

Mono   NK 
 495    5 

$`(500,1e+03]`

Mono   NK 
   5  495 

raw score
$`(0.001,500]`

Mono   NK 
 352  148 

$`(500,1e+03]`

Mono   NK 
 148  352 

0.4
here0
here1
[1] 1000
Binarizing matrix...
Computing Term Frequency IDF...
Computing SVD using irlba...
[info] fast knn
[1] 10
Stationary step: 129
Stationary Delta: 9.72483358510218e-06
[1] 2
np score
$`(0.001,500]`

Mono   NK 
 495    5 

$`(500,1e+03]`

Mono   NK 
   5  495 

raw score
$`(0.001,500]`

Mono   NK 
 348  152 

$`(500,1e+03]`

Mono   NK 
 152  348 
> auc_v
[1] 0.999996 0.993608 0.996492 0.998352
> auc_z
[1] 0.814540 0.799432 0.783644 0.773496
auROC=(matrix(c(auc_z, auc_v), ncol=2) %>% t %>% c)
df3 <- data.frame(noise=c(0.25, 0.25, 0.3, 0.3, 0.35, 0.35, 0.4, 0.4)*100, 
                  TPR=c(0.726, 0.998, 0.722, 0.972, 0.705, 0.99, 0.696, 0.99)*100,
                  auROC=auROC,
                  type=c("Z score", "SCAVENGE", "Z score", "SCAVENGE", "Z score", "SCAVENGE", "Z score", "SCAVENGE"))
df3$type <- factor(df3$type, levels=c("Z score", "SCAVENGE") )
p <- ggplot(df3, aes(x=noise, y=TPR, fill=type)) + 
   geom_bar(stat="identity", position=position_dodge()) 
  
p1 <- p + scale_fill_brewer(palette="Paired") + pretty_plot()
p2 <- p + scale_fill_brewer(palette="Reds") + pretty_plot()
p <- ggplot(df3, aes(x=noise, y=TPR, fill=type)) + 
   geom_bar(stat="identity", position=position_dodge())+ #coord_cartesian(ylim=c(50, 100)) + 
   scale_fill_brewer(palette="Reds") + pretty_plot()
ggsave(p, file = "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-19.pdf", width = 4, height = 3)
p <- ggplot(df3, aes(x=noise, y=auROC, fill=type)) + 
   geom_bar(stat="identity", position=position_dodge())+ #coord_cartesian(ylim=c(50, 100)) + 
   scale_fill_brewer(palette="Reds") + pretty_plot()
ggsave(p, file = "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/remade-20.pdf", width = 4, height = 3)


# ------
# sparsity of simulated data
# ------
test_q0_n200_c2 <- simulate_scatac(n_cells=500, which_celltypes=c("Mono", "NK"), n_frags_per_cell = 10000, 
                                    rate_noise = 0.3, seed = 9527, shuffle = FALSE, bulk=counts)
peakbycellmat <- test_q0_n200_c2 # [1] 91200  4628
peakbycellmat_bin <- peakbycellmat
peakbycellmat_bin@x[peakbycellmat_bin@x>=1] <- 1
peak_detected_ratio_C2 <- data.frame(rowSums(peakbycellmat_bin)*100/ncol(peakbycellmat_bin), "C2")
cell_coverage_ratio_C2 <- data.frame(colSums(peakbycellmat_bin)*100/nrow(peakbycellmat_bin), "C2")

peak_detected_ratio_C2[, 1] = 100 - peak_detected_ratio_C2[, 1]
cell_coverage_ratio_C2[, 1] = 100 - cell_coverage_ratio_C2[, 1]
colnames(peak_detected_ratio_C2) <- colnames(cell_coverage_ratio_C2) <- c("percent", "condition")
# library(ggpubr)

# p <- ggplot(peak_detected_ratio, aes(x=percent, fill=condition)) +
#         xlab("Sparsity of peaks (%)") +
#         geom_density(alpha=0.5) + xlim(75, 100) #  + theme_pubr()
# p = p + theme_pubclean(flip=T)
# ggsave(p, file = "percentage_of_peak_sparsity_simulation.pdf", width = 5, height = 5)


# p <- ggplot(cell_coverage_ratio, aes(x=percent, fill=condition)) +
#         xlab("Sparsity of cells (%)") +
#         geom_density(alpha=0.5) + xlim(75, 100) #  + theme_pubr()
# p = p + theme_pubclean(flip=T)
# ggsave(p, file = "percentage_of_cell_sparsity_simulation.pdf", width = 5, height = 5)


## c9
c9 <- c("Mono", "mDC", "GMP-C", "CLP", "B", "pDC", "CD8", "CD4", "NK")
test_q30_n200_c9 <- simulate_scatac(n_cells=200, which_celltypes=c9, n_frags_per_cell = 10000, 
                            rate_noise = 0.3, seed = 9527, shuffle = FALSE, bulk=counts)
columns <- paste(rep(c9, each=200), rep(1:200, times=9), sep="_")
print(length(columns))
dimnames(test_q40_n200_c9)[[2]] = columns
idx_no0 <- rowSums(test_q40_n200_c9)!=0; sum(!idx_no0)
        # Create objects for g-chromVAR
test_q40_n200_c9 = test_q40_n200_c9[idx_no0, ]
peakbycellmat <- test_q40_n200_c9 # [1] 91200  4628
peakbycellmat_bin <- peakbycellmat
peakbycellmat_bin@x[peakbycellmat_bin@x>=1] <- 1
peak_detected_ratio_C9 <- data.frame(rowSums(peakbycellmat_bin)*100/ncol(peakbycellmat_bin), "C9")
cell_coverage_ratio_C9 <- data.frame(colSums(peakbycellmat_bin)*100/nrow(peakbycellmat_bin), "C9")

peak_detected_ratio_C9[, 1] = 100 - peak_detected_ratio_C9[, 1]
cell_coverage_ratio_C9[, 1] = 100 - cell_coverage_ratio_C9[, 1]
colnames(peak_detected_ratio_C9) <- colnames(cell_coverage_ratio_C9) <- c("percent", "condition")
library(ggpubr)

peak_detected_ratio <- rbind.data.frame(peak_detected_ratio_C2, peak_detected_ratio_C9)
cell_coverage_ratio <- rbind.data.frame(cell_coverage_ratio_C2, cell_coverage_ratio_C9)

p <- ggplot(peak_detected_ratio, aes(x=percent, fill=condition)) +
        xlab("Sparsity of peaks (%)") +
        geom_density(alpha=0.3) + xlim(75, 100) #  + theme_pubr()
p = p + theme_pubclean(flip=T)
ggsave(p, file = "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/percentage_of_peak_sparsity_simulation_.pdf", width = 5, height = 5)


p <- ggplot(cell_coverage_ratio, aes(x=percent, fill=condition)) +
        xlab("Sparsity of cells (%)") +
        geom_density(alpha=0.3) + xlim(75, 100) #  + theme_pubr()
p = p + theme_pubclean(flip=T)
ggsave(p, file = "/Users/fyu/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Vijay/project/0817-singleVar/data/20220220-nbt_revision-1/percentage_of_cell_sparsity_simulation.pdf", width = 5, height = 5)


####----------------------------------
# semi simulation with 10x scdata




####----------------------------------
# real 10x scata data 
# we have set this dataset as the vignette of SCAVENGE R package, please go through it to get the result
trait_file <- "finemappedtraits_hg19/mono.PP001.bed"
pbmc5krda <- paste0("pbmc5k_SE.rda")
load(pbmc5krda)

### gchromVAR analysis
SE_pbmc5k <- addGCBias(SE_pbmc5k, genome = BSgenome.Hsapiens.UCSC.hg19)
SE_pbmc5k_bg <- getBackgroundPeaks(SE_pbmc5k, niterations=200)
trait_import <- importBedScore(rowRanges(SE_pbmc5k), trait_file, colidx=5)
SE_pbmc5k_DEV <- computeWeightedDeviations(SE_pbmc5k, trait_import, background_peaks = SE_pbmc5k_bg)

z_score_mat <- data.frame(colData(SE_pbmc5k), z_score=t(assays(SE_pbmc5k_DEV)[["z"]]) %>% c)
### Generate the seed cell index (using the top 5% if too many cells are eligible)
seed_idx <- seedindex(z_score_mat$z_score, 0.05)
scale_factor <- cal_scalefactor(z_score=z_score_mat$z_score, 0.01)

# **Calculate tfidf-mat**
peak_by_cell_mat <- assay(SE_pbmc5k)
tfidf_mat <- tfidf(bmat=peak_by_cell_mat, mat_binary=TRUE, TF=TRUE, log_TF=TRUE)
# **Calculate lsi-mat**
lsi_mat <- do_lsi(tfidf_mat, dims=30)
# **Calculate m-knn graph**
mutualknn50 <- getmutualknn(lsi_mat, 50)

### Network propagation  
np_score <- randomWalk_sparse(intM=mutualknn50, queryCells=rownames(mutualknn50)[seed_p0.05_top], gamma=0.05)
# **Trait relevant score (TRS) with scaled and normalized**  
A few cells are singletons are removed from further analysis
omit_idx <- np_score==0
sum(omit_idx)
mutualknn50 <- mutualknn50[!omit_idx, !omit_idx]
np_score <- np_score[!omit_idx]
TRS <- np_score %>% capOutlierQuantile(., 0.99) %>% max_min_scale
TRS <- TRS * scale_factor
mono_mat <- data.frame(z_score_mat[!omit_idx, ], seed_idx[!omit_idx], np_score, TRS)

### UMAP plots of cell type annotation and cell-to-cell graph 
# **Cell type annotation**
p <- ggplot(data=mono_mat, aes(x, y, color=color)) + geom_point(size=1, na.rm = TRUE) + 
    pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
p

# **Visualize cell-to-cell graph if you have low-dimensional coordinates such as UMAP1 and UMAP2**
mutualknn50_graph <- graph_from_adjacency_matrix(mutualknn50, mode = "undirected", diag = F)
plot.igraph(mutualknn50_graph, vertex.size=0.8, vertex.label=NA, vertex.color=adjustcolor("#c7ce3d", alpha.f = 1), vertex.frame.color=NA, 
            edge.color=adjustcolor("#443dce", alpha.f = 1), edge.width=0.3, edge.curved=.5, 
            layout=as.matrix(data.frame(mono_mat$x, mono_mat$y)))

### Comparsion before and after SCAVENGE analysis  
# - Z score based visualization   
# **Scatter plot**
viridis = c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF", "#5DC863FF", "#AADC32FF", "#FDE725FF")
p1 <- ggplot(data=mono_mat, aes(x, y, color=z_score)) + geom_point(size=1, na.rm = TRUE, alpha = 0.6) + 
scale_color_gradientn(colors = viridis) + scale_alpha()+
    pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
p1
# **Bar plot**  
pp1 <- ggplot(data=mono_mat,  aes(x=color, y=z_score))  +
    geom_boxplot(aes(fill=color, color=color), outlier.shape=NA) + 
    guides(fill=FALSE) + pretty_plot(fontsize = 10) +
    stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + theme(legend.position = "none")
pp1

# - SCAVENGE TRS based visualization  
# **Scatter plot**  
p2 <- ggplot(data=mono_mat, aes(x, y, color=TRS)) + geom_point(size=1, na.rm = TRUE, alpha = 0.6) + 
scale_color_gradientn(colors = viridis) + scale_alpha()+
    pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
p2
# **Bar plot**
pp2 <- ggplot(data=mono_mat,  aes(x=color, y=TRS))  +
    geom_boxplot(aes(fill=color, color=color), outlier.shape=NA) + 
    guides(fill=FALSE) + pretty_plot(fontsize = 10) +
    stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + theme(legend.position = "none")
pp2

write.csv(mono_mat, "scavenge-mono-pbmc5k.csv", row.names=F)

# please request more data (ArchR obj) for gene score of canonical markers for cell types
