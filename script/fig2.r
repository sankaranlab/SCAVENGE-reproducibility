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
setwd("/Users/fyu/Documents/GitHub/SCAVENGE-reproducibility/data")


# this function is adapted from https://github.com/pinellolab/scATAC-benchmarking/tree/master/Synthetic_Data
simulate_scatac <- function(n_cells, which_celltypes, n_frags_per_cell = 1000, 
                            rate_noise = 0, seed = 100, shuffle = FALSE, bulk=xx, peaks=peaks){
  # Reproducibility
  set.seed(seed)
#   which_celltypes <- sort(which_celltypes)
  stopifnot(rate_noise < 1) 
  stopifnot(n_frags_per_cell 100)
  n_peaks <- dim(bulk)[1]
  #--
  # Set up cell labels
  #--
  message("here0")
  if(length(n_cells) 1){
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
  # Generate cell-type specific peaks
  lapply(which_celltypes, function(celltype){
    
    # Apply different rates per cell depending on group label for generating cell-type specific peaks
    n_cells_this_celltype <- sum(cell_labels == celltype)
    counts_celltype <- bulk[, celltype]
    
    # Define probabilities
    #                        Prob observting frag                Total number of fragments epxpected; the 0.5s are for two alleles that will be simulated/added later
    prob_per_peaks <- counts_celltype/sum(counts_celltype) * (n_frags_per_cell*0.5 * (1-rate_noise)) + ((rate_noise*n_frags_per_cell)/n_peaks*0.5) 
    
    # Cap probabilities at something sensible
    prob_per_peaks <- ifelse(prob_per_peaks 0.9, 0.9, prob_per_peaks)
    
    # Represent the two haplotypes as two random draws
    mat1 <- (matrix(rbinom(n_peaks*n_cells_this_celltype, size = 1, prob = prob_per_peaks),
                    ncol = n_cells_this_celltype, byrow = FALSE) )
    mat2 <- (matrix(rbinom(n_peaks*n_cells_this_celltype, size = 1, prob = prob_per_peaks),
                    ncol = n_cells_this_celltype, byrow = FALSE) )
    
    mat <- mat1 + mat2
    Matrix(mat)
  }) %>% do.call(what = "cbind") -sparse_matrix
  return(sparse_matrix)
  message("here2")
  colnames(sparse_matrix) <- final_names
  message("here2")
  peaknames = paste(peaks$V1,peaks$V2,peaks$V3,sep = "_")
  rownames(sparse_matrix) <- peaknames
  message("here2")
  sparse_matrix
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
idx_no0 <- rowSums(test_q40_n200_c9)!=0; sum(!idx_no0) # remove peaks with no signals
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
colnames(zscoreWeighted) <- c("Celltype", "Trait", "zscore")
zscoreWeighted$Celltype2 <- stringr::str_split(zscoreWeighted$Celltype, "_", simplify=T)[, 1]
zscoreWeighted$rank <- rank(-zscoreWeighted$zscore)
# lsi mat
lsimat30 <- test_q40_n200_c9[idx_no0, ] %>% tfidf %>% do_lsi
# mknn graph with different k
mutualknn30 <- getmutualknn(lsimat30, 30)
mutualknn30_graph <- graph_from_adjacency_matrix(mutualknn30, mode = "undirected", diag = F)
mutualknn5 <- getmutualknn(lsimat30, 5)
mutualknn5_graph <- graph_from_adjacency_matrix(mutualknn5, mode = "undirected", diag = F)
mutualknn10 <- getmutualknn(lsimat30, 10)
mutualknn10_graph <- graph_from_adjacency_matrix(mutualknn10, mode = "undirected", diag = F)
mutualknn20 <- getmutualknn(lsimat30, 20)
mutualknn20_graph <- graph_from_adjacency_matrix(mutualknn20, mode = "undirected", diag = F)
mutualknn50 <- getmutualknn(lsimat30, 50)
mutualknn50_graph <- graph_from_adjacency_matrix(mutualknn50, mode = "undirected", diag = F)


seedindex()
qq <- 0.05
seed_p0.05_top <- seedindex(zscoreWeighted$zscore, qq)
zscore_topseed_np <- randomWalk_sparse(intM=mutualknn30, queryGenes=rownames(mutualknn30)[seed_p0.05_top], gamma=0.05)
# only top seed times z factor
zscore_topseed_np_scaled <- zscore_topseed_np %>% capOutlierQuantile(., 0.99) %>% max_min_scale
zscoreWeighted2 <- zscoreWeighted
zscoreWeighted2$zscore_topseed_np_scaled <- zscore_topseed_np_scaled
sum(zscoreWeighted2$zscore_topseed_np_scaled==0) # remove singletons isolated in the graph
zscoreWeighted2 <- zscoreWeighted2[zscoreWeighted2$zscore_topseed_np_scaled!=0, ]
zscoreWeighted2$rank_np2 <- BiocGenerics::rank(-zscoreWeighted2$zscore_topseed_np_scaled, ties.method="first")

# comparsion of ranks 
# TRS 
pp <- ggplot(zscoreWeighted2, aes(x = rank_np2, y = 0.3, fill = Celltype2)) +
  geom_col() + pretty_plot() + 
  scale_fill_manual(values = c("#ec4a9a", "#dddddd")) +
  xlab("rank")
pp

# z score
pp <- ggplot(zscoreWeighted2, aes(x = rank, y = 0.3, fill = Celltype2)) +
  geom_col() + pretty_plot() + 
  scale_fill_manual(values = c("#ec4a9a", "#dddddd")) +
  xlab("rank")
pp

# look at different quarters
# donutplots
quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank_np2, 4), table)
$`(0.061,236]`

Mono   NK 
 231    4 

$`(236,470]`

Mono   NK 
 211   24 

$`(470,705]`

Mono   NK 
   1  234 

$`(705,941]`

Mono   NK 
   7  228

96.8
# Hole size
hsize <- 2
(df1 <- data.frame(value=c(231, 4), group=c("Mono", "NK"), hsize))
df2 <- data.frame(value=c(211, 24), group=c("Mono", "NK"), hsize)
df3 <- data.frame(value=c(1, 234), group=c("Mono", "NK"), hsize)
df4 <- data.frame(value=c(7, 228), group=c("Mono", "NK"), hsize)
pp <- ggplot(df1, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
pp <- ggplot(df2, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
pp <- ggplot(df3, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
pp <- ggplot(df4, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())


(quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank, 4), table))
$`(0.001,251]`

Mono   NK 
 185   46 

$`(251,500]`

Mono   NK 
 142   89 

$`(500,750]`

Mono   NK 
  85  149 

$`(750,1e+03]`

Mono   NK 
  38  206 
72.6
(df1 <- data.frame(value=c(186, 46), group=c("Mono", "NK"), hsize))
df2 <- data.frame(value=c(142, 89), group=c("Mono", "NK"), hsize)
df3 <- data.frame(value=c(85, 149), group=c("Mono", "NK"), hsize)
df4 <- data.frame(value=c(38, 206), group=c("Mono", "NK"), hsize)
pp <- ggplot(df1, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
pp <- ggplot(df2, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
pp <- ggplot(df3, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
pp <- ggplot(df4, aes(x = hsize, y = value, fill = group)) +
  geom_col(color = "black") + geom_text(aes(label = value), position = position_stack(vjust = 0.5)) + coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ec4a9a", "#e0ecf4"))+ xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# investigation of different proportion of seed cells

qq <- (1:50)/100
mm <- rep(0, 50)
nn <- rep(0, 50)
aa <- list(); length(aa) <- 50
for (i in 1:length(qq)){
    seed_p0.05_top <- zscoreWeighted$rank<=1000*qq[i]
    zscore_topseed_np <- randomWalk_sparse(intM=mutualknn30, queryGenes=(1:nrow(mutualknn30))[seed_p0.05_top], gamma=0.05)
    # only top seed times z factor
    zscore_topseed_np_scaled <- zscore_topseed_np %>% capOutlierQuantile(., 0.99) %>% max_min_scale
    zscoreWeighted2 <- zscoreWeighted
    zscoreWeighted2$zscore_topseed_np_scaled <- zscore_topseed_np_scaled
    sum(zscoreWeighted2$zscore_topseed_np_scaled==0) %>% print
    zscoreWeighted2 <- zscoreWeighted2[zscoreWeighted2$zscore_topseed_np_scaled!=0, ]
    zscoreWeighted2$rank_np2 <- BiocGenerics::rank(-zscoreWeighted2$zscore_topseed_np_scaled, ties.method="first")

    quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank_np2, 2), table)
    mm[i] <- sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2))[1]
    quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank, 2), table)
    nn[i] <- sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2))[1]
    aa[[i]] <- table(zscoreWeighted$Celltype2[seed_p0.05_top])
}

terminal_x <- 30
seed_mono <- sapply(aa, "[", 1)
seed_nk <- sapply(aa, "[", 2)
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
  scale_fill_manual(values = c("#94c954", "#168a53")) + pretty_plot() 
ggsave(pp, file = "./fig_ok/seedpercent-mononk-30.pdf", width = 5, height = 3)

# ---- for accuracy
terminal_x <- 30
df2 <- data.frame(accuracy=mm, percent=1:50)
df2 <- df2[1:terminal_x, ]
pp <- ggplot(data=df2, aes(x=percent, y=accuracy)) +
        geom_bar(stat="identity", fill="#2ea196") + pretty_plot()
pp <- pp + geom_hline(yintercept=72, linetype="dashed", 
                color = "#ec4a9a", size=0.8)
ggsave(pp, file = "./fig_ok/seedpercent-accuracy-30.pdf", width = 5, height = 3)


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank_np2, 8), table)
sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2) )
quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank, 8), table)
sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2) )
qq <- 0.03
seed_p0.05_top <- zscoreWeighted$rank<=1000*qq
zscore_topseed_np <- randomWalk_sparse(intM=mutualknn30, queryGenes=(1:nrow(mutualknn30))[seed_p0.05_top], gamma=0.05)
# only top seed times z factor
zscore_topseed_np_scaled <- zscore_topseed_np %>% capOutlierQuantile(., 0.99) %>% max_min_scale
plot(1:1000, zscore_topseed_np_scaled)
sum(zscore_topseed_np[1:500]==0)
sum(zscore_topseed_np[501:1000]==0)

zscoreWeighted2 <- zscoreWeighted
zscoreWeighted2$zscore_topseed_np_scaled <- zscore_topseed_np_scaled
sum(zscoreWeighted2$zscore_topseed_np_scaled==0)
zscoreWeighted2 <- zscoreWeighted2[zscoreWeighted2$zscore_topseed_np_scaled!=0, ]
zscoreWeighted2$rank_np2 <- BiocGenerics::rank(-zscoreWeighted2$zscore_topseed_np_scaled, ties.method="first")

quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank_np2, 2), table)
sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2) )
quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank, 2), table)
sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2) )
quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank_np2, 4), table)
sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2) )
quantile4tab <- tapply(zscoreWeighted2$Celltype2, cut(zscoreWeighted2$rank, 4), table)
sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2) )


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# test for k used for knn-graph construction
source("/Users/fyu/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vijay/script/general_script/scVAR/general_utils.r")
mutualknn30 <- getmutualknn(lsimat30, 30)
mutualknn30_graph <- graph_from_adjacency_matrix(mutualknn30, mode = "undirected", diag = F)

mutualknn5 <- getmutualknn(lsimat30, 5)
mutualknn5_graph <- graph_from_adjacency_matrix(mutualknn5, mode = "undirected", diag = F)

mutualknn10 <- getmutualknn(lsimat30, 10)
mutualknn10_graph <- graph_from_adjacency_matrix(mutualknn10, mode = "undirected", diag = F)

mutualknn20 <- getmutualknn(lsimat30, 20)
mutualknn20_graph <- graph_from_adjacency_matrix(mutualknn20, mode = "undirected", diag = F)

mutualknn50 <- getmutualknn(lsimat30, 50)
mutualknn50_graph <- graph_from_adjacency_matrix(mutualknn50, mode = "undirected", diag = F)

qq <- 0.05
seed_p0.05_top <- zscoreWeighted$rank<=1000*qq
zscore_topseed_np <- randomWalk_sparse(intM=mutualknn30, queryGenes=(1:nrow(mutualknn30))[seed_p0.05_top], gamma=0.05)
zscore_topseed_np_scaled <- zscore_topseed_np %>% capOutlierQuantile(., 0.99) %>% max_min_scale
plot(1:1000, zscore_topseed_np_scaled)
sum(zscore_topseed_np[1:500]==0)
sum(zscore_topseed_np[501:1000]==0)

graph_list <- list(mutualknn5, mutualknn10, mutualknn20, mutualknn30, mutualknn50)
mm <- rep(0, 50)
nn <- rep(0, 50)

aa <- list(); length(aa) <- 50
for (i in 1:length(graph_list)){
    seed_p0.05_top <- zscoreWeighted$rank<=1000*0.05
    zscore_topseed_np <- randomWalk_sparse(intM=graph_list[[i]], queryGenes=(1:nrow(graph_list[[i]]))[seed_p0.05_top], gamma=0.05)
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
    # nn[i] <- sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2))[1]
    # aa[[i]] <- table(zscoreWeighted$Celltype2[seed_p0.05_top])
}
# k=5 (355/(355+75))=0.826; (309/(355+121))=0.650
np score
$`(0.212,395]`

Mono   NK 
 355   40 

$`(395,790]`

Mono   NK 
  75  319 

raw score
$`(0.001,500]`

Mono   NK 
 309  104 

$`(500,1e+03]`

Mono   NK 
 121  255 

# k=10 (407/(407+49))=0.893; (327/(327+129))=0.717
np score
$`(0.093,454]`

Mono   NK 
 407   47 

$`(454,909]`

Mono   NK 
  49  405 

raw score
$`(0.001,500]`

Mono   NK 
 327  126 

$`(500,1e+03]`

Mono   NK 
 129  326 

# k=20 (423/(423+26))=0.942; (324/(324+125))=0.72
np score
$`(0.066,468]`

Mono   NK 
 423   45 

$`(468,936]`

Mono   NK 
  26  441 

raw score
$`(0.001,500]`

Mono   NK 
 324  135 

$`(500,1e+03]`

Mono   NK 
 125  351 

# k=30 (442/(442+8))=0.972; (327/(327+123))=0.72
np score
$`(0.061,470]`

Mono   NK 
 442   28 

$`(470,941]`

Mono   NK 
   8  462 

raw score
$`(0.001,500]`

Mono   NK 
 327  135 

$`(500,1e+03]`

Mono   NK 
 123  355 

# k=50 (441/(441+8))=0.97; (327/(327+123))=0.72
np score
$`(0.054,474]`

Mono   NK 
 441   33 

$`(474,948]`

Mono   NK 
   8  465 

raw score
$`(0.001,500]`

Mono   NK 
 324  138 

$`(500,1e+03]`

Mono   NK 
 125  360 


library(ggplot2)
df <- data.frame(Numberofk=1:5, TPR=c(0.826, 0.893, 0.942, 0.972, 0.97)*100)
p<-ggplot(data=df, aes(x=Numberofk, y=TPR)) +
  geom_bar(stat="identity", fill="steelblue") + pretty_plot() + geom_text(aes(label=TPR), vjust=1.6, color="white", size=3.5)+ coord_flip()
ggsave(p, file = "./fig_ok/robustness-knnk.pdf", width = 4, height = 3)


source("/Users/fyu/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vijay/script/general_script/scVAR/general_utils.r")
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

    source("/Users/fyu/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vijay/script/general_script/scVAR/general_utils.r")
    lsimat30=matSVD(mat=test_q40_n200_c9[idx_no0, ], nComponents = 50, binarize = TRUE) # 1000    50
    mutualknn30 <- getmutualknn(lsimat30, 30)
    mutualknn30_graph <- graph_from_adjacency_matrix(mutualknn30, mode = "undirected", diag = F)
    seed_p0.05_top <- zscoreWeighted$rank<=1000*0.05
    zscore_topseed_np <- randomWalk_sparse(intM=mutualknn30, queryGenes=(1:nrow(mutualknn30))[seed_p0.05_top], gamma=0.05)
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
    # nn[i] <- sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2))[1]
    # aa[[i]] <- table(zscoreWeighted$Celltype2[seed_p0.05_top])
}
1000

Mono   NK 
 233  235 

$`(468,937]`

Mono   NK 
 213  255 

raw score
$`(0.001,500]`

Mono   NK 
 267  195 

$`(500,1e+03]`

Mono   NK 
 179  295 

2500
297/(297+154)=0.659
np score
$`(0.06,471]`
Mono   NK 
 297  174 
$`(471,942]`
Mono   NK 
 154  316 
275/(275+176)=0.610
raw score
$`(0.001,500]`
Mono   NK 
 275  189 
$`(500,1e+03]`

Mono   NK 
 176  301 

5000
361/(361+76)=0.82
np score
$`(0.074,464]`

Mono   NK 
 361  103 

$`(464,928]`

Mono   NK 
  76  387 
283/(283+154)=0.648
raw score
$`(0.001,500]`

Mono   NK 
 283  165 

$`(500,1e+03]`

Mono   NK 
 154  325 

7500
405/(405+31)=0.929
np score
$`(0.075,464]`

Mono   NK 
 405   58 

$`(464,927]`

Mono   NK 
  31  432 
298/(298+138)=0.683
raw score
$`(0.001,500]`

Mono   NK 
 298  150 

$`(500,1e+03]`

Mono   NK 
 138  340 

10000
442/(442+8)=0.982
np score
$`(0.061,470]`

Mono   NK 
 442   28 

$`(470,941]`

Mono   NK 
   8  462 
327/(327+123)=0.727
raw score
$`(0.001,500]`

Mono   NK 
 327  135 

$`(500,1e+03]`

Mono   NK 
 123  355 

438/(438+12)=0.973
np score
$`(0.061,470]`

Mono   NK 
 438   32 

$`(470,941]`

Mono   NK 
  12  458 
371/(371+79)=0.824
raw score
$`(0.001,500]`

Mono   NK 
 371   88 

$`(500,1e+03]`

Mono   NK 
  79  402 

50000
445/(445+13)=0.972
np score
$`(0.053,474]`

Mono   NK 
 445   29 

$`(474,949]`

Mono   NK 
  13  461 

raw score
$`(0.001,500]`
423/(423+35)=0.924
Mono   NK 
 423   37 

$`(500,1e+03]`

Mono   NK 
  35  453 

75000
[1] 1000
[info] fast knn
Stationary step: 205
Stationary Delta: 9.86049399289457e-06
[1] 58
np score
$`(0.059,472]`

Mono   NK 
 442   29 

$`(472,943]`

Mono   NK 
  10  461 

raw score
$`(0.002,500]`

Mono   NK 
 435   20 

$`(500,1e+03]`

Mono   NK 
  17  470 

1e+05
[1] 1000
[info] fast knn
Stationary step: 205
Stationary Delta: 9.80026580950041e-06
[1] 567
np score
(0.568,217]   (217,433] 
        217         216 
raw score
(0.477,262]   (262,525] 
        229         204 

# c(2500, 2500, 5000, 5000, 7500, 7500, 10000, 10000, 25000, 25000, 50000, 50000)
df3 <- data.frame(reads=c("2500", "2500", "5000", "5000", "7500", "7500", "10000", "10000", "25000", "25000", "50000", "50000"), 
                  TPR=c(0.610, 0.659, 0.648, 0.82, 0.683, 0.929, 0.727, 0.982, 0.824, 0.973, 0.924, 0.972)*100,
                  type=c("Z score", "SCAVENGE", "Z score", "SCAVENGE", "Z score", "SCAVENGE", "Z score", "SCAVENGE", "Z score", "SCAVENGE", "Z score", "SCAVENGE"))
df3$type <- factor(df3$type, levels=c("Z score", "SCAVENGE") )
df3$reads <- factor(df3$reads, levels=c("2500", "5000", "7500", "10000", "25000", "50000") )

p <- ggplot(df3, aes(x=reads, y=TPR, fill=type)) + 
   geom_bar(stat="identity", position=position_dodge())
  
p1 <- p + scale_fill_brewer(palette="Paired") + pretty_plot()
p2 <- p + scale_fill_brewer(palette="Reds") + pretty_plot()
p3 <- p + scale_fill_brewer(palette="Greens") + pretty_plot()


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# We set q to various values including 0.1, 0.3, 0.5 to test the robustness to noise
q_number <- c(0.3, 0.35, 0.4, 0.45, 0.5)
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

    source("/Users/fyu/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vijay/script/general_script/scVAR/general_utils.r")
    lsimat30=matSVD(mat=test_q40_n200_c9[idx_no0, ], nComponents = 50, binarize = TRUE) # 1000    50
    mutualknn30 <- getmutualknn(lsimat30, 30)
    mutualknn30_graph <- graph_from_adjacency_matrix(mutualknn30, mode = "undirected", diag = F)
    seed_p0.05_top <- zscoreWeighted$rank<=1000*0.05
    zscore_topseed_np <- randomWalk_sparse(intM=mutualknn30, queryGenes=(1:nrow(mutualknn30))[seed_p0.05_top], gamma=0.05)
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
    # nn[i] <- sapply(quantile4tab, function(x) round(100*x[1]/(x[2]+x[1]), 2))[1]
    # aa[[i]] <- table(zscoreWeighted$Celltype2[seed_p0.05_top])
}
0.25
$`(0.061,470]`
Mono   NK 
 442   28 
$`(470,941]`
 NK 
470 
324/(324+118)=0.733
$`(0.001,500]`
Mono   NK 
 324  137 
$`(500,1e+03]`
Mono   NK 
 118  361 

0.3
np score
$`(0.063,470]`
Mono   NK 
 440   29 
$`(470,939]`
Mono   NK 
   8  461 

raw score
$`(0.001,500]`
Mono   NK 
 326  135 
$`(500,1e+03]`
Mono   NK 
 122  355 

0.35 
364/(364+74)=0.831
Mono   NK 
 364   98 
$`(462,925]`
Mono   NK 
  74  388 
309/(309+129)=0.705
raw score
$`(0.001,500]`
Mono   NK 
 309  141 
$`(500,1e+03]`
Mono   NK 
 129  345 

0.4
320/(320+97)=0.767
Mono   NK 
 320  136 
$`(456,912]`
Mono   NK 
  97  358 

301/(301+116)=0.72
raw score
$`(0.001,500]`
Mono   NK 
 301  149 

$`(500,1e+03]`
Mono   NK 
 116  345 

282/(282+138)=0.671
0.45
np score
$`(0.096,453]`

Mono   NK 
 282  171 
$`(453,906]`
Mono   NK 
 138  314 

288/(288+132)=0.685
raw score
$`(0.001,500]`
Mono   NK 
 288  155 
$`(500,1e+03]`
Mono   NK 
 132  330 

0.5
np score
$`(0.105,448]`

Mono   NK 
 202  246 

$`(448,897]`

Mono   NK 
 198  250 

raw score
$`(0.001,500]`
Mono   NK 
 280  161 
$`(500,1e+03]`
Mono   NK 
 120  335 

df3 <- data.frame(noise=c(0.25, 0.25, 0.3, 0.3, 0.35, 0.35, 0.4, 0.4), 
                  TPR=c(0.733, 1, 0.72, 0.972, 0.705, 0.831, 0.71, 0.767)*100,
                  type=c("Z score", "SCAVENGE", "Z score", "SCAVENGE", "Z score", "SCAVENGE", "Z score", "SCAVENGE"))
df3$type <- factor(df3$type, levels=c("Z score", "SCAVENGE") )
p <- ggplot(df3, aes(x=noise, y=TPR, fill=type)) + 
   geom_bar(stat="identity", position=position_dodge())
  
p1 <- p + scale_fill_brewer(palette="Paired") + pretty_plot()
p2 <- p + scale_fill_brewer(palette="Reds") + pretty_plot()


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
mutualknn30 <- getmutualknn(lsi_mat, 30)

### Network propagation  
np_score <- randomWalk_sparse(intM=mutualknn30, queryGenes=rownames(mutualknn30)[seed_idx], gamma=0.05)
# **Trait relevant score (TRS) with scaled and normalized**  
A few cells are singletons are removed from further analysis
omit_idx <- np_score==0
sum(omit_idx)
mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
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
mutualknn30_graph <- graph_from_adjacency_matrix(mutualknn30, mode = "undirected", diag = F)
plot.igraph(mutualknn30_graph, vertex.size=0.8, vertex.label=NA, vertex.color=adjustcolor("#c7ce3d", alpha.f = 1), vertex.frame.color=NA, 
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
