library(readr)
library(ggfortify)
library(decoupleR)
library(factoextra)
library(GSEABase)
library(readxl)
library(tibble)
library(proxy)
library(data.table)
library(extrafont)
library(OmnipathR)
# Identify number of cell line clusters -----

# Load transcriptomics
RNA <- as.data.frame(read_csv("data/RNA/RNA_log2_FPKM_clean.csv"))
rownames(RNA) <- RNA[,1]
RNA <- RNA[,-1]

# Keep genes where less than 1/3 of values are missing
RNA <- RNA[rowSums(is.na(RNA)) < (dim(RNA)[2]/3),]

# Normalize gene data across cell lines 
SDs <- apply(RNA,1,function(x){sd(x,na.rm = T)})
means <- rowMeans(RNA, na.rm = T)
RNA <- (RNA - means) / SDs

# Load metadata
RNA_metadata <- read_csv(file = "data/metadata/RNA_metadata.csv")

# TF clustering with normalized weighted mean approach (dorothea + decoupleR)
## First load ressource
# dorothea_df <- decoupleR::get_collectri()
load("support/dorothea_df.RData")
## Calculate TF activities per cell line
TF_activity <- apply(RNA,2,function(x){
  x <- as.data.frame(x[which(!is.na(x))])
  TFs <- run_ulm(as.matrix(x), network = dorothea_df, minsize = 20)
  TFs <- as.data.frame(TFs)
  TFs <- TFs[which(TFs$statistic == "ulm"),c(2,4)]
  as_input <- TFs[,2]
  names(as_input) <- TFs[,1]
  return(as_input)
})

## Convert into data.frame
TF_activity_df <- data.frame(t(data.table::rbindlist(
  lapply(TF_activity, function(x) data.table::data.table(t(x))),
  fill = TRUE)))
colnames(TF_activity_df) <- names(TF_activity)

# ## Scale the TF activities (obsolete)
# SDs <- apply(TF_activity_df,1,function(x){sd(x,na.rm = T)})
# means <- rowMeans(TF_activity_df, na.rm = T)
# TF_activity_scaled_df <- (TF_activity_df - means) / SDs
# 

#We don't scale the TF themeselves since the RNA is already scaled
TF_activity_scaled_df <- TF_activity_df

# ## Convert values to ranks (from lowest to highest activity) with only complete cases
TF_activity_scaledrank_df <- as.data.frame(apply(TF_activity_df[complete.cases(TF_activity_df),], 2, base::rank))
# 
## Visualize correlation of TFs between different cell lines
TF_activity_scaledrank_df_corr <- TF_activity_scaledrank_df %>%
  as.matrix %>%
  cor %>%
  as.data.frame %>%
  rownames_to_column(var = 'var1') %>%
  gather(var2, value, -var1) %>%
  pivot_wider(names_from = var1, values_from = value) %>%
  as.data.frame()
rownames(TF_activity_scaledrank_df_corr) <- TF_activity_scaledrank_df_corr$var2
TF_activity_scaledrank_df_corr <- TF_activity_scaledrank_df_corr[,-1]
pheatmap(TF_activity_scaledrank_df_corr)

## Determine optimal number of clusters
df <- as.data.frame(t(TF_activity_df[complete.cases(TF_activity_df),]))
### Silhouette method
fviz_nbclust(df, kmeans, method = "silhouette", k.max = 10)+
  labs(title="") +
  geom_vline(xintercept = 4, linetype = 2) +
  theme_bw() + 
  theme(axis.title.x = element_text(family = "Times New Roman"), 
        axis.title.y = element_text(family = "Times New Roman"), 
        axis.text = element_text(family = "Times New Roman"))

### Elbow method
fviz_nbclust(df, kmeans, method = "wss", k.max = 30) +
  geom_vline(xintercept = 3, linetype = 2) +
  geom_vline(xintercept = 4, linetype = 2) +
  geom_vline(xintercept = 6, linetype = 2) +
  geom_vline(xintercept = 9, linetype = 2)

## --> 3 clusters seems to be an alright number based on silhouet, elbow and interpretability needs.

## Perform k-mean clustering and assign cell lines to clusters (here: explore different numbers based on peaks in Silhouette plot)
RNA_kmeans <- list()
SD <- list()
clusters_df <- list()
heatmap_avgrank_reduced <- list()
heatmap_avgNES_reduced <- list()
heatmap_cluster <- list()
pc_plot <- list()
RNA_metadata_cluster <-  list()
TF_activity_cluster_avgNES <-  list()
TF_activity_cluster_avgNES_reduced <- list()
TF_activity_cluster_avgrank <- list()
TF_activity_cluster_avgrank_reduced <- list()

# for(i in c(3,4,6,9)) {
for(i in c(3)) {  
  #Set seed so assignment to clusters is consistent across multiple runs
  set.seed(12) 
  print(i)
  RNA_kmeans[[i]] <- kmeans(df, i, iter.max = 10000, nstart = 100)
  
  ## Extract clusters
  clusters <- data.frame("cluster" = RNA_kmeans[[i]]$cluster)
  
  ## Inspect clustering
  RNA_kmeans[[i]]$cluster <- data.frame("Cell_line"=names(RNA_kmeans[[i]]$cluster),
                                        "Cluster"=as.character(RNA_kmeans[[i]]$cluster))
  
  pc_plot[[i]] <- autoplot(prcomp(t(TF_activity_scaled_df[complete.cases(TF_activity_scaled_df),])), 
           data = RNA_kmeans[[i]]$cluster,
           colour = "Cluster") +
    theme_bw() +
    scale_color_manual(values=c("1"="red", "2"="cyan", "3"="orange")) +
    scale_fill_manual(values=c("1"="red", "2"="cyan", "3"="orange")) + 
    theme(axis.title.x = element_text(family = "Times New Roman"), 
          axis.title.y = element_text(family = "Times New Roman"), 
          axis.text = element_text(family = "Times New Roman"), 
          legend.text = element_text(family = "Times New Roman"), 
          legend.title = element_text(family = "Times New Roman"))
  
  ## Size of clusters (amount of cell lines per cluster)
  table(RNA_kmeans[[i]]$cluster)
  table(clusters)
  
  ## Cell line TF activity with cluster assignment
  anno_colors <- list(
    "Cluster"=c("1"="red", "2"="cyan", "3"="orange")
  )
  anno_col <- data.frame("Cluster"=as.character(RNA_kmeans[[i]]$cluster$Cluster), 
                         row.names = RNA_kmeans[[i]]$cluster$Cell_line)
  
  heatmap_cluster[[i]] <- pheatmap(TF_activity_scaled_df[,order(RNA_kmeans[[i]]$cluster$Cluster)], 
                                   annotation_col = anno_col, 
                                   annotation_colors = anno_colors,
                                   cluster_cols = F, 
                                   show_rownames = F, 
                                   show_colnames = F,
                                   fontfamily = "Times New Roman",
                                   fontsize = 12, 
                                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
                                   breaks = seq(-4,4,8/100))
  
  ## Compare average of TF values per cluster
  ggpubr::compare_means(value ~ Cluster,
                merge(cbind(pivot_longer(TF_activity_scaled_df,
                                         cols = colnames(TF_activity_scaled_df),
                                         names_to = "sample",
                                         values_to = "value"), 
                            "feature"=rep(rownames(TF_activity_scaled_df),60)), 
                      RNA_kmeans[[i]]$cluster, by.x = "sample", by.y = "Cell_line"),
                method="t.test",
                p.adjust.method = "BH",
                paired = F,
                alternative = "two.sided")
  
  ## Compare to cell line metadata
  clusters_df[[i]] <- cbind("cell_line"=rownames(clusters),clusters)
  RNA_metadata_cluster[[i]] <- merge(clusters_df[[i]], RNA_metadata, by="cell_line")
  
  ## TF activity per cluster (calculate mean of cell lines per cluster)
  ### Based on TF activity ranks
  TF_activity_cluster_avgrank[[i]] <- as.data.frame(sapply(sort(unique(clusters$cluster)), function(x, TF_activity_scaledrank_df, clusters){
    df <- as.data.frame(TF_activity_scaledrank_df[,clusters[clusters$cluster == x, 'cell_line']])
    rownames(df) <- rownames(TF_activity_scaledrank_df)
    ifelse(!is.null(dim(df)), return(rowMeans(df, na.rm = T)), return(df))
  }, TF_activity_scaledrank_df = TF_activity_scaledrank_df, clusters = clusters_df[[i]], USE.NAMES = T))
  
  ### Based on TF activity values
  TF_activity_cluster_avgNES[[i]] <- as.data.frame(sapply(sort(unique(clusters$cluster)), function(x, TF_activity_scaled_df, clusters){
    df <- as.data.frame(TF_activity_scaled_df[,clusters[clusters$cluster == x, 'cell_line']])
    rownames(df) <- rownames(TF_activity_scaled_df)
    ifelse(!is.null(dim(df)), return(rowMeans(df, na.rm = T)), return(df))
  }, TF_activity_scaled_df = TF_activity_scaled_df, clusters = clusters_df[[i]], USE.NAMES = T))
  
  ### Keep top 25 variable genes (based on ranks = scale-free)
  SDs <- apply(TF_activity_cluster_avgrank[[i]], 1, sd)
  SDs <- sort(SDs, decreasing = T)[1:25]
  SD[[i]] <- SDs
  
  ## Reorder based on SD value
  TF_activity_cluster_avgrank_reduced[[i]] <- TF_activity_cluster_avgrank[[i]][names(SDs),]
  TF_activity_cluster_avgNES_reduced[[i]] <- TF_activity_cluster_avgNES[[i]][names(SDs),]
  
  ## Inspect results
  heatmap_avgrank_reduced[[i]] <- pheatmap(TF_activity_cluster_avgrank_reduced[[i]], cluster_rows = F, cluster_cols = F, display_numbers = T)
  
  colnames(TF_activity_cluster_avgNES_reduced[[i]]) <- c("Cluster 1","Cluster 2","Cluster 3")
  heatmap_avgNES_reduced[[i]] <- pheatmap(TF_activity_cluster_avgNES_reduced[[i]], 
                                          cluster_rows = T, 
                                          cluster_cols = F, 
                                          display_numbers = F,
                                          fontfamily = "Times New Roman",
                                          fontsize = 12,
                                          angle_col = 0,
                                          color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100))
  
  ggpubr::compare_means(value ~ Cluster,
                pivot_longer(TF_activity_cluster_avgNES_reduced[[i]],
                             cols=colnames(TF_activity_cluster_avgNES_reduced[[i]]),
                             values_to = "value", 
                             names_to = "Cluster"),
                method="t.test",
                p.adjust.method = "BH",
                paired = F,
                alternative = "two.sided")
}

# Print number of samples per cluster
lapply(clusters_df[!unlist(lapply(clusters_df,is.null))], function(x) table(x$cluster))

# Compare clusters by top TFs by similarity (Jaccard distance)
# overlap <- data.frame(rbind(NA,data.table::rbindlist(
#   lapply(SD[!unlist(lapply(SD,is.null))], function(x) data.table::data.table(t(x))),
#   fill = TRUE
# ), fill=T))[-1,-1]
# overlap[is.na(overlap)] <- 0
# overlap[overlap>0] <- 1
# rownames(overlap) <- c("3","4","6","9")
# colnames(overlap) <- c("3","4","6","9")
# pheatmap(simil(overlap, overlap, method = "Jaccard"), cluster_cols = F, cluster_rows = F, display_numbers = T)

# Compare clusters by cell line assignment
cluster_number <- 3
tissue_origin_cluster <- RNA_metadata_cluster[[cluster_number]][,c(2:3,8)]
colnames(tissue_origin_cluster) <- c("cluster", "tissue", "histology")
rownames(tissue_origin_cluster) <- RNA_metadata_cluster[[cluster_number]][,1]
pheatmap(TF_activity_scaled_df[,order(RNA_metadata_cluster[[cluster_number]]$cluster)], annotation_col = tissue_origin_cluster, cluster_cols = F, show_rownames = F, show_colnames = T)

cell_clusters <- lapply(clusters_df[!unlist(lapply(clusters_df,is.null))], function(x) lapply(clusters_df[!unlist(lapply(clusters_df,is.null))], function(y){
  d <- simil(table(reshape2::melt(x, id.var = -1)[-2]),
      table(reshape2::melt(y, id.var = -1)[-2]),
      method = "Jaccard")
  if(dim(d)[1]!=dim(d)[2]) pheatmap(d, cluster_rows = F, cluster_cols = F)
}))

# Inspect cluster by heatmap
lapply(pc_plot[!unlist(lapply(pc_plot,is.null))], function(x) x)

plot.new()
heatmap_cluster[[3]]
# plot.new()
# heatmap_cluster[[4]]
# plot.new()
# heatmap_cluster[[6]]
# plot.new()
# heatmap_cluster[[9]]

plot.new()
heatmap_avgNES_reduced[[3]]
# plot.new()
# heatmap_avgNES_reduced[[4]]
# plot.new()
# heatmap_avgNES_reduced[[6]]
# plot.new()
# heatmap_avgNES_reduced[[9]]

# Save optimal cluster result (here: 3)
write_csv(RNA_metadata_cluster[[3]], file = "data/metadata/RNA_metadata_cluster.csv")


# Justify reduction of genes -----
library(readr)
library(ggfortify)
library(decoupleR)
library(factoextra)
library(GSEABase)
library(readxl)

# Load transcriptomics
RNA <- as.data.frame(read_csv("data/RNA/RNA_log2_FPKM_clean.csv"))
rownames(RNA) <- RNA[,1]
RNA <- RNA[,-1]

# Remove genes with excessive amount of NAs (only keep genes with max. amount of NAs = 33.3 % across cell lines)
RNA <- RNA[rowSums(is.na(RNA))<(dim(RNA)[2]/3),]


# TF clustering with normalized weighted mean approach (dorothea + decoupleR)
## First load ressource
dorothea_df <- decoupleR::get_dorothea(levels = c("A","B","C"))

# Extract number of TFs and mean number of targets per TF depending on gene data set size
results_table <- list()
for(i in c(1:dim(RNA)[1])){
  
  RNA_sd <- sort(apply(RNA, 1, function(x) sd(x,na.rm = T)), decreasing = T)
  RNA_sd <- RNA_sd[1:i]
  RNA_reduced <- RNA[names(RNA_sd),]
  
  ## Normalize gene data across cell lines 
  SDs <- apply(RNA_reduced,1,function(x){sd(x,na.rm = T)})
  means <- rowMeans(RNA_reduced, na.rm = T)
  RNA_reduced <- (RNA_reduced - means) / SDs
  
  ## Calculate TF activities per cell line
  TF_activity_reduced <- apply(RNA_reduced,2,function(x){
    x <- as.data.frame(x[which(!is.na(x))])
    dorothea_df_intersected <- as.data.frame(intersect_regulons(mat = x, network = dorothea_df, minsize = 20, .source = source, .target = target))[,c(1,3)]
  })
  
  ## Extract overlapping TFs
  overlap_TF <- lapply(c(1:length(TF_activity_reduced)),function(x) unique(TF_activity_reduced[[x]]$source)) %>%
    melt() %>%
    count(value)
  overlap_TF <- unique(overlap_TF[overlap_TF$n==length(TF_activity_reduced),1])
  TF_activity_reduced <- lapply(c(1:length(TF_activity_reduced)),function(x) TF_activity_reduced[[x]][TF_activity_reduced[[x]]$source %in% overlap_TF,])
  
  ## Calculate mean of TF targets and SD considering all cell lines 
  TF_activity_reduced_df <- data.frame(t(mapply(function(x) data.frame("mean"=mean(table(TF_activity_reduced[[x]]$source)), "sd"=sd(table(TF_activity_reduced[[x]]$source))), c(1:length(TF_activity_reduced)))))
  
  ## Save how many genes were considered, mean and SD of TF targets as well as numer of identified TFs
  results_table[[i]] <- data.frame("mean"=mean(melt(TF_activity_reduced_df$mean)[,1]),
                                   "sd"=mean(melt(TF_activity_reduced_df$sd)[,1]),
                                   "Number_of_TFs"=length(overlap_TF),
                                   "Number_of_genes"=i)
  results_table[[i]][is.na(results_table[[i]])] <- 0 
  print(i)
}

#save(results_table,file = "results/RNA/TF_TFtarget_gene_reduction.Rdata")
load(file = "results/RNA/TF_TFtarget_gene_reduction.Rdata")

## From list to data frame
results_table_df <- data.frame()
for (i in c(1:length(results_table))) {
  results_table_df <- rbind(results_table_df, results_table[[i]])
}
## Plot result with two y axis (function of TF targets and TFs depending on number of genes)

plot <- ggplot(results_table_df, aes(Number_of_genes,Number_of_TFs)) +
  geom_point(color = "blue", size = 0.1) +
  geom_line(aes(color = "Number of TFs")) +
  geom_line(aes(y=mean+sd, color = "SD of TF targets per TF")) +
  geom_point(aes(y = mean, color = "red"), size = 0.1) +
  geom_line(aes(y = mean, color = "Average TF targets per TF")) +
  scale_y_continuous("Number of TFs", breaks = seq(0,230,10), sec.axis = sec_axis(~ ., name = "Average TF targets per TF", breaks = seq(0,230,10))) +
  scale_x_reverse("Number of genes") +
  scale_color_manual(name = "", values = c("Number of TFs" = "blue", "Average TF targets per TF" = "red", "SD of TF targets per TF" = "#FFCCCB")) +
  ggtitle("Number of identified transcription factors and targets (mean) depending on number of genes") +
  theme_bw() +
  theme(axis.title.y.left =element_text(colour = "blue"),axis.title.y.right =element_text(colour = "red"),
        axis.title.x = element_text(family = "Times New Roman"), 
        axis.title.y = element_text(family = "Times New Roman"), 
        axis.text = element_text(family = "Times New Roman"), 
        legend.text = element_text(family = "Times New Roman"), 
        legend.title = element_text(family = "Times New Roman"))

ggsave(file="results/RNA/TF_TFtarget_gene_reduction.pdf", device = "pdf", plot = plot, width = 25, height = 20)



