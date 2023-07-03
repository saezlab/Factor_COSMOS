In order to evaluate the data obtained from the MOFA-COSMOS pipeline
downstream, further analyses can be performed.

``` r
library(decoupleR)
library(liana)
library(moon)
library(readr)
library(pheatmap)
library(dplyr)
library(reshape2)
library(MOFA2)
library(ggplot2)
library(RCy3)
library(RColorBrewer)
library(cosmosR)
```

## Calculate cell-line specific transcriptomics

Firstly, we can calculate cell-line specific transcriptomic values by
multiplying the MOFA feature weights with the MOFA factor weights (here:
Factor 4). While this can seemlessly be extended to metabolomic data (or
any other omics for which we have MOFA weights), we will focus on the
RNA layer for the sake of simplicity.

``` r
# RNA cell-line-specific values (RNA_specific = MOFA_feature_weight * MOFA_factor_weight)

## Get RNA raw values
RNA_raw <- as.data.frame(read_csv("data/RNA/RNA_log2_FPKM_clean.csv"))
```

    ## Rows: 11265 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (1): Genes
    ## dbl (60): 786-0, A498, A549/ATCC, ACHN, BT-549, CAKI-1, CCRF-CEM, COLO 205, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
rownames(RNA_raw) <- RNA_raw[,1]
RNA_raw <- RNA_raw[,-1]
RNA_raw <- melt(as.data.frame(cbind(RNA_raw,row.names(RNA_raw))))
```

    ## Using row.names(RNA_raw) as id variables

``` r
RNA_raw$view <- "RNA"
RNA_raw <- RNA_raw[,c(2,1,4,3)]                 
names(RNA_raw) <- c("sample","feature","view","value")

## Get MOFA feature and factor weights from Factor 4  
model <- load_model('results/mofa/mofa_res_10factor.hdf5')
```

    ## Warning in load_model("results/mofa/mofa_res_10factor.hdf5"): There are duplicated features names across different views. We will add the suffix *_view* only for those features 
    ##             Example: if you have both TP53 in mRNA and mutation data it will be renamed to TP53_mRNA, TP53_mutation

    ## Warning in .quality_control(object, verbose = verbose): Factor(s) 1 are strongly correlated with the total number of expressed features for at least one of your omics. Such factors appear when there are differences in the total 'levels' between your samples, *sometimes* because of poor normalisation in the preprocessing steps.

``` r
meta_data <- read_csv("data/metadata/RNA_metadata_cluster.csv")[,c(1,2)]
```

    ## Rows: 60 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (12): cell_line, tissue of origin a, sex a, prior treatment a,b, Epithel...
    ## dbl  (4): cluster, age a, mdr f, doubling time g
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
colnames(meta_data) <- c("sample","cluster")
samples_metadata(model) <- meta_data[meta_data$sample %in% samples_metadata(model)$sample,]

weights <- get_weights(model, views = "all", factors = "all")
RNA_MOFA <- data.frame(weights$RNA[,4])
RNA_MOFA <- cbind(rownames(RNA_MOFA), RNA_MOFA)
colnames(RNA_MOFA) <- c("feature", "value")
RNA_MOFA$feature <- gsub("_RNA","",RNA_MOFA$feature)

weights <- data.frame(get_factors(model, factors = "all")$single_group[,4])
Factor_MOFA <- cbind(rownames(weights), weights)
colnames(Factor_MOFA) <- c("sample", "value")

## Calculate score
RNA_weight_cell_line <- merge(RNA_raw, RNA_MOFA, by.x='feature', by.y='feature')
RNA_weight_cell_line <- merge(RNA_weight_cell_line, Factor_MOFA, by.x='sample', by.y='sample')

RNA_weight_cell_line <- cbind(RNA_weight_cell_line, "final" = (RNA_weight_cell_line$value.y * RNA_weight_cell_line$value))
colnames(RNA_weight_cell_line) <- c("sample","feature","view","RNA_raw","RNA_MOFA","sample_MOFA","Final")

nodes_cellline_weight <- dcast(RNA_weight_cell_line, feature~sample, value.var = "Final")
rownames(nodes_cellline_weight) <- nodes_cellline_weight$feature
nodes_cellline_weight <- nodes_cellline_weight[,-1]
```

Then we can further select a subset of nodes to focus on, here we select
interesting nodes that are well representing transcriptional regulations
(e.g.  HIF1A, MYC, STAT1), PTM regulations (e.g. JAK1, ABL1),
ligand-receptor interactions (e.g. VEGFA-ITGB1) and metabolic
regulations (e.g. ACO1 and citrate).

``` r
## Select nodes 
interesting_nodes <- data.frame("name" = c("VEGFA","HIF1A","PRKCA","FYN","ILK","ACO1","JAK1","ABL1","ITGB1","STAT1","AKT1","MYC","TKT"))
nodes_cellline_weight_filtered <- nodes_cellline_weight[rownames(nodes_cellline_weight) %in% interesting_nodes$name,]
```

And add the meta data information.

``` r
RNA_metadata_cluster <- as.data.frame(read_csv(file = "data/metadata/RNA_metadata_cluster.csv"))
```

    ## Rows: 60 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (12): cell_line, tissue of origin a, sex a, prior treatment a,b, Epithel...
    ## dbl  (4): cluster, age a, mdr f, doubling time g
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
tissue_origin_cluster <- as.data.frame(RNA_metadata_cluster[,c(2,3)]) #:2,8
colnames(tissue_origin_cluster) <- c("cluster","Tissue") #", , "histology")
rownames(tissue_origin_cluster) <- RNA_metadata_cluster[,1]
tissue_origin_cluster$Tissue <- as.character(tissue_origin_cluster$Tissue)

tissue_origin_cluster <- tissue_origin_cluster[,-1,drop = F]
# tissue_origin_cluster$cluster <- as.character(tissue_origin_cluster$cluster)
```

To visualize the results, heatmaps are shown.

``` r
var1 <- rev(brewer.pal(9, "Set1"))
names(var1) <- unique(tissue_origin_cluster$Tissue)
# var2 <- c("red","cyan","orange")
# names(var2) <- unique(tissue_origin_cluster$cluster)[order(unique(tissue_origin_cluster$cluster))]
# anno_colors <- list(Tissue=var1,
#                     cluster=var2)
anno_colors <- list(Tissue=var1)


pheatmap(nodes_cellline_weight, 
         display_numbers = F, 
         cluster_rows = F, 
         show_rownames =  F, 
         annotation_col = tissue_origin_cluster, 
         annotation_colors =  anno_colors)
```

![](Further_MOFA_COSMOS_analyses_files/figure-gfm/Visualization%20transcriptomics-1.png)<!-- -->

``` r
pheatmap(nodes_cellline_weight_filtered, 
         display_numbers = F,
         annotation_col = tissue_origin_cluster, 
         annotation_colors =  anno_colors)
```

![](Further_MOFA_COSMOS_analyses_files/figure-gfm/Visualization%20transcriptomics-2.png)<!-- -->
Since, we are also interested in comparing the two cell lines with the
most differing features, we use the MOFA factor values to determine
these by calculating the highest difference.

``` r
# Calculate activity scores of most different cell lines

## Get cell line factor weights
factor_weights <- as.data.frame(get_factors(model,factors = 4))
factor_weights$sample <- rownames(factor_weights)
factor_weights <- merge(factor_weights,meta_data,by="sample")

## Identify cell lines with highest difference in factor values
diff_matrix <- data.frame()
for(x in c(1:58)) {
  for(y in c(1:58)) {
    diff_matrix[x,y] <- abs(factor_weights$Factor4[x]-factor_weights$Factor4[y])
  }
}
rownames(diff_matrix) <- factor_weights$sample
colnames(diff_matrix) <- factor_weights$sample

max(diff_matrix) #sample 28 and 44 = MOLT-4 & SF-539
```

    ## [1] 3.564818

``` r
max_cell_lines <- rownames(diff_matrix)[which(diff_matrix == max(diff_matrix), arr.ind = T)][c(1,2)]

plot_factor(model, 
            factors = 4,
            color_by = "cluster",
            dot_size = 3,
            dodge = T,           
            legend = T,          
            add_violin = T,     
            violin_alpha = 0.25) +
  scale_color_manual(values=c("1"="red", "2"="cyan", "3"="orange")) +
  scale_fill_manual(values=c("1"="red", "2"="cyan", "3"="orange")) + 
  geom_label(label=c(ifelse(get_factors(model)$single_group[,4]==max(get_factors(model)$single_group[,4]) | get_factors(model)$single_group[,4]==min(get_factors(model)$single_group[,4]),model@samples_metadata$sample, NA)),
                     show.legend = F, size = 3, position=position_dodge(width=0.9), family = "Times New Roman", label.size = 0.1) +
  theme_bw() +
  ylab("Weight") +
  xlab("") +
  theme(strip.text = element_text(family = "Times New Roman"), axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y.right =element_text(colour = "red"),
        axis.title.x = element_text(family = "Times New Roman"), axis.title.y = element_text(family = "Times New Roman"), axis.text = element_text(family = "Times New Roman"), legend.text = element_text(family = "Times New Roman"), legend.title = element_text(family = "Times New Roman"))
```

    ## Warning: `fct_explicit_na()` was deprecated in forcats 1.0.0.
    ## ℹ Please use `fct_na_value_to_level()` instead.
    ## ℹ The deprecated feature was likely used in the MOFA2 package.
    ##   Please report the issue at <https://github.com/bioFAM/MOFA2>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: Removed 56 rows containing missing values (`geom_label()`).

![](Further_MOFA_COSMOS_analyses_files/figure-gfm/Identify%20distance%20between%20cell%20lines-1.png)<!-- -->
Let’s focus on the interesting node values for these cell lines.

``` r
pheatmap(nodes_cellline_weight_filtered[,max_cell_lines], 
         display_numbers = T,
         annotation_col = tissue_origin_cluster, 
         annotation_colors =  anno_colors,
         cluster_cols = F)
```

![](Further_MOFA_COSMOS_analyses_files/figure-gfm/Visualization%20transcriptomics%20most%20distant%20cell%20line-1.png)<!-- -->

Interestingly, we can see the clear pattern differences for most genes
here at the RNA level, except for FYN, AKT1 and TKT. However, these are
still interesting candidate from a topological perspective, as they can
explain connections between ITGB1, PRKCA and ABL1. Thus it’s possible
that those mechanisms are primarilly driven at the PTM level, instead of
the transcriptional regulation level.

## Use cell-line specific transcriptomics to infer biological activities

Further, we can use the cell-line specific transcriptomic values to
calculate cell-line specific activity scores with decoupleR and prior
knowledge.

``` r
## Calculate cell-line specific activity scores
### LIANA, mOOn, Dorothea

# Load LIANA (receptor and ligand) consensus network
# ligrec_ressource <- distinct(liana::decomplexify(liana::select_resource("Consensus")[[1]]))
load("support/ligrec_ressource.RData")
ligrec_geneset <- ligrec_ressource[,c("source_genesymbol","target_genesymbol")]
ligrec_geneset$set <- paste(ligrec_geneset$source_genesymbol, ligrec_geneset$target_genesymbol, sep = "___")
ligrec_geneset <- reshape2::melt(ligrec_geneset, id.vars = "set")[,c(3,1)]
names(ligrec_geneset)[1] <- "gene"
ligrec_geneset$mor <- 1
ligrec_geneset$likelihood <- 1
ligrec_geneset <- distinct(ligrec_geneset)

# Load Dorothea (TF) network
load("support/dorothea_df.RData")


lig_rec_scores <- list()
dorothea_scores <- list()

for(i in unique(RNA_weight_cell_line$sample)) {
  RNA_sample <- RNA_weight_cell_line[RNA_weight_cell_line$sample == i,7]
  names(RNA_sample) <- RNA_weight_cell_line[RNA_weight_cell_line$sample == i,2]
  
  # Don't use inferred values
  features <- RNA_weight_cell_line[RNA_weight_cell_line$sample==i,"RNA_raw"]
  names(features) <- RNA_weight_cell_line[RNA_weight_cell_line$sample==i,"feature"]
  features <- features[!is.na(features)]
  RNA_sample <- RNA_sample[names(features)]
  RNA_sample <- RNA_sample[!is.na(RNA_sample)]
  
  # Calculate regulatory activities from Receptor and Ligand network
  ligrec_high_vs_low <- run_ulm(mat = as.matrix(RNA_sample), network = ligrec_geneset, .source = set, .target = gene, minsize = 2) 
  ligrec_high_vs_low <- ligrec_high_vs_low[ligrec_high_vs_low$statistic == "ulm",]
  ligrec_high_vs_low_vector <- ligrec_high_vs_low$score
  names(ligrec_high_vs_low_vector) <- ligrec_high_vs_low$source
  
  rec_inputs <- ligrec_high_vs_low_vector
  names(rec_inputs) <- gsub(".+___","",names(rec_inputs))
  rec_inputs <- tapply(rec_inputs, names(rec_inputs), mean)
  receptors <- names(rec_inputs)
  rec_inputs <- as.numeric(rec_inputs)
  names(rec_inputs) <- receptors
  
  lig_rec_scores[[i]] <- rec_inputs
  
  lig_inputs <- ligrec_high_vs_low_vector
  names(lig_inputs) <- gsub("___.+","",names(lig_inputs))
  lig_inputs <- tapply(lig_inputs, names(lig_inputs), mean)
  ligands <- names(lig_inputs)
  lig_inputs <- as.numeric(lig_inputs)
  names(lig_inputs) <- ligands
  
  lig_rec_scores[[i]] <- c(lig_rec_scores[[i]],lig_inputs)
  
  # Calculate regulatory activities from TF network
  TF_high_vs_low <- run_ulm(mat = as.matrix(RNA_sample), network = dorothea_df, minsize = 10)
  TF_high_vs_low <- TF_high_vs_low[TF_high_vs_low$statistic == "ulm",]
  TF_high_vs_low_vector <- TF_high_vs_low$score
  names(TF_high_vs_low_vector) <- TF_high_vs_low$source
  
  dorothea_scores[[i]] <- TF_high_vs_low_vector 
  
  # Prepare moon analysis input
  TF_high_vs_low <- as.data.frame(TF_high_vs_low[,c(2,4)])
  row.names(TF_high_vs_low) <- TF_high_vs_low[,1]
  
}
```

After that, we create a dataframe that includes the activity values for
our nodes inside the sub signaling pathway.

``` r
## Filter interesting nodes 
lig_rec_scores_df <- data.frame(t(data.table::rbindlist(
  lapply(lig_rec_scores, function(x) data.table::data.table(t(x))),
  fill = TRUE)))
colnames(lig_rec_scores_df) <- names(lig_rec_scores)
lig_rec_scores_filtered_df <- lig_rec_scores_df[rownames(lig_rec_scores_df) %in% interesting_nodes$name,]

dorothea_scores_df <- data.frame(t(data.table::rbindlist(
  lapply(dorothea_scores, function(x) data.table::data.table(t(x))),
  fill = TRUE)))
colnames(dorothea_scores_df) <- names(dorothea_scores)
dorothea_scores_filtered_df <- dorothea_scores_df[rownames(dorothea_scores_df) %in% interesting_nodes$name,]

## Create dataframe with activity scores
activity_scores <- rbind(lig_rec_scores_filtered_df,dorothea_scores_filtered_df[!(rownames(dorothea_scores_filtered_df) %in%  rownames(lig_rec_scores_filtered_df)),])

all_activity_scores <- rbind(lig_rec_scores_df,dorothea_scores_df[!(rownames(dorothea_scores_df) %in% rownames(lig_rec_scores_df)),])

activity_scores <- cbind(rownames(activity_scores),activity_scores)
write_csv(activity_scores, file = "data/mofa/activity_scores.csv")

all_activity_scores <- cbind(rownames(all_activity_scores),all_activity_scores)
write_csv(all_activity_scores, file = "data/mofa/all_activity_scores.csv")
```

We can visualize the different activity scores with heatmaps.

``` r
## Visualize results
activity_scores <- as.data.frame(read_csv(file = "data/mofa/activity_scores.csv"))
```

    ## Rows: 5 Columns: 59
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (1): rownames(activity_scores)
    ## dbl (58): 786-0, A498, A549/ATCC, ACHN, BT-549, CAKI-1, CCRF-CEM, COLO 205, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
rownames(activity_scores) <- activity_scores[,1]
activity_scores <- activity_scores[,-1]

all_activity_scores <- as.data.frame(read_csv(file = "data/mofa/all_activity_scores.csv"))
```

    ## Rows: 452 Columns: 59
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (1): rownames(all_activity_scores)
    ## dbl (58): 786-0, A498, A549/ATCC, ACHN, BT-549, CAKI-1, CCRF-CEM, COLO 205, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
rownames(all_activity_scores) <- all_activity_scores[,1]
all_activity_scores <- all_activity_scores[,-1]

pheatmap(activity_scores, 
         display_numbers = F, 
         cluster_rows = T, 
         cluster_cols = T,
         annotation_col = tissue_origin_cluster, 
         annotation_colors =  anno_colors)
```

![](Further_MOFA_COSMOS_analyses_files/figure-gfm/Visualize%20activity%20scores-1.png)<!-- -->

``` r
pheatmap(activity_scores[,colnames(activity_scores) %in% max_cell_lines], 
         display_numbers = T, 
         cluster_rows = T, cluster_cols = F,
         annotation_col = tissue_origin_cluster, 
         annotation_colors =  anno_colors)
```

![](Further_MOFA_COSMOS_analyses_files/figure-gfm/Visualize%20activity%20scores%20for%20specific%20cell%20lines-1.png)<!-- -->

Since we would like to see relative differences between our cell lines,
we calculate the relative activity values here through scaling
($\frac{feature_{single} - mean(feature)}{sd(feature)}$).

``` r
## Transform activity values to relative ones
SDs <- apply(activity_scores,1,function(x){sd(x,na.rm = T)})
means <- rowMeans(activity_scores, na.rm = T)
activity_scores_scaled <- (activity_scores - means) / SDs

SDs <- apply(all_activity_scores,1,function(x){sd(x,na.rm = T)})
means <- rowMeans(all_activity_scores, na.rm = T)
all_activity_scores_scaled <- (all_activity_scores - means) / SDs
```

Let’s visualize the scaled activity scores.

``` r
## Visualize results
pheatmap(activity_scores_scaled, 
         display_numbers = F, 
         cluster_rows = T, 
         cluster_cols = T,
         annotation_col = tissue_origin_cluster, 
         annotation_colors =  anno_colors)
```

![](Further_MOFA_COSMOS_analyses_files/figure-gfm/Visualization%20of%20relative%20activities-1.png)<!-- -->

``` r
pheatmap(activity_scores_scaled[,colnames(activity_scores_scaled) %in% max_cell_lines], 
         display_numbers = T, 
         cluster_rows = T, 
         cluster_cols = F,
         annotation_col = tissue_origin_cluster, 
         annotation_colors =  anno_colors)
```

![](Further_MOFA_COSMOS_analyses_files/figure-gfm/Visualization%20of%20relative%20activities-2.png)<!-- -->

## Visualize results of cell-line specific values with CytoScape

Finally, we can add the activity scores to our network in CytoScape
allowing us to visualize all results nicely. For more flexibility, you
can directly import the SIF and ATT files in cytoscpae and use the
graphical user interface to make your own style.

``` r
## Add weighted, transformed activity scores to network
combined_SIF_reduced <- read_csv(file = "results/cosmos/SIF_mofamoon_combined.csv")
combined_ATT_reduced <- read_csv(file = "results/cosmos/ATT_mofamoon_combined.csv")

activity_scores_scaled <- cbind("id"=rownames(activity_scores_scaled),activity_scores_scaled)

style.name <- "dataStyle_activity"
createVisualStyle(style.name)

activity_weights <- data.frame()
cosmos_weights <- data.frame()

for(i in max_cell_lines) {
  
  nodes <- data.frame(id = combined_ATT_reduced$Nodes, 
                      NodeType = as.integer(combined_ATT_reduced$NodeType), 
                      mofa_weights = as.numeric(combined_ATT_reduced$mofa_weights), 
                      Activity = combined_ATT_reduced$Activity, 
                      stringsAsFactors = F)
  
  nodes <- merge(nodes, activity_scores_scaled[,c("id",i)], by = "id", all = T)
  
  colnames(nodes) <- c("id","NodeType","mofa_weights","Activity","Activity_Scores")
  
  nodes$Activity_Scores[is.na(nodes$Activity_Scores)] <- 0
  nodes$mofa_weights[is.na(nodes$mofa_weights)] <- 0
    
  edges <- data.frame(source=combined_SIF_reduced$Node1, target = combined_SIF_reduced$Node2, sign = combined_SIF_reduced$Sign, weight = combined_SIF_reduced$Weight,stringsAsFactors = F)
  
  ##https://www.genome.jp/entry/N00001
  EGF_EGFR_RAS_ERK_pathway_extended <- c("VEGFA","HIF1A","PRKCA","FYN","ILK","ACO1","JAK1","ABL1","ITGB1","STAT1","AKT1","MYC","TKT")
  
  edges_MAPK_reduced_selected <- edges[edges$target %in% EGF_EGFR_RAS_ERK_pathway_extended & edges$source %in% EGF_EGFR_RAS_ERK_pathway_extended,]
  
  edges_MAPK_reduced_selected_start_end <- edges_MAPK_reduced_selected
  # edges_MAPK_reduced_selected_start_end <- rbind(edges_MAPK_reduced_selected, edges[edges$target %in% "EGFR",],edges[edges$source %in% c("STAT1","MYC"),])
  
  nodes_MAPK_reduced_selected_start_end <- nodes[nodes$id %in% edges_MAPK_reduced_selected_start_end$source | nodes$id %in% edges_MAPK_reduced_selected_start_end$target,]
  
  createNetworkFromDataFrames(nodes_MAPK_reduced_selected_start_end,edges_MAPK_reduced_selected_start_end,title = paste0(i, "_EGF_EGFR_RAS_ERK_pathway_start_end"), table.key.column = "name", collection = "Sub-Network per cell line")
  setVisualStyle(style.name, network = paste0(i, "_EGF_EGFR_RAS_ERK_pathway_start_end"))
  
  # selectNodes(list("VEGFA","HIF1A","PRKCA","FYN","ILK","ACO1","JAK1","ABL1","ITGB1","STAT1","AKT1","MYC","TKT"), by.col = "name", network = paste0(i, "_EGF_EGFR_RAS_ERK_pathway_start_end"))
  # deleteSelectedNodes(network = paste0(i, "_EGF_EGFR_RAS_ERK_pathway_start_end"))
  
  # print("lol")
  
  setEdgeTargetArrowShapeMapping(table.column = "sign",table.column.values =  list("-1","1"), shapes = list("T","ARROW"), style.name = style.name, network = paste0(i, "_EGF_EGFR_RAS_ERK_pathway_start_end"))
  setEdgeLineWidthDefault(2, style.name)
  
  # print("plop")
  
  # duplicated_edges <- as.data.frame(getTableColumns('edge','sign', network = paste0(i, "_EGF_EGFR_RAS_ERK_pathway_start_end")))
  # duplicated_edges[,2] <- rownames(duplicated_edges)
  # duplicated_edges <- duplicated_edges[is.na(duplicated_edges$sign),2]
  # selectEdges(duplicated_edges, network = paste0(i, "_EGF_EGFR_RAS_ERK_pathway_start_end"))
  # deleteSelectedEdges(network = paste0(i, "_EGF_EGFR_RAS_ERK_pathway_start_end"))

  activity_weights = rbind(activity_weights,getTableColumns('node','Activity_Scores', network = paste0(i, "_EGF_EGFR_RAS_ERK_pathway_start_end")))
  cosmos_weights = rbind(cosmos_weights,getTableColumns('node','Activity', network = paste0(i, "_EGF_EGFR_RAS_ERK_pathway_start_end")))
}  

setNodeShapeDefault("ROUND_RECTANGLE", style.name) #remember to specify your style.name!
setNodeShapeMapping(table.column = "NodeType", table.column.values = list("1","2","3","4","5"),shapes = list("Octagon","VEE","Triangle","Diamond","Ellipse"), style.name = style.name, network = paste0(i, "_EGF_EGFR_RAS_ERK_pathway_start_end"))
setNodeColorDefault("#89D0F5", style.name)
setNodeLabelMapping("name", style.name, network = paste0(i, "_EGF_EGFR_RAS_ERK_pathway_start_end"))
setNodeHeightDefault(35,style.name)
setNodeFontSizeDefault(12,style.name)
setNodeFontFaceDefault("TimesNewRomanPSMT",style.name)
setNodeWidthDefault(75,style.name)
setNodeBorderWidthDefault(12, style.name = style.name)

activity_weights_min = min(activity_weights[,1],na.rm=TRUE)
activity_weights_max = max(activity_weights[,1],na.rm=TRUE)
activity_weights.values = c(activity_weights_min,0,activity_weights_max)
border.colors <- c(rev(brewer.pal(length(activity_weights.values), "RdBu")))
setNodeBorderColorMapping('Activity_Scores', activity_weights.values, border.colors, style.name = style.name, network = paste0(i, "_EGF_EGFR_RAS_ERK_pathway_start_end"))

cosmos_weights_min = min(cosmos_weights[,1],na.rm=TRUE)
cosmos_weights_max = max(cosmos_weights[,1],na.rm=TRUE)
cosmos_data.values = c(cosmos_weights_min,0,cosmos_weights_max)
node.colors <- c(rev(brewer.pal(length(cosmos_data.values), "RdBu")))
setNodeColorMapping('Activity',cosmos_data.values, node.colors, style.name = style.name, network = paste0(i, "_EGF_EGFR_RAS_ERK_pathway_start_end"))
```

``` r
for (cell_line in max_cell_lines) {
  exportImage(filename = paste0("results/cytoscape/subnetwork_",cell_line,"_analysis.pdf"), type = "PDF", resolution = 900, network = paste0(cell_line, "_EGF_EGFR_RAS_ERK_pathway_start_end"))
}

saveSession(filename = "results/cytoscape/cell_line_specific_cosmos_network_analysis.cys")
```

## Session info

``` r
sessionInfo()
```

    ## R version 4.2.0 (2022-04-22)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Monterey 12.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] cosmosR_1.5.2      RColorBrewer_1.1-3 RCy3_2.16.0        ggplot2_3.4.0     
    ##  [5] MOFA2_1.6.0        reshape2_1.4.4     dplyr_1.1.1        pheatmap_1.0.12   
    ##  [9] readr_2.1.4        moon_0.1.0         liana_0.1.5        decoupleR_2.5.2   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] uuid_1.1-0                  readxl_1.4.2               
    ##   [3] backports_1.4.1             circlize_0.4.15            
    ##   [5] corrplot_0.92               plyr_1.8.8                 
    ##   [7] igraph_1.4.2                repr_1.1.6                 
    ##   [9] sp_1.6-0                    BiocParallel_1.30.4        
    ##  [11] listenv_0.9.0               GenomeInfoDb_1.32.4        
    ##  [13] digest_0.6.31               foreach_1.5.2              
    ##  [15] htmltools_0.5.5             OmnipathR_3.9.6            
    ##  [17] fansi_1.0.4                 magrittr_2.0.3             
    ##  [19] checkmate_2.1.0             base64url_1.4              
    ##  [21] ScaledMatrix_1.4.1          cluster_2.1.4              
    ##  [23] doParallel_1.0.17           tzdb_0.3.0                 
    ##  [25] limma_3.52.4                ComplexHeatmap_2.12.1      
    ##  [27] globals_0.16.2              matrixStats_0.63.0         
    ##  [29] vroom_1.6.1                 prettyunits_1.1.1          
    ##  [31] colorspace_2.1-0            ggrepel_0.9.2              
    ##  [33] rvest_1.0.3                 rappdirs_0.3.3             
    ##  [35] xfun_0.38                   crayon_1.5.2               
    ##  [37] RCurl_1.98-1.10             jsonlite_1.8.4             
    ##  [39] graph_1.74.0                progressr_0.13.0           
    ##  [41] iterators_1.0.14            glue_1.6.2                 
    ##  [43] gtable_0.3.1                zlibbioc_1.42.0            
    ##  [45] XVector_0.36.0              GetoptLong_1.0.5           
    ##  [47] DelayedArray_0.22.0         BiocSingular_1.12.0        
    ##  [49] Rhdf5lib_1.18.2             future.apply_1.10.0        
    ##  [51] shape_1.4.6                 SingleCellExperiment_1.18.1
    ##  [53] BiocGenerics_0.42.0         HDF5Array_1.24.2           
    ##  [55] scales_1.2.1                edgeR_3.38.4               
    ##  [57] Rcpp_1.0.10                 progress_1.2.2             
    ##  [59] clue_0.3-63                 reticulate_1.28            
    ##  [61] dqrng_0.3.0                 bit_4.0.5                  
    ##  [63] rsvd_1.0.5                  stats4_4.2.0               
    ##  [65] metapod_1.4.0               httr_1.4.5                 
    ##  [67] dir.expiry_1.4.0            farver_2.1.1               
    ##  [69] XML_3.99-0.13               pkgconfig_2.0.3            
    ##  [71] scuttle_1.6.3               uwot_0.1.14                
    ##  [73] RJSONIO_1.3-1.8             locfit_1.5-9.7             
    ##  [75] utf8_1.2.3                  labeling_0.4.2             
    ##  [77] tidyselect_1.2.0            rlang_1.1.0                
    ##  [79] later_1.3.0                 munsell_0.5.0              
    ##  [81] cellranger_1.1.0            tools_4.2.0                
    ##  [83] cli_3.6.1                   generics_0.1.3             
    ##  [85] evaluate_0.20               stringr_1.5.0              
    ##  [87] fastmap_1.1.1               yaml_2.3.7                 
    ##  [89] bit64_4.0.5                 fs_1.6.1                   
    ##  [91] knitr_1.42                  purrr_1.0.1                
    ##  [93] future_1.30.0               sparseMatrixStats_1.8.0    
    ##  [95] scran_1.24.1                xml2_1.3.3                 
    ##  [97] compiler_4.2.0              rstudioapi_0.14            
    ##  [99] filelock_1.0.2              curl_5.0.0                 
    ## [101] png_0.1-8                   tibble_3.2.1               
    ## [103] statmod_1.5.0               stringi_1.7.12             
    ## [105] highr_0.10                  basilisk.utils_1.8.0       
    ## [107] forcats_1.0.0               logger_0.2.2               
    ## [109] IRdisplay_1.1               lattice_0.20-45            
    ## [111] bluster_1.6.0               Matrix_1.5-3               
    ## [113] vctrs_0.6.1                 pillar_1.9.0               
    ## [115] lifecycle_1.0.3             rhdf5filters_1.8.0         
    ## [117] GlobalOptions_0.1.2         BiocNeighbors_1.14.0       
    ## [119] cowplot_1.1.1               bitops_1.0-7               
    ## [121] irlba_2.3.5.1               GenomicRanges_1.48.0       
    ## [123] R6_2.5.1                    IRanges_2.30.1             
    ## [125] parallelly_1.34.0           codetools_0.2-18           
    ## [127] rhdf5_2.40.0                SummarizedExperiment_1.26.1
    ## [129] rjson_0.2.21                withr_2.5.0                
    ## [131] SeuratObject_4.1.3          uchardet_1.1.1             
    ## [133] S4Vectors_0.34.0            GenomeInfoDbData_1.2.8     
    ## [135] parallel_4.2.0              hms_1.1.3                  
    ## [137] IRkernel_1.3.2              grid_4.2.0                 
    ## [139] beachmat_2.12.0             tidyr_1.3.0                
    ## [141] basilisk_1.8.1              rmarkdown_2.21             
    ## [143] DelayedMatrixStats_1.18.2   MatrixGenerics_1.8.1       
    ## [145] Rtsne_0.16                  pbdZMQ_0.3-9               
    ## [147] base64enc_0.1-3             Biobase_2.56.0
