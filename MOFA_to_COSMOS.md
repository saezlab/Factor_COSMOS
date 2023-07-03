## MOFA-COSMOS pipeline

### Installation and dependency

COSMOS is dependent on CARNIVAL for exhibiting the signalling pathway
optimization. CARNIVAL requires the interactive version of IBM Cplex or
CBC-COIN solver as the network optimizer. The IBM ILOG Cplex is freely
available through Academic Initiative
[here](https://community.ibm.com/community/user/datascience/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students).
As an alternative, the CBC solver is open source and freely available
for any user, but has a significantly lower performance than CPLEX. The
CBC executable can be find under cbc/. Alternatively for small networks,
users can rely on the freely available lpSolve R-package, which is
automatically installed with the package.

In this tutorial we use *CPLEX* which is strongly recommended.

``` r
# Install cosmosR from bioconductor (stable version)
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("cosmosR")

# Or install cosmosR from github directly to obtain the newest version
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("cosmosR", quietly = TRUE)) devtools::install_github("saezlab/cosmosR")
```

We are using MOFA2 (Argelaguet et al., 2018) to find correlation between
different omics and use its output (factors) to restrict COSMOS input.
However, since we are using the python version of MOFA2, please make
sure to have [mofapy2](https://github.com/bioFAM/mofapy2) as well as
panda and numpy installed in your (conda) environment (e.g. by using
reticulate: reticulate::use_condaenv(“base”, required = T) %\>%
reticulate::conda_install(c(“mofapy2”,“panda”,“numpy”)). For downstream
analysis, the R version of MOFA2 is used.

``` r
# Install MOFA2 from bioconductor (stable version)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("MOFA2", quietly = TRUE)) BiocManager::install("MOFA2")
```

General information (tutorials) on how to use COSMOS and MOFA2, see
[COSMOS](https://saezlab.github.io/cosmosR/articles/tutorial.html) and
[MOFA2](https://biofam.github.io/MOFA2/tutorials.html).

For data manipulation, visualization, etc. further packages are needed
which are loaded here along with MOFA2 and COSMOS.

``` r
library(cosmosR)
library(MOFA2)
library(readr)
library(ggplot2)
library(ggfortify)
library(dplyr)
library(reshape2)
library(liana)
library(decoupleR)
library(moon)
library(pheatmap)
library(gridExtra)
library(liana)
library(GSEABase)
library(tidyr)
library(RCy3)
```

### From omics data to MOFA ready input

Since extensive pre-processing is beyond the scope of this tutorial,
here only the transformation of pre-processed single omics data to a
long data.frame (columns: “sample”, “group”, “feature”, “view”, “value”)
is shown. For more details and information regarding appropriate input,
please refer to the MOFA tutorial [MOFA2: training a model in
R](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/getting_started_R.html).
The raw data is accessed through [NCI60
cellminer](https://discover.nci.nih.gov/cellminer/home.do) and the
scripts for pre-processing can be found in the associated folder.

Here, we have the following views (with samples as columns and genes as
rows): transcriptomics (RNA-seq
[<PMID:31113817>](https://pubmed.ncbi.nlm.nih.gov/31113817/)),
proteomics (SWATH (Mass spectrometry)
[<PMID:31733513>](https://pubmed.ncbi.nlm.nih.gov/31733513/%3E)} and
metabolomics (LC/MS & GC/MS (Mass spectrometry) [DTP NCI60
data](https://wiki.nci.nih.gov/display/NCIDTPdata/Molecular+Target+Data)).
The group information is stored inside the RNA metadata and based on
transcription factor clustering (see scripts/RNA/Analyze_RNA.R).
Further, only samples are kept for which each omic is available (however
this is not necessary since MOFA2 is able to impute data).

``` r
# Transcriptomics
## Load data
RNA <- as.data.frame(read_csv("data/RNA/RNA_log2_FPKM_clean.csv"))
rownames(RNA) <- RNA[,1]
RNA <- RNA[,-1]

## Remove genes with excessive amount of NAs (only keep genes with max. amount of NAs = 33.3 % across cell lines)
RNA <- RNA[rowSums(is.na(RNA))<(dim(RNA)[2]/3),]

## Only keep highly variable genes (as suggested by MOFA): here top ≈70% genes to keep >75% of TFs (103 out of 133). For stronger selection, we can further reduce the number of genes (e.g. top 50%)
RNA_sd <- sort(apply(RNA, 1, function(x) sd(x,na.rm = T)), decreasing = T)
hist(RNA_sd, breaks = 100)
```

![](MOFA_to_COSMOS_files/figure-gfm/Preparing%20MOFA%20input-1.png)<!-- -->

``` r
hist(RNA_sd[1:6000], breaks = 100)
```

![](MOFA_to_COSMOS_files/figure-gfm/Preparing%20MOFA%20input-2.png)<!-- -->

``` r
RNA_sd <- RNA_sd[1:6000]
RNA <- RNA[names(RNA_sd),]

# Proteomics
## Load data
proteo <- as.data.frame(read_csv("data/proteomic/Prot_log10_SWATH_clean.csv"))
rownames(proteo) <- proteo[,1]
proteo <- proteo[,-1]

## Only keep highly variable proteins (as suggested by MOFA): here top ≈60%. For stronger selection, we can further reduce the number of proteins (e.g. only keep top 50%)
proteo_sd <- sort(apply(proteo, 1, function(x) sd(x,na.rm = T)), decreasing = T)
hist(proteo_sd, breaks = 100)
```

![](MOFA_to_COSMOS_files/figure-gfm/Preparing%20MOFA%20input-3.png)<!-- -->

``` r
hist(proteo_sd[1:round(dim(proteo)[1]*3/5)], breaks = 100)
```

![](MOFA_to_COSMOS_files/figure-gfm/Preparing%20MOFA%20input-4.png)<!-- -->

``` r
proteo_sd <- proteo_sd[1:round(dim(proteo)[1]*3/5)]
proteo <- proteo[names(proteo_sd),]

# Metabolomics
## Load data
metab <- as.data.frame(read_csv("data/metabolomic/metabolomic_clean_vsn.csv"))
rownames(metab) <- metab[,1]
metab <- metab[,-1]

## Since we have only a limited number of metabolite measurements, all measurements are kept here
metab_sd <- apply(metab, 1, function(x) sd(x,na.rm=T))
hist(metab_sd, breaks = 100)
```

![](MOFA_to_COSMOS_files/figure-gfm/Preparing%20MOFA%20input-5.png)<!-- -->

``` r
# Create long data frame
## Only keep samples with each view present 
overlap_patients <- intersect(intersect(names(RNA),names(proteo)),names(metab))
RNA <- RNA[,overlap_patients]
proteo <- proteo[,overlap_patients]
metab <- metab[,overlap_patients]

## Create columns required for MOFA
### RNA
RNA <- melt(as.data.frame(cbind(RNA,row.names(RNA))))
RNA$view <- "RNA"
RNA <- RNA[,c(2,1,4,3)]                 
names(RNA) <- c("sample","feature","view","value")

### Proteomics
proteo <- melt(as.data.frame(cbind(proteo,row.names(proteo))))
proteo$view <- "proteo"
proteo <- proteo[,c(2,1,4,3)]                 
names(proteo) <- c("sample","feature","view","value")

### Metabolomics
metab <- melt(as.data.frame(cbind(metab,row.names(metab))))
metab$view <- "metab"
metab <- metab[,c(2,1,4,3)]                 
names(metab) <- c("sample","feature","view","value")

## Merge long data frame to one
mofa_ready_data <- as.data.frame(do.call(rbind,list(RNA,proteo,metab)))

## Add metadata information (cluster assignment)
meta_data <- read_csv("data/metadata/RNA_metadata_cluster.csv")[,c(1,2)]
colnames(meta_data) <- c("sample","cluster")
mofa_ready_data <- merge(mofa_ready_data, meta_data, by = "sample")
mofa_ready_data <- mofa_ready_data[,c(1,5,2,3,4)]
names(mofa_ready_data) <- c("sample","group","feature","view","value")

# Rename clusters
mofa_ready_data[grepl("1",mofa_ready_data$group),2] <- "cluster_1"
mofa_ready_data[grepl("2",mofa_ready_data$group),2] <- "cluster_2"
mofa_ready_data[grepl("3",mofa_ready_data$group),2] <- "cluster_3"

# Optional: Only keep cluster 1 and 3
#mofa_ready_data <- mofa_ready_data[which(mofa_ready_data$group %in% c("1","3")),]

# Here: Remove group column for unsupervised analysis
mofa_ready_data <- mofa_ready_data[,-2]

# Save MOFA metadata information and MOFA input
write_csv(mofa_ready_data, file = "data/mofa/mofa_data.csv")
```

The final data frame has the following structure:

``` r
head(mofa_ready_data)
```

    ##   sample  feature view  value
    ## 1  786-0     KRT8  RNA  7.261
    ## 2  786-0      FN1  RNA  7.823
    ## 3  786-0   S100A4  RNA     NA
    ## 4  786-0 RNA5-8S5  RNA 11.045
    ## 5  786-0     MT1E  RNA  7.010
    ## 6  786-0      VIM  RNA  9.366

We have successfully combined single omic data sets to a long
multi-omics data frame.

### MOFA run

After pre-processing the data into MOFA appropriate input, let’s perform
the MOFA analysis. To change specific MOFA options the function
write_MOFA_options_file is used. For more information regarding option
setting see [MOFA option
tutorial](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/getting_started_R.html).
Since we don’t know how many factors are needed to explain our data
well, various maximum numbers of factors are tested. The final results
can be found in the results/mofa/ folder.

Do not run this chunk, unless you want to repeat the mofa optimisation,
because it takes a lot of time ot run. The results are already available
in the results/mofa folder.

``` r
for(i in c(5:15,20,length(unique(mofa_ready_data$sample)))){
wd <- getwd()
write_MOFA_options_file(input_file = '/data/mofa/mofa_data.csv', output_file = paste0('/results/mofa/mofa_res_',i,'factor.hdf5'), factors = i, convergence_mode = "slow", likelihoods = c('gaussian','gaussian','gaussian'))
system(paste0(paste('python3', wd, sep = " "), paste('/scripts/mofa/mofa2.py', "options_MOFA.csv")))
}
```

### Investigate MOFA output

We first load the different MOFA models containing the MOFA results
together with the metadata information.

``` r
# Load MOFA output
for(i in c(5:15,20,length(unique(mofa_ready_data$sample)))){
  filepath <- paste0('results/mofa/mofa_res_',i,'factor.hdf5')
  model <- load_model(filepath)
  meta_data <- read_csv("data/metadata/RNA_metadata_cluster.csv")[,c(1,2)]
  colnames(meta_data) <- c("sample","cluster")
  mofa_metadata <- merge(samples_metadata(model), meta_data, by = "sample")
  names(mofa_metadata) <- c("sample","group","cluster")
  samples_metadata(model) <- mofa_metadata
  assign(paste0("model",i),model)
}
```

Then, we can compare the factors from different MOFA models and check on
consistency (here: model 7-13).

``` r
compare_factors(list(model7,model8,model9,model10,model11,model12,model13), cluster_rows = F, cluster_cols = F,)
```

![](MOFA_to_COSMOS_files/figure-gfm/Compare%20MOFA%20models-1.png)<!-- -->

Here, we can see that the inferred MOFA factor weights do not change
drastically around 10. For the sake of this tutorial, we choose the
10-factor MOFA model.

``` r
model <- model10
```

The next part deals with the analysis and interpretation of the MOFA
factors with specific focus on the factor weights. A more detailed
tutorial on how to perform downstream analysis after MOFA model training
can be found here [MOFA+: downstream analysis (in
R)](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/downstream_analysis.html).
In this context, it is important to emphasize that the factors can be
interpreted similarly to principal components (PCs) in a principal
component analysis (PCA).

An initial metric we can explore is the total variance ($R^2$) explained
per view by our MOFA model. This helps us understand how well the model
represents the data.

``` r
# Investigate factors: Explained total variance per view
plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
```

![](MOFA_to_COSMOS_files/figure-gfm/Total%20variance%20per%20factor-1.png)<!-- -->

The bar plot demonstrates that the nine factors for the view RNA explain
more than 50% of the variance. In contrast, for both metabolomics as
well as proteomics around 20% of the variance can be explained by all
factors.

Since we would like to use the found correlation for downstream COSMOS
analysis, we first investigate the variance ($R^2$) each factor can
explain per view as well as investigate the factor explaining the most
variance across all views.

``` r
# Investigate factors: Explained variance per view for each factor 
pheatmap(model@cache$variance_explained$r2_per_factor[[1]], display_numbers = T, angle_col = "0", legend_labels = c("0","10", "20", "30", "40", "Variance\n\n"), legend = T, main = "", legend_breaks = c(0,10, 20, 30, 40, max(model@cache$variance_explained$r2_per_factor[[1]])), cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white","red"))(100), fontsize_number = 10)
```

![](MOFA_to_COSMOS_files/figure-gfm/Variance%20per%20view%20per%20factor-1.png)<!-- -->

By analyzing this heatmap, we can observe that Factor 1 and 3 explain a
large proportion of variance of the RNA, Factor 2 mainly explains the
variance of the metabolomics and Factor 3 highlights variability from
the proteomics view. Consequently, no factor explains a large portion of
the variance across all views. To potentially justify to use different
factor weights downstream for different views, we can analyze the
correlation between factors. Thus, the next plot highlights the spearman
correlation between each factor.

``` r
# Investigate factors: Correlation between factors 
plot_factor_cor(model, type = 'upper', method = "spearman", addCoef.col = "grey")
```

![](MOFA_to_COSMOS_files/figure-gfm/Correlation%20between%20factors-1.png)<!-- -->

Since the correlation between each factor is relatively low, we have to
choose a single factor. In this case, we are going to use Factor 4,
since a decent amount of variability from each view is explained.

Additionally, we can use a beeswarm plot to examine each sample’s factor
value and to see if we are able to separate the samples based on their
transcription factor clustering.

``` r
# Investigate factors: Beeswarm plots of individual factors
plot_factor(model, 
                 factors = 'all',
                 color_by = "cluster",
                 dot_size = 3,
                 dodge = T,           
                 legend = T,          
                 add_violin = T,     
                 violin_alpha = 0.25) + 
   scale_color_manual(values=c("1"="red", "2"="cyan", "3"="orange")) +
   scale_fill_manual(values=c("1"="red", "2"="cyan", "3"="orange"))
```

    ## Warning: `fct_explicit_na()` was deprecated in forcats 1.0.0.
    ## ℹ Please use `fct_na_value_to_level()` instead.
    ## ℹ The deprecated feature was likely used in the MOFA2 package.
    ##   Please report the issue at <https://github.com/bioFAM/MOFA2>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](MOFA_to_COSMOS_files/figure-gfm/Factor%20value%20per%20sample%20and%20factor-1.png)<!-- -->

``` r
plot_factor(model, 
                 factors = 4,
                 color_by = "cluster",
                 dot_size = 3,
                 dodge = T,           
                 legend = T,          
                 add_violin = T,     
                 violin_alpha = 0.25) + 
   scale_color_manual(values=c("1"="red", "2"="cyan", "3"="orange")) +
   scale_fill_manual(values=c("1"="red", "2"="cyan", "3"="orange"))
```

![](MOFA_to_COSMOS_files/figure-gfm/Factor%20value%20per%20sample%20and%20factor-2.png)<!-- -->

``` r
plot_factor(model, 
                 factors = 'all',
                 dot_size = 3,
                 dodge = T,           
                 legend = T,          
                 add_violin = T,     
                 violin_alpha = 0.25)
```

![](MOFA_to_COSMOS_files/figure-gfm/Factor%20value%20per%20sample%20and%20factor-3.png)<!-- -->

Interestingly, with Factor 4, we have a separation of the samples
belonging to cluster 2 from samples of cluster 1 and 3. Taking into
account that this factor mainly explains the variation in the RNA view,
the samples are classifiable into the clusters especially by
transcriptomics. This is expected since the assignment of the samples
was done based on transcription factor clustering via transcriptomics.
Thus, it was shown that the original structure of the data was conserved
in the low-dimensional representation via factors.

Further, we can inspect the composition of the factors by plotting the
top 10 weights per factor and view.

``` r
## Top 20 factor weights per view
grid.arrange(
  ### RNA
  plot_top_weights(model, 
                   view = "RNA",
                   factor = 4,
                   nfeatures = 10) +
    ggtitle("RNA"),
  ### Proteomics
  plot_top_weights(model, 
                   view = "proteo",
                   factor = 4,
                   nfeatures = 10) + 
    ggtitle("Proteomics"),
  ### Metabolomics
  plot_top_weights(model, 
                   view = "metab",
                   factor = 4,
                   nfeatures = 10) + 
    ggtitle("Metabolomics")
)
```

![](MOFA_to_COSMOS_files/figure-gfm/Factor%20weights%20per%20view-1.png)<!-- -->

Using this visualization, features with strong association with the
factor (large absolute values) can be easily identified and further
inspected by literature investigation.

Moreover, we can use heatmaps to highlight the coordinated heterogeneity
that MOFA captures in the original data and define clusters by
hierarchical clustering. Here, for example, the heatmap for Factor 4 of
the RNA view using the top 25 feature is shown and hierarchical
clustering with complete linkage is performed. Additionally, the cluster
assignment of each sample is depicted.

``` r
anno_colors <- list(
  "1"="red", "2"="cyan", "3"="orange"
)

model@samples_metadata$cluster <- as.character(model@samples_metadata$cluster)
### RNA
plot_data_heatmap(model,
                  view = "RNA",       
                  factor = 4,             
                  features =25,          
                  cluster_rows = T, cluster_cols = T,
                  show_rownames = T, show_colnames = F,
                  annotation_samples = "cluster",
                  annotation_colors = anno_colors)
```

![](MOFA_to_COSMOS_files/figure-gfm/Heatmap%20RNA%20view%20factor%202-1.png)<!-- -->

Further tools to investigate the latent factors are available inside the
MOFA framework ([MOFA+: downstream analysis (in
R)](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/downstream_analysis.html)).

### From MOFA output to TF activity, ligand-receptor activity

In the next step, using the weights of Factor 4 (highest shared $R^2$
across views), different databases (liana and dorothea) and decoupleR,
the transcription factor activities, ligand-receptor scores are
estimated.

First, we have to extract the weights from our MOFA model and keep the
weights from Factor 4 in the RNA view.

``` r
weights <- get_weights(model, views = "all", factors = "all")

# Extract Factor with highest explained variance in RNA view
RNA <- data.frame(weights$RNA[,4]) 
row.names(RNA) <- gsub("_RNA","",row.names(RNA))
```

Then we load the consensus networks from our databases.

``` r
# Load LIANA (receptor and ligand) consensus network
# ligrec_ressource <- distinct(liana::decomplexify(liana::select_resource("Consensus")[[1]]))
# save(ligrec_ressource, file = "support/ligrec_ressource.RData")
load(file = "support/ligrec_ressource.RData")
ligrec_geneset <- ligrec_ressource[,c("source_genesymbol","target_genesymbol")]
ligrec_geneset$set <- paste(ligrec_geneset$source_genesymbol, ligrec_geneset$target_genesymbol, sep = "___")
ligrec_geneset <- reshape2::melt(ligrec_geneset, id.vars = "set")[,c(3,1)]
names(ligrec_geneset)[1] <- "gene"
ligrec_geneset$mor <- 1
ligrec_geneset$likelihood <- 1
ligrec_geneset <- distinct(ligrec_geneset)

# Load Dorothea (TF) network
# dorothea_df <- decoupleR::get_collectri()
# save(dorothea_df, file = "support/dorothea_df.RData")
load(file = "support/dorothea_df.RData")
```

Then, by using decoupleR and prior knowledge networks, we calculate the
different regulatory activities inferred by the normalized weighted mean
approach. Depending on the task, the minimum number of targets per
source varies required to maintain the source in the output (e.g. two
targets by definition for ligand-receptor interactions). To calculate
activities of master regulons through moon, transcription factor
activities estimated by decoupleR and dorothea are used as an input.

``` r
# Calculate regulatory activities from Receptor and Ligand network
ligrec_high_vs_low <- run_ulm(mat = as.matrix(RNA), network = ligrec_geneset, .source = set, .target = gene, minsize = 2) 
ligrec_high_vs_low <- ligrec_high_vs_low[ligrec_high_vs_low$statistic == "ulm",]
ligrec_high_vs_low_vector <- ligrec_high_vs_low$score
names(ligrec_high_vs_low_vector) <- ligrec_high_vs_low$source

ligrec_high_vs_low_top <- ligrec_high_vs_low[order(abs(ligrec_high_vs_low$score),decreasing = T)[1:15],c(2,4)]
ligrec_high_vs_low_top$source <- factor(ligrec_high_vs_low_top$source, levels = ligrec_high_vs_low_top$source)
ggplot(ligrec_high_vs_low_top, aes(x=source, y = score)) + geom_bar(stat= "identity",position = "dodge") + theme_minimal() + theme(axis.text.x = element_text(angle = 315, vjust = 0.5, hjust=0))
```

![](MOFA_to_COSMOS_files/figure-gfm/Activity%20estimations-1.png)<!-- -->

``` r
# Calculate regulatory activities from TF network
TF_high_vs_low <- run_ulm(mat = as.matrix(RNA), network = dorothea_df, minsize = 10)
TF_high_vs_low <- TF_high_vs_low[TF_high_vs_low$statistic == "ulm",]
TF_high_vs_low_vector <- TF_high_vs_low$score
names(TF_high_vs_low_vector) <- TF_high_vs_low$source

# Prepare moon analysis input
TF_high_vs_low <- as.data.frame(TF_high_vs_low[,c(2,4)])
row.names(TF_high_vs_low) <- TF_high_vs_low[,1]

TF_high_vs_low_top_10 <- TF_high_vs_low[order(abs(TF_high_vs_low$score), decreasing = T)[1:10],]
TF_high_vs_low_top_10$source <- factor(TF_high_vs_low_top_10$source, levels = TF_high_vs_low_top_10$source)
ggplot(TF_high_vs_low_top_10, aes(x=source, y = score)) + geom_bar(stat= "identity",position = "dodge") + theme_minimal()
```

![](MOFA_to_COSMOS_files/figure-gfm/Activity%20estimations-2.png)<!-- -->

``` r
# Combine results
ligrec_TF_moon_inputs <- list("ligrec" = ligrec_high_vs_low_vector,
                              "TF" = TF_high_vs_low_vector)

save(ligrec_TF_moon_inputs, file = "data/cosmos/ligrec_TF_moon_inputs.Rdata")
```

We have now successfully used the MOFA weights (Factor 4, RNA view) to
infer not only TF activities but also the ligand-receptor interactions.

### Extracting network with COSMOS to formulate mechanistic hypotheses

In this part the COSMOS run is performed. We can use the output of
decoupleR (the activities) next to the factor weights as an input for
COSMOS to specifically focus the analysis on pre-computed data-driven
results.

Here, in the first step, the data is loaded. Besides the activity
estimations and the factor weights, the previous filtered out expressed
gene names are loaded.

``` r
# Load data
load(file = "data/cosmos/ligrec_TF_moon_inputs.Rdata")

signaling_input <- ligrec_TF_moon_inputs$TF
ligrec_input <- ligrec_TF_moon_inputs$ligrec

RNA_input <- weights$RNA[,4] 
prot_input <- weights$proteo[,4]
metab_inputs <- weights$metab[,4]

names(RNA_input) <- gsub("_RNA","",names(RNA_input))
names(prot_input) <- gsub("_proteo","",names(prot_input))

expressed_genes <- as.data.frame(read_csv("data/RNA/RNA_log2_FPKM_clean.csv"))$Genes #since only complete cases were considered in the first part, here the full gene list is loaded
```

Next, we perform filtering steps in order to identify top deregulated
features. Here, the individual filtering is based on absolute factor
weight or activity score. The first step involves the visualization of
the respective distribution to define an appropriate threshold.

First, the RNA factor weights are filtered based on a threshold defined
by their distribution as well a threshold defined by the distribution of
the protein factor weights.

``` r
##RNA
{plot(density(RNA_input))
abline(v = -0.2)
abline(v = 0.2)}
```

![](MOFA_to_COSMOS_files/figure-gfm/Plot:%20RNA%20weights%20&%20protein%20weights-1.png)<!-- -->

``` r
{plot(density(prot_input))
abline(v = -0.05)
abline(v = 0.05)}
```

![](MOFA_to_COSMOS_files/figure-gfm/Plot:%20RNA%20weights%20&%20protein%20weights-2.png)<!-- -->

Based on the plots and indicated by the straight lines, the RNA weights
threshold is set to -0.2 and 0.2, and the protein weights threshold to
-0.05 and 0.05. RNA factor weight values lying inside this threshold are
set to 0 and previously filtered out genes (before MOFA analysis) are
also included with value 0.

``` r
for(gene in names(RNA_input)) {
  if (RNA_input[gene] > -0.2 & RNA_input[gene] < 0.2)
  {
    RNA_input[gene] <- 0
  } else
  {
    RNA_input[gene] <- sign(RNA_input[gene]) * 10
  }
  if (gene %in% names(prot_input))
  {
    if (prot_input[gene] > -0.05 & prot_input[gene] < 0.05)
    {
      RNA_input[gene] <- 0
    } else
    {
      RNA_input[gene] <- sign(RNA_input[gene]) * 10
    }
  }
}

expressed_genes <- expressed_genes[!(expressed_genes %in% names(RNA_input))]
genes <- expressed_genes
expressed_genes <- rep(0,length(expressed_genes))
names(expressed_genes) <- genes

RNA_input <- c(RNA_input, expressed_genes)
```

The same procedure is repeated for the transcription factor activities.

``` r
## Transcription factors
{plot(density(signaling_input))
abline(v = -0.5)
abline(v = 3.5)}
```

![](MOFA_to_COSMOS_files/figure-gfm/Plot:%20TF%20activity%20estimates-1.png)<!-- -->

The threshold is set to -0.5 and 3.5 respectively. Further, the
TF_to_remove variable is later used to remove transcription factors with
activities below this threshold from the prior knowledge network.

``` r
TF_to_remove <- signaling_input[signaling_input > -0.5 & signaling_input < 3.5]
signaling_input <- signaling_input[signaling_input < -0.5 | signaling_input > 3.5]
```

The same procedure is repeated for the ligand-receptor activities.

``` r
## Ligand-receptor interactions
{plot(density(ligrec_input))
abline(v = -0.5)
abline(v = 2.5)}
```

![](MOFA_to_COSMOS_files/figure-gfm/Plot:%20Ligand-receptor%20activity%20estimates-1.png)<!-- -->

Here, only activities higher than 2.5 or lower than -0.5 are kept (see
straight lines in plot). Since at this point we are interested in the
receptors, the names are adjusted accordingly and if there are multiple
mentions of a receptor, its activity is calculated by the mean of the
associated ligand-receptor interactions.

``` r
ligrec_input <- ligrec_input[ligrec_input > 2.5 | ligrec_input < -0.5]

rec_inputs <- ligrec_input
names(rec_inputs) <- gsub(".+___","",names(rec_inputs))
rec_inputs <- tapply(rec_inputs, names(rec_inputs), mean)
receptors <- names(rec_inputs)
rec_inputs <- as.numeric(rec_inputs)
names(rec_inputs) <- receptors
```

Finally, the same procedure is repeated for the metabolite factor
weights.

``` r
## Metabolites
{plot(density(metab_inputs))
abline(v = -0.2)
abline(v = 0.2)}
```

![](MOFA_to_COSMOS_files/figure-gfm/Plot:%20Metabolite%20factor%20weights-1.png)<!-- -->

The threshold is set to 0.2 and -0.2 respectively. Further, the
metab_to_exclude variable is later used to remove metabolites with
activities below this threshold from the prior knowledge network.
Moreover, metabolite names are translated into HMDB format and metabolic
compartment codes (m = mitochondria, c = cytosol) as well as metab\_\_
prefix is added to HMDB IDs. More information regarding compartment
codes can be found [here](http://bigg.ucsd.edu/compartments).

``` r
metab_to_HMDB <- as.data.frame(
  read_csv("data/metabolomic/MetaboliteToHMDB.csv"))
metab_to_HMDB <- metab_to_HMDB[metab_to_HMDB$common %in% names(metab_inputs),]
metab_inputs <- metab_inputs[metab_to_HMDB$common]
names(metab_inputs) <- metab_to_HMDB$HMDB

metab_inputs <- cosmosR::prepare_metab_inputs(metab_inputs, compartment_codes = c("m","c"))

metab_to_exclude <- metab_inputs[abs(metab_inputs) < 0.2]
metab_inputs <- metab_inputs[abs(metab_inputs) > 0.2]
```

Next, the prior knowledge network (PKN) is loaded. To see how the full
meta PKN was assembled, see
[PKN](https://github.com/saezlab/meta_PKN_BIGG.git). Since we are not
interested in including transcription factors and metabolites that were
previously filtered out, we also remove the respective nodes from the
PKN.

``` r
## Load and filter meta network
data("meta_network")

meta_network <- meta_network_cleanup(meta_network)

meta_network <- meta_network[!(meta_network$source %in% names(TF_to_remove)) & !(meta_network$target %in% names(TF_to_remove)),]
meta_network <- meta_network[!(meta_network$source %in% names(metab_to_exclude)) & !(meta_network$target %in% names(metab_to_exclude)),]
```

Then the datasets are merged, entries that are not included in the PKN
are removed and the filtering steps are completed. The merging in this
case is based on the idea of deriving the forward network via the
receptor to transcription factor/metabolite direction.

``` r
upstream_inputs <- c(rec_inputs)
downstream_inputs <- c(metab_inputs*10, signaling_input)

upstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(upstream_inputs, meta_network)
```

    ## [1] "COSMOS: 5 input/measured nodes are not in PKN any more: BCAM, CD63, LRP10, NOTCH1, RXRA and 0 more."

``` r
downstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(downstream_inputs, meta_network)
```

    ## [1] "COSMOS: 33 input/measured nodes are not in PKN any more: Metab__HMDB0011747_m, Metab__HMDB0000272_m, Metab__HMDB0001173_m, Metab__HMDB0001548_m, Metab__HMDB0001893_m, Metab__HMDB0001123_m and 27 more."

``` r
rec_to_TF_cosmos_input <- list(upstream_inputs_filtered, downstream_inputs_filtered)
save(rec_to_TF_cosmos_input, file = "data/cosmos/rec_to_TF_cosmos_input.Rdata")
```

At this point, we can perform the COSMOS analysis. Firstly, the options
for the CARNIVAL run such as the time limit and min gap tolerance are
set and the user must provide a path to its CPLEX executable. You can
check the CARNIVAL_options variable to see all possible options that can
be adjusted.

``` r
## Set CARNIVAL options to solve pre-optimization problem
my_options <- default_CARNIVAL_options(solver = "cplex")
my_options$solverPath <- "cplex_macos/cplex"
my_options$solver <- "cplex"
my_options$mipGAP <- 0.05
my_options$threads <- 2
my_options$timelimit <- 360
my_options$limitPop <- 100

data("HMDB_mapper_vec")
```

Here, we are trying to find the “forward” network causally connecting
the signaling data to metabolic data by using Carnival’s ILP solution
and the PKN. To further simplify the prior knowledge network (and
perform checks on the input data), the pre-processing function of COSMOS
is used. Details about the options can be found in the documentation as
well as in the [COSMOS
tutorial](https://saezlab.github.io/cosmosR/articles/tutorial.html).

``` r
pre_run_rec_to_TF <- preprocess_COSMOS_signaling_to_metabolism(meta_network = meta_network,
                                                               signaling_data = upstream_inputs_filtered,
                                                               metabolic_data = downstream_inputs_filtered,
                                                               diff_expression_data = RNA_input,
                                                               maximum_network_depth = 4,
                                                               remove_unexpressed_nodes = T,
                                                               filter_tf_gene_interaction_by_optimization = F,
                                                               CARNIVAL_options = my_options)
```

In this part, we can set up the options for the optimization run. The
running time should be much higher here than in pre-optimization. You
can increase the number of threads to use if you have many available
CPUs.

``` r
## Set CARNIVAL options to solve optimization problem
my_options$mipGAP <- 0.05
my_options$threads <- 6
my_options$timelimit <- 720
my_options$limitPop <- 100
```

After pre-optimization, we can perform the actual COSMOS run using the
pre-optimized solution.

``` r
run_rec_to_TF <- run_COSMOS_signaling_to_metabolism(data = pre_run_rec_to_TF,
                                                    CARNIVAL_options = my_options)
run_rec_to_TF <- cosmosR::format_COSMOS_res(run_rec_to_TF, metab_mapping = HMDB_mapper_vec)
save(run_rec_to_TF, file = "results/cosmos/run_rec_to_TF.Rdata")
```

In order to analyze the gained sub-network, we first process the COSMOS
result by extracting the list of interactions (nodes, sign, weight) in
the simple interaction format (SIF) and by extracting a list of
accompanying attributes (ATT). Here, we can also remove nodes with
average activity of 0 as well as with weight of 0.

``` r
## COSMOS output evaluation
load("results/cosmos/run_rec_to_TF.Rdata")

SIF_rec_to_TF <- as.data.frame(run_rec_to_TF[[1]])
SIF_rec_to_TF <- SIF_rec_to_TF[which(SIF_rec_to_TF$Weight != 0),]

ATT_rec_to_TF <- as.data.frame(run_rec_to_TF[[2]])
colnames(ATT_rec_to_TF)[1] <- "Nodes"
ATT_rec_to_TF$measured <- ifelse(ATT_rec_to_TF$NodeType %in% c("M","T","S","P"),1,0)
ATT_rec_to_TF$Activity <- ATT_rec_to_TF$AvgAct
ATT_rec_to_TF <- ATT_rec_to_TF[ATT_rec_to_TF$AvgAct != 0,]

res_mofa_rec_to_TF <- list(SIF_rec_to_TF,ATT_rec_to_TF)
save(res_mofa_rec_to_TF, file = "results/cosmos/res_mofamoon_rec_to_TFmet.RData")
write_csv(SIF_rec_to_TF, file = "results/cosmos/SIF_res_mofamoon_rec_to_TFmet.csv")
write_csv(ATT_rec_to_TF, file = "results/cosmos/ATT_res_mofamoon_rec_to_TFmet.csv")
```

Since we are also interested to derive the forward network via the
transcription factor/moon to ligand direction, this COSMOS analysis is
additionally performed.

We first assign the activities of the ligand-receptor interactions to
the ligands and calculate the mean activity of the ligands in case of
multiple entries.

``` r
lig_inputs <- ligrec_input
names(lig_inputs) <- gsub("___.+","",names(lig_inputs))
lig_inputs <- tapply(lig_inputs, names(lig_inputs), mean)
ligands <- names(lig_inputs)
lig_inputs <- as.numeric(lig_inputs)
names(lig_inputs) <- ligands
```

Then the datasets are merged appropriately, entries that are not
included in the PKN are removed and the filtering steps are completed.

``` r
upstream_inputs <- c(signaling_input)
downstream_inputs <- c(lig_inputs, metab_inputs*10) #factor 10 to scale up the metabs to the ulm scores

upstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(upstream_inputs, meta_network)
```

    ## [1] "COSMOS: 9 input/measured nodes are not in PKN any more: HR, ZNF143, TFE3, FOXN1, BMAL1, TBX15 and 3 more."

``` r
downstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(downstream_inputs, meta_network)
```

    ## [1] "COSMOS: 29 input/measured nodes are not in PKN any more: FABP5, LGALS1, LGALS3BP, LTBP3, YBX1, Metab__HMDB0011747_m and 23 more."

``` r
TF_to_lig_cosmos_input <- list(upstream_inputs_filtered, downstream_inputs_filtered)
save(TF_to_lig_cosmos_input, file = "data/cosmos/TF_to_lig_cosmos_input.Rdata")
```

``` r
## Set CARNIVAL options to solve optimization problem
my_options$mipGAP <- 0.05
my_options$threads <- 6
my_options$timelimit <- 360
my_options$limitPop <- 100
```

After setting the variables for the CPLEX solver, we can perform the
pre-optimization run …

``` r
pre_run_TF_to_lig <- preprocess_COSMOS_signaling_to_metabolism(meta_network = meta_network, 
                                                               signaling_data = upstream_inputs_filtered,
                                                               metabolic_data = downstream_inputs_filtered,
                                                               diff_expression_data = RNA_input,
                                                               maximum_network_depth = 5,
                                                               remove_unexpressed_nodes = T,
                                                               filter_tf_gene_interaction_by_optimization = F,
                                                               CARNIVAL_options = my_options)
```

… set the options for the actual run …

``` r
## Set CARNIVAL options to solve optimization problem
my_options$mipGAP <- 0.05
my_options$threads <- 6
my_options$timelimit <- 720
my_options$limitPop <- 100
```

… and perform the optimization.

``` r
run_TF_to_lig <- run_COSMOS_signaling_to_metabolism(data = pre_run_TF_to_lig,
                                                    CARNIVAL_options = my_options)

run_TF_to_lig <- cosmosR::format_COSMOS_res(run_TF_to_lig, metab_mapping = HMDB_mapper_vec)
save(run_TF_to_lig, file = "results/cosmos/run_TF_to_lig.Rdata")
```

Again, we extract the information from the determined network and filter
out inactive nodes.

``` r
load("results/cosmos/run_TF_to_lig.Rdata")


SIF_TF_to_lig <- as.data.frame(run_TF_to_lig[[1]])
SIF_TF_to_lig <- SIF_TF_to_lig[which(SIF_TF_to_lig$Weight != 0),]


ATT_TF_to_lig <- as.data.frame(run_TF_to_lig[[2]])
colnames(ATT_TF_to_lig)[1] <- "Nodes"
ATT_TF_to_lig$measured <- ifelse(ATT_TF_to_lig$NodeType %in% c("M","T","S","P"),1,0)
ATT_TF_to_lig$Activity <- ATT_TF_to_lig$AvgAct
ATT_TF_to_lig <- ATT_TF_to_lig[ATT_TF_to_lig$AvgAct != 0,]

res_mofa_TF_to_lig <- list(SIF_TF_to_lig,ATT_TF_to_lig)
save(res_mofa_TF_to_lig, file = "results/cosmos/res_mofamoon_TF_to_lig.RData")
write_csv(SIF_TF_to_lig, file = "results/cosmos/SIF_mofamoon_TF_to_lig.csv")
write_csv(ATT_TF_to_lig, file = "results/cosmos/ATT_mofamoon_TF_to_lig.csv")
```

The next step connects both networks, resulting in a complete network of
interactions between receptors, ligands, transcription factor
regulators, transcription factors, and metabolites.

After merging both inferred networks ATT file, we can identify whether
the activity of the node is positive or negative and save the result
under the variable sign_activity. Here, to avoid multiple mentions of
one specific node, we can combine multiple mentions by calculating the
mean.

``` r
combined_ATT <- as.data.frame(rbind(ATT_rec_to_TF, ATT_TF_to_lig))
combined_ATT$sign_activity <- sign(combined_ATT$AvgAct)

combined_ATT$NodeType <- ifelse(combined_ATT$NodeType == "", 0, 1) #if NodeType is available, set to 1, else set to 0
combined_ATT <- as.data.frame(combined_ATT %>% 
  group_by(Nodes) %>%
  summarise(across(.fns= ~ mean(.x, na.rm = TRUE))))
```

Further, we can also add the MOFA weight and activity score information
to the list of attributes.

``` r
## MOFA weights
MOFA_weights <- get_weights(model, factors = 4, as.data.frame = T)
MOFA_weights <- MOFA_weights %>%
  arrange(feature, view) %>%
  filter(!duplicated(feature))
MOFA_weights <- MOFA_weights[,c(1,3)]
colnames(MOFA_weights) <- c("Nodes", "mofa_weights")

#remove RNA and proteo tags and integrate values
MOFA_weights$Nodes <- gsub("_RNA$","",MOFA_weights$Nodes)
MOFA_weights$Nodes <- gsub("_proteo$","",MOFA_weights$Nodes)
MOFA_weights$sign <- sign(MOFA_weights$mofa_weights)
MOFA_weights$mofa_weights <- abs(MOFA_weights$mofa_weights)

MOFA_weights <- MOFA_weights %>% group_by(Nodes) %>% summarise_each(funs(max(., na.rm = TRUE)))
```

    ## Warning: `funs()` was deprecated in dplyr 0.8.0.
    ## ℹ Please use a list of either functions or lambdas:
    ## 
    ## # Simple named list: list(mean = mean, median = median)
    ## 
    ## # Auto named with `tibble::lst()`: tibble::lst(mean, median)
    ## 
    ## # Using lambdas list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
MOFA_weights <- as.data.frame(MOFA_weights)
MOFA_weights$mofa_weights <- MOFA_weights$mofa_weights * MOFA_weights$sign
MOFA_weights <- MOFA_weights[,-3]

#make metabolite names coherent between mofa and cosmos
#for now i just add the metab input that is already formatted
#I will miss some weights that were not input but will fix that at later time
MOFA_weights_metabaddon <- as.data.frame(metab_inputs)
names(MOFA_weights_metabaddon) <- "mofa_weights"
MOFA_weights_metabaddon$Nodes <- row.names(MOFA_weights_metabaddon)

MOFA_weights_metabaddon[, 2] <- sapply(MOFA_weights_metabaddon[, 2], function(x, HMDB_mapper_vec) {
        x <- gsub("Metab__", "", x)
        suffixe <- stringr::str_extract(x, "_[a-z]$")
        x <- gsub("_[a-z]$", "", x)
        if (x %in% names(HMDB_mapper_vec)) {
            x <- HMDB_mapper_vec[x]
            x <- paste("Metab__", x, sep = "")
        }
        if (!is.na(suffixe)) {
            x <- paste(x, suffixe, sep = "")
        }
        return(x)
    }, HMDB_mapper_vec = HMDB_mapper_vec)

MOFA_weights <- as.data.frame(rbind(MOFA_weights, MOFA_weights_metabaddon))


## Ligands
lig_weights <- ligrec_TF_moon_inputs$ligrec
names(lig_weights) <- gsub("___.+","",names(lig_weights))
lig_weights <- tapply(lig_weights, names(lig_weights), mean)
ligands <- names(lig_weights)
lig_weights <- as.numeric(lig_weights)
names(lig_weights) <- ligands
lig_weights <- data.frame(Nodes = names(lig_weights), feature_weights = lig_weights)

## Receptors
rec_weights <- ligrec_TF_moon_inputs$ligrec
names(rec_weights) <- gsub(".+___","",names(rec_weights))
rec_weights <- tapply(rec_weights, names(rec_weights), mean)
receptors <- names(rec_weights)
rec_weights <- as.numeric(rec_weights)
names(rec_weights) <- receptors
rec_weights <- data.frame(Nodes = names(rec_weights), feature_weights = rec_weights)

## TF
TF_weights <- ligrec_TF_moon_inputs$TF
names(TF_weights) <- gsub("_TF","",names(TF_weights))
TF_weights <- TF_weights[!(names(TF_weights) %in% c(rec_weights$Nodes,lig_weights$Nodes))]
TF_weights <- data.frame(Nodes = names(TF_weights), feature_weights = TF_weights)



## Combine
feature_weights <- as.data.frame(rbind(TF_weights, rec_weights, lig_weights))

## Add weights to data
combined_ATT <- merge(combined_ATT, MOFA_weights, all.x = T)
combined_ATT <- merge(combined_ATT, feature_weights, all.x = T)
```

To identify whether a node is a ligand, receptor or TFr, we first load
the consensus PKNs and only keep entries which represent an interaction
in our network.

``` r
# LIANA
load("support/ligrec_ressource.RData")
ligrec_df <- ligrec_ressource[,c("source_genesymbol","target_genesymbol")]
ligrec_df <- distinct(ligrec_df)
names(ligrec_df) <- c("Node1","Node2")

ligrec_df$Node1 <- gsub("-","_",ligrec_df$Node1)
ligrec_df$Node2 <- gsub("-","_",ligrec_df$Node2)
ligrec_df$Sign <- 1
ligrec_df$Weight <- 1

ligrec_df <- ligrec_df[ligrec_df$Node1 %in% combined_ATT$Nodes & ligrec_df$Node2 %in% combined_ATT$Nodes, ]

# Load Dorothea (TF) network
load("support/dorothea_df.RData")
dorothea_df <- dorothea_df[dorothea_df$source %in% combined_ATT$Nodes,]
```

Then, we add the information to the network (ligand = 2, receptor = 3,
transcription factor = 4, other = 1, no NodeType available = 0).

``` r
combined_ATT$NodeType <- ifelse(combined_ATT$Nodes %in% ligrec_df$Node1, 2, ifelse(combined_ATT$Nodes %in% ligrec_df$Node2, 3, ifelse(combined_ATT$Nodes %in% dorothea_df$source, 4, combined_ATT$NodeType))) #if node is a ligand set NodeType to 2, if node is a receptor set NodeType to 3, if node is a TF set NodeType to 4, else keep 0 or 1

write_csv(combined_ATT, file = "results/cosmos/ATT_mofamoon_combined.csv")
```

Finally, we merge both SIF together by summarizing multiple entries via
mean calculation and add found interactions from the filtered LIANA
consensus PKN.

``` r
combined_SIF <- as.data.frame(rbind(SIF_rec_to_TF, SIF_TF_to_lig))
combined_SIF <- as.data.frame(combined_SIF %>%
  group_by(Node1, Node2) %>%
  summarise(across(.fns= ~ mean(.x, na.rm = TRUE))))
```

    ## `summarise()` has grouped output by 'Node1'. You can override using the
    ## `.groups` argument.

``` r
combined_SIF <- as.data.frame(rbind(combined_SIF,ligrec_df))

res_mofa_combined <- list(combined_SIF, combined_ATT)
save(res_mofa_combined, file = "results/cosmos/res_mofamoon_combined.RData")
write_csv(combined_SIF, file = "results/cosmos/SIF_mofamoon_combined.csv")
```

To filter the final network based on node properties, the strategy shown
here can be used to find interesting ligand-receptor interactions. The
final result is saved and later inspected (“Network visualization”).

``` r
ligands <- combined_ATT[combined_ATT$NodeType == 2,"Nodes"] # 2 for ligand
receptors <- combined_ATT[combined_ATT$NodeType == 3,"Nodes"] # 3 for receptors

combined_SIF_reduced <- combined_SIF[!(combined_SIF$Node2 %in% receptors & !(combined_SIF$Node2 %in% combined_SIF$Node1)),] #no node that is a receptor but not a source
combined_SIF_reduced <- combined_SIF_reduced[!(combined_SIF_reduced$Node2 %in% ligands & !(combined_SIF_reduced$Node2 %in% combined_SIF_reduced$Node1)),] #no node that is a ligand but not a source
combined_SIF_reduced <- combined_SIF_reduced[!(combined_SIF_reduced$Node1 %in% receptors & !(combined_SIF_reduced$Node1 %in% combined_SIF_reduced$Node2)),] #no node that is a receptor but not a regulated target
combined_SIF_reduced <- combined_SIF_reduced[!(combined_SIF_reduced$Node1 %in% ligands & !(combined_SIF_reduced$Node1 %in% combined_SIF_reduced$Node2)),] #no node that is a ligand but not a regulated target 

combined_ATT_reduced <- combined_ATT[combined_ATT$Nodes %in% combined_SIF_reduced$Node1 | combined_ATT$Nodes %in% combined_SIF_reduced$Node2,]

res_mofa_combined_reduced <- list(combined_SIF_reduced, combined_ATT_reduced)
save(res_mofa_combined_reduced, file = "results/cosmos/res_mofamoon_combined_reduced.RData")
write_csv(combined_SIF_reduced, file = "results/cosmos/SIF_mofamoon_combined_reduced.csv")
write_csv(combined_ATT_reduced, file = "results/cosmos/ATT_mofamoon_combined_reduced.csv")
```

Also the specific sub-network of a gene (e.g. MYC) can be extracted by
taking into account n-step-distant nodes (here: 4 steps).

``` r
combined_SIF_reduced <- read_csv(file = "results/cosmos/SIF_mofamoon_combined_reduced.csv")
combined_ATT_reduced <- read_csv(file = "results/cosmos/ATT_mofamoon_combined_reduced.csv")
sig_prots <- combined_SIF_reduced[,c(1,4,2,3)]
names(sig_prots) <- c("source","interaction","target","sign")
background <- unique(c(sig_prots$source,sig_prots$target))

subnetwork_per_gene <- list()

for(node in background){
  
  SIF <- cosmosR:::keep_controllable_neighbours(sig_prots, n_steps = 4, input_nodes = node) #keeps the nodes in network that are no more then n_steps away from the starting nodes
  ATT <- combined_ATT_reduced[combined_ATT_reduced$Nodes %in% SIF$source | combined_ATT_reduced$Nodes %in% SIF$target,] #keep supplementing information for nodes
  
  write_csv(ATT, file = paste0("results/cosmos/subnetwork/ATT_subnetwork_",node,".csv"))
  write_csv(SIF, file = paste0("results/cosmos/subnetwork/SIF_subnetwork_",node,".csv"))
  
  subnetwork_per_gene[[node]] <- list("SIF" = SIF, "ATT" = ATT)
  save(subnetwork_per_gene, file = "results/cosmos/subnetwork/subnetwork_per_gene.RData")
}
```

### Pathway enrichment analysis

A potential downstream analysis of the biological network could be to
search for over-represented pathways. The tutorial
./Pathway_enrichment_analysis.Rmd explains possible steps.

### Network visualization via CytoScape

Further, we can analyze the biological network by visualization through
CytoScape. A tutorial for that is given under ./RCytoScape_tutorial.Rmd.

Here, a potential output of this analysis highlighting a sub pathway is
shown.

<div class="figure" style="text-align: center">

<embed src="results/cytoscape/sub_pathway_figure.pdf" title="COSMOS sub signaling pathway" width="800px" height="460px" type="application/pdf" />
<p class="caption">
COSMOS sub signaling pathway
</p>

</div>

### Further analyses

Further MOFA-COSMOS downstream analyses focusing on analyzing cell line
specific MOFA transcriptomics are given in
./Further_MOFA_COSMOS_analyses.Rmd.

### Session info

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
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] RCy3_2.16.0          tidyr_1.3.0          GSEABase_1.58.0     
    ##  [4] graph_1.74.0         annotate_1.74.0      XML_3.99-0.13       
    ##  [7] AnnotationDbi_1.58.0 IRanges_2.30.1       S4Vectors_0.34.0    
    ## [10] Biobase_2.56.0       BiocGenerics_0.42.0  gridExtra_2.3       
    ## [13] pheatmap_1.0.12      moon_0.1.0           decoupleR_2.5.2     
    ## [16] liana_0.1.5          reshape2_1.4.4       dplyr_1.1.1         
    ## [19] ggfortify_0.4.15     ggplot2_3.4.0        readr_2.1.4         
    ## [22] MOFA2_1.6.0          cosmosR_1.5.2       
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] utf8_1.2.3                  reticulate_1.28            
    ##   [3] tidyselect_1.2.0            RSQLite_2.2.20             
    ##   [5] htmlwidgets_1.6.1           grid_4.2.0                 
    ##   [7] BiocParallel_1.30.4         lpSolve_5.6.17             
    ##   [9] Rtsne_0.16                  devtools_2.4.5             
    ##  [11] base64url_1.4               munsell_0.5.0              
    ##  [13] ScaledMatrix_1.4.1          codetools_0.2-18           
    ##  [15] statmod_1.5.0               scran_1.24.1               
    ##  [17] pbdZMQ_0.3-9                future_1.30.0              
    ##  [19] miniUI_0.1.1.1              withr_2.5.0                
    ##  [21] colorspace_2.1-0            progressr_0.13.0           
    ##  [23] filelock_1.0.2              highr_0.10                 
    ##  [25] uuid_1.1-0                  logger_0.2.2               
    ##  [27] knitr_1.42                  rstudioapi_0.14            
    ##  [29] SingleCellExperiment_1.18.1 listenv_0.9.0              
    ##  [31] labeling_0.4.2              MatrixGenerics_1.8.1       
    ##  [33] repr_1.1.6                  GenomeInfoDbData_1.2.8     
    ##  [35] farver_2.1.1                bit64_4.0.5                
    ##  [37] rhdf5_2.40.0                basilisk_1.8.1             
    ##  [39] parallelly_1.34.0           vctrs_0.6.1                
    ##  [41] generics_0.1.3              xfun_0.38                  
    ##  [43] R6_2.5.1                    doParallel_1.0.17          
    ##  [45] GenomeInfoDb_1.32.4         clue_0.3-63                
    ##  [47] rsvd_1.0.5                  RJSONIO_1.3-1.8            
    ##  [49] locfit_1.5-9.7              bitops_1.0-7               
    ##  [51] rhdf5filters_1.8.0          cachem_1.0.7               
    ##  [53] DelayedArray_0.22.0         vroom_1.6.1                
    ##  [55] promises_1.2.0.1            scales_1.2.1               
    ##  [57] gtable_0.3.1                beachmat_2.12.0            
    ##  [59] globals_0.16.2              processx_3.8.0             
    ##  [61] rlang_1.1.0                 CARNIVAL_2.7.2             
    ##  [63] GlobalOptions_0.1.2         checkmate_2.1.0            
    ##  [65] BiocManager_1.30.19         yaml_2.3.7                 
    ##  [67] backports_1.4.1             httpuv_1.6.8               
    ##  [69] tools_4.2.0                 usethis_2.1.6              
    ##  [71] ellipsis_0.3.2              RColorBrewer_1.1-3         
    ##  [73] sessioninfo_1.2.2           Rcpp_1.0.10                
    ##  [75] plyr_1.8.8                  visNetwork_2.1.2           
    ##  [77] base64enc_0.1-3             sparseMatrixStats_1.8.0    
    ##  [79] progress_1.2.2              zlibbioc_1.42.0            
    ##  [81] purrr_1.0.1                 RCurl_1.98-1.10            
    ##  [83] ps_1.7.2                    basilisk.utils_1.8.0       
    ##  [85] prettyunits_1.1.1           GetoptLong_1.0.5           
    ##  [87] cowplot_1.1.1               urlchecker_1.0.1           
    ##  [89] SeuratObject_4.1.3          SummarizedExperiment_1.26.1
    ##  [91] ggrepel_0.9.2               cluster_2.1.4              
    ##  [93] fs_1.6.1                    magrittr_2.0.3             
    ##  [95] circlize_0.4.15             matrixStats_0.63.0         
    ##  [97] pkgload_1.3.2               hms_1.1.3                  
    ##  [99] mime_0.12                   evaluate_0.20              
    ## [101] xtable_1.8-4                readxl_1.4.2               
    ## [103] shape_1.4.6                 compiler_4.2.0             
    ## [105] tibble_3.2.1                crayon_1.5.2               
    ## [107] htmltools_0.5.5             later_1.3.0                
    ## [109] tzdb_0.3.0                  DBI_1.1.3                  
    ## [111] corrplot_0.92               ComplexHeatmap_2.12.1      
    ## [113] rappdirs_0.3.3              Matrix_1.5-3               
    ## [115] cli_3.6.1                   uchardet_1.1.1             
    ## [117] parallel_4.2.0              metapod_1.4.0              
    ## [119] igraph_1.4.2                GenomicRanges_1.48.0       
    ## [121] forcats_1.0.0               pkgconfig_2.0.3            
    ## [123] dir.expiry_1.4.0            OmnipathR_3.9.6            
    ## [125] sp_1.6-0                    IRdisplay_1.1              
    ## [127] scuttle_1.6.3               xml2_1.3.3                 
    ## [129] foreach_1.5.2               dqrng_0.3.0                
    ## [131] XVector_0.36.0              rvest_1.0.3                
    ## [133] stringr_1.5.0               callr_3.7.3                
    ## [135] digest_0.6.31               Biostrings_2.64.1          
    ## [137] rmarkdown_2.21              cellranger_1.1.0           
    ## [139] uwot_0.1.14                 edgeR_3.38.4               
    ## [141] DelayedMatrixStats_1.18.2   curl_5.0.0                 
    ## [143] shiny_1.7.4                 rjson_0.2.21               
    ## [145] lifecycle_1.0.3             jsonlite_1.8.4             
    ## [147] Rhdf5lib_1.18.2             BiocNeighbors_1.14.0       
    ## [149] limma_3.52.4                fansi_1.0.4                
    ## [151] pillar_1.9.0                lattice_0.20-45            
    ## [153] KEGGREST_1.36.3             fastmap_1.1.1              
    ## [155] httr_1.4.5                  pkgbuild_1.4.0             
    ## [157] glue_1.6.2                  remotes_2.4.2              
    ## [159] png_0.1-8                   iterators_1.0.14           
    ## [161] bit_4.0.5                   bluster_1.6.0              
    ## [163] stringi_1.7.12              profvis_0.3.7              
    ## [165] HDF5Array_1.24.2            blob_1.2.3                 
    ## [167] BiocSingular_1.12.0         memoise_2.0.1              
    ## [169] IRkernel_1.3.2              irlba_2.3.5.1              
    ## [171] future.apply_1.10.0
