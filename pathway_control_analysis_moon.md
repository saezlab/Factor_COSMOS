``` r
library(readr)
library(pheatmap)
library(cosmosR)
library(decoupleR)
library(GSEABase)
```

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: annotate

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: XML

    ## Loading required package: graph

    ## 
    ## Attaching package: 'graph'

    ## The following object is masked from 'package:XML':
    ## 
    ##     addNode

``` r
full_moon_res_combined <- as.data.frame(
  read_csv("results/cosmos/moon/full_moon_res_combined.csv"))
```

    ## Rows: 3122 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (1): source
    ## dbl (2): score, level
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
combined_meta_network_translated <- as.data.frame(
  read_csv("results/cosmos/moon/combined_meta_network_translated.csv"))
```

    ## Rows: 14027 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (2): source, target
    ## dbl (1): interaction
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
background_nodes <- full_moon_res_combined[abs(full_moon_res_combined$score) > 1,"source"]
```

``` r
top_nodes <- full_moon_res_combined[abs(full_moon_res_combined$score) > 2,"source"]

pathway_control_set <- list()
i <- 1
for(node_of_interest in top_nodes)
{
  downstream_nodes <- unique(unlist(cosmosR:::keep_controllable_neighbours(combined_meta_network_translated, n_steps = 2, input_nodes = node_of_interest)[,c(1,2)]))
  
  if(length(downstream_nodes) > 0)
  {
    downstream_nodes <- downstream_nodes[-which(downstream_nodes == node_of_interest)]
    downstream_nodes <- downstream_nodes[which(downstream_nodes %in% background_nodes)]
    if(length(downstream_nodes) > 0)
    {
      res_ORA <- as.data.frame(piano::runGSAhyper(genes = downstream_nodes, universe = background_nodes, gsc = piano::loadGSC(pathways_NABA_KEGG))$resTab)
      res_ORA$log2fold_ratio <- log2((res_ORA[,3]/(res_ORA[,3]+res_ORA[,4])) / (res_ORA[,5]/(res_ORA[,5]+res_ORA[,6])))
      res_ORA$node_of_interest <- node_of_interest
      res_ORA$pathway <- row.names(res_ORA)
  
      pathway_control_set[[i]] <- res_ORA 
      i <- i + 1
    }
  } 
}
```

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14025 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 2 genes of interest in 146 gene sets, using a background of 1101 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 10010 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 225 genes of interest in 146 gene sets, using a background of 878 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13803 from  14027 interactions are removed from the PKN"
    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13875 from  14027 interactions are removed from the PKN"
    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14006 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 12 genes of interest in 146 gene sets, using a background of 1091 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13834 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 37 genes of interest in 146 gene sets, using a background of 1066 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13345 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 30 genes of interest in 146 gene sets, using a background of 1073 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13947 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 11 genes of interest in 146 gene sets, using a background of 1092 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13975 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 12 genes of interest in 146 gene sets, using a background of 1091 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13968 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 13 genes of interest in 146 gene sets, using a background of 1090 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14022 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 3 genes of interest in 146 gene sets, using a background of 1100 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12828 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 85 genes of interest in 146 gene sets, using a background of 1018 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13949 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 15 genes of interest in 146 gene sets, using a background of 1088 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12867 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 67 genes of interest in 146 gene sets, using a background of 1036 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12978 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 59 genes of interest in 146 gene sets, using a background of 1044 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13960 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 14 genes of interest in 146 gene sets, using a background of 1089 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13970 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 11 genes of interest in 146 gene sets, using a background of 1092 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13980 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 12 genes of interest in 146 gene sets, using a background of 1091 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13601 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 20 genes of interest in 146 gene sets, using a background of 1083 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13048 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 45 genes of interest in 146 gene sets, using a background of 1058 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12521 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 83 genes of interest in 146 gene sets, using a background of 1020 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13406 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 71 genes of interest in 146 gene sets, using a background of 1032 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14012 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 5 genes of interest in 146 gene sets, using a background of 1098 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13961 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12872 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 77 genes of interest in 146 gene sets, using a background of 1026 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12583 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 52 genes of interest in 146 gene sets, using a background of 1051 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 11364 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 139 genes of interest in 146 gene sets, using a background of 964 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13885 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 23 genes of interest in 146 gene sets, using a background of 1080 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13197 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 39 genes of interest in 146 gene sets, using a background of 1064 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14025 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 2 genes of interest in 146 gene sets, using a background of 1101 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13328 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 62 genes of interest in 146 gene sets, using a background of 1041 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12208 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 72 genes of interest in 146 gene sets, using a background of 1031 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13728 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 18 genes of interest in 146 gene sets, using a background of 1085 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13977 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14025 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 2 genes of interest in 146 gene sets, using a background of 1101 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14025 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 2 genes of interest in 146 gene sets, using a background of 1101 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14025 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 2 genes of interest in 146 gene sets, using a background of 1101 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14013 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 3 genes of interest in 146 gene sets, using a background of 1100 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13998 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 5 genes of interest in 146 gene sets, using a background of 1098 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13996 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13982 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13963 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 12 genes of interest in 146 gene sets, using a background of 1091 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12966 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 75 genes of interest in 146 gene sets, using a background of 1028 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 11788 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 93 genes of interest in 146 gene sets, using a background of 1010 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 11833 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 128 genes of interest in 146 gene sets, using a background of 975 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14019 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13506 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 26 genes of interest in 146 gene sets, using a background of 1077 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13930 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 21 genes of interest in 146 gene sets, using a background of 1082 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 10960 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 140 genes of interest in 146 gene sets, using a background of 963 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14021 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 4 genes of interest in 146 gene sets, using a background of 1099 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 6 genes of interest in 146 gene sets, using a background of 1097 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12350 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 56 genes of interest in 146 gene sets, using a background of 1047 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13932 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 18 genes of interest in 146 gene sets, using a background of 1085 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14002 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 2 genes of interest in 146 gene sets, using a background of 1101 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14026 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 1 genes of interest in 146 gene sets, using a background of 1102 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14009 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14014 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 12 genes of interest in 146 gene sets, using a background of 1091 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14019 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14026 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 1 genes of interest in 146 gene sets, using a background of 1102 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14022 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 5 genes of interest in 146 gene sets, using a background of 1098 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14022 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 5 genes of interest in 146 gene sets, using a background of 1098 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14022 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 5 genes of interest in 146 gene sets, using a background of 1098 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14022 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 5 genes of interest in 146 gene sets, using a background of 1098 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14022 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 5 genes of interest in 146 gene sets, using a background of 1098 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14025 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 2 genes of interest in 146 gene sets, using a background of 1101 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14019 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14026 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 1 genes of interest in 146 gene sets, using a background of 1102 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14026 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 1 genes of interest in 146 gene sets, using a background of 1102 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 3 genes of interest in 146 gene sets, using a background of 1100 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14019 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14019 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14019 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14013 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 3 genes of interest in 146 gene sets, using a background of 1100 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 3 genes of interest in 146 gene sets, using a background of 1100 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 3 genes of interest in 146 gene sets, using a background of 1100 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14026 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 1 genes of interest in 146 gene sets, using a background of 1102 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14026 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 1 genes of interest in 146 gene sets, using a background of 1102 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14020 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 3 genes of interest in 146 gene sets, using a background of 1100 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13946 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 21 genes of interest in 146 gene sets, using a background of 1082 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13985 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14021 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 6 genes of interest in 146 gene sets, using a background of 1097 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13882 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 32 genes of interest in 146 gene sets, using a background of 1071 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13656 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 44 genes of interest in 146 gene sets, using a background of 1059 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12575 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 74 genes of interest in 146 gene sets, using a background of 1029 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 11888 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 101 genes of interest in 146 gene sets, using a background of 1002 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12960 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 64 genes of interest in 146 gene sets, using a background of 1039 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13945 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 15 genes of interest in 146 gene sets, using a background of 1088 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13184 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 40 genes of interest in 146 gene sets, using a background of 1063 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13955 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 14 genes of interest in 146 gene sets, using a background of 1089 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13590 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 36 genes of interest in 146 gene sets, using a background of 1067 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13616 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 32 genes of interest in 146 gene sets, using a background of 1071 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13318 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 67 genes of interest in 146 gene sets, using a background of 1036 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12517 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 57 genes of interest in 146 gene sets, using a background of 1046 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13953 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 16 genes of interest in 146 gene sets, using a background of 1087 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12675 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 69 genes of interest in 146 gene sets, using a background of 1034 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14021 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 4 genes of interest in 146 gene sets, using a background of 1099 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14011 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 5 genes of interest in 146 gene sets, using a background of 1098 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13165 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 81 genes of interest in 146 gene sets, using a background of 1022 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13306 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 72 genes of interest in 146 gene sets, using a background of 1031 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13389 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 67 genes of interest in 146 gene sets, using a background of 1036 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13389 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 67 genes of interest in 146 gene sets, using a background of 1036 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12360 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 81 genes of interest in 146 gene sets, using a background of 1022 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 11663 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 99 genes of interest in 146 gene sets, using a background of 1004 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 9897 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 193 genes of interest in 146 gene sets, using a background of 910 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13588 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 27 genes of interest in 146 gene sets, using a background of 1076 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14024 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 3 genes of interest in 146 gene sets, using a background of 1100 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13601 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 20 genes of interest in 146 gene sets, using a background of 1083 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14027 from  14027 interactions are removed from the PKN"
    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13941 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 11 genes of interest in 146 gene sets, using a background of 1092 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 11989 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 130 genes of interest in 146 gene sets, using a background of 973 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12796 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 87 genes of interest in 146 gene sets, using a background of 1016 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12907 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 84 genes of interest in 146 gene sets, using a background of 1019 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 6839 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 384 genes of interest in 146 gene sets, using a background of 719 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 8294 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 293 genes of interest in 146 gene sets, using a background of 810 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 7870 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 301 genes of interest in 146 gene sets, using a background of 802 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13966 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 14 genes of interest in 146 gene sets, using a background of 1089 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14015 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 5 genes of interest in 146 gene sets, using a background of 1098 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14015 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12902 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 84 genes of interest in 146 gene sets, using a background of 1019 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12763 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 89 genes of interest in 146 gene sets, using a background of 1014 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14011 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 10 genes of interest in 146 gene sets, using a background of 1093 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14011 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 12 genes of interest in 146 gene sets, using a background of 1091 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 8257 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 316 genes of interest in 146 gene sets, using a background of 787 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14009 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13989 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 14 genes of interest in 146 gene sets, using a background of 1089 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14016 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14015 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 8 genes of interest in 146 gene sets, using a background of 1095 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 6 genes of interest in 146 gene sets, using a background of 1097 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14015 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 7 genes of interest in 146 gene sets, using a background of 1096 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13856 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 28 genes of interest in 146 gene sets, using a background of 1075 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13880 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 24 genes of interest in 146 gene sets, using a background of 1079 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13340 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 36 genes of interest in 146 gene sets, using a background of 1067 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13075 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 77 genes of interest in 146 gene sets, using a background of 1026 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13790 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 15 genes of interest in 146 gene sets, using a background of 1088 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13276 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 34 genes of interest in 146 gene sets, using a background of 1069 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14021 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 4 genes of interest in 146 gene sets, using a background of 1099 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 11997 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 105 genes of interest in 146 gene sets, using a background of 998 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13356 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 68 genes of interest in 146 gene sets, using a background of 1035 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13911 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 17 genes of interest in 146 gene sets, using a background of 1086 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13531 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 26 genes of interest in 146 gene sets, using a background of 1077 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12232 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 123 genes of interest in 146 gene sets, using a background of 980 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12905 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 82 genes of interest in 146 gene sets, using a background of 1021 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13894 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 21 genes of interest in 146 gene sets, using a background of 1082 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13644 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 45 genes of interest in 146 gene sets, using a background of 1058 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13910 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 17 genes of interest in 146 gene sets, using a background of 1086 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13933 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 18 genes of interest in 146 gene sets, using a background of 1085 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12901 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 82 genes of interest in 146 gene sets, using a background of 1021 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 4 genes of interest in 146 gene sets, using a background of 1099 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 6 genes of interest in 146 gene sets, using a background of 1097 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 6 genes of interest in 146 gene sets, using a background of 1097 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14018 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 6 genes of interest in 146 gene sets, using a background of 1097 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13956 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 15 genes of interest in 146 gene sets, using a background of 1088 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 9854 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 204 genes of interest in 146 gene sets, using a background of 899 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13941 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 11 genes of interest in 146 gene sets, using a background of 1092 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 10365 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 199 genes of interest in 146 gene sets, using a background of 904 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 8098 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 298 genes of interest in 146 gene sets, using a background of 805 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 7982 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 327 genes of interest in 146 gene sets, using a background of 776 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 8970 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 255 genes of interest in 146 gene sets, using a background of 848 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 8853 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 262 genes of interest in 146 gene sets, using a background of 841 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 10764 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 184 genes of interest in 146 gene sets, using a background of 919 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13549 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 33 genes of interest in 146 gene sets, using a background of 1070 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 11008 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 122 genes of interest in 146 gene sets, using a background of 981 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13977 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 11 genes of interest in 146 gene sets, using a background of 1092 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12992 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 78 genes of interest in 146 gene sets, using a background of 1025 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13091 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 98 genes of interest in 146 gene sets, using a background of 1005 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13964 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 11 genes of interest in 146 gene sets, using a background of 1092 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13182 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 40 genes of interest in 146 gene sets, using a background of 1063 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 9592 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 217 genes of interest in 146 gene sets, using a background of 886 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12840 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 89 genes of interest in 146 gene sets, using a background of 1014 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12271 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 95 genes of interest in 146 gene sets, using a background of 1008 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12666 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 65 genes of interest in 146 gene sets, using a background of 1038 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12547 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 52 genes of interest in 146 gene sets, using a background of 1051 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12193 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 66 genes of interest in 146 gene sets, using a background of 1037 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14016 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 4 genes of interest in 146 gene sets, using a background of 1099 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12689 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 94 genes of interest in 146 gene sets, using a background of 1009 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12626 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 49 genes of interest in 146 gene sets, using a background of 1054 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13877 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 44 genes of interest in 146 gene sets, using a background of 1059 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 10264 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 146 genes of interest in 146 gene sets, using a background of 957 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12244 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 65 genes of interest in 146 gene sets, using a background of 1038 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12812 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 88 genes of interest in 146 gene sets, using a background of 1015 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13840 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 26 genes of interest in 146 gene sets, using a background of 1077 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 9050 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 282 genes of interest in 146 gene sets, using a background of 821 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 10743 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 137 genes of interest in 146 gene sets, using a background of 966 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 8416 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 296 genes of interest in 146 gene sets, using a background of 807 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13950 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 10 genes of interest in 146 gene sets, using a background of 1093 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13955 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 15 genes of interest in 146 gene sets, using a background of 1088 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 10406 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 165 genes of interest in 146 gene sets, using a background of 938 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14022 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 4 genes of interest in 146 gene sets, using a background of 1099 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13970 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 13 genes of interest in 146 gene sets, using a background of 1090 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14024 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 3 genes of interest in 146 gene sets, using a background of 1100 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13997 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12190 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 86 genes of interest in 146 gene sets, using a background of 1017 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14027 from  14027 interactions are removed from the PKN"
    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14021 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 4 genes of interest in 146 gene sets, using a background of 1099 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13881 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 24 genes of interest in 146 gene sets, using a background of 1079 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13777 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 33 genes of interest in 146 gene sets, using a background of 1070 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13197 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 39 genes of interest in 146 gene sets, using a background of 1064 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14023 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 3 genes of interest in 146 gene sets, using a background of 1100 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 11679 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 100 genes of interest in 146 gene sets, using a background of 1003 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12443 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 70 genes of interest in 146 gene sets, using a background of 1033 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13087 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 43 genes of interest in 146 gene sets, using a background of 1060 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 11905 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 149 genes of interest in 146 gene sets, using a background of 954 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14027 from  14027 interactions are removed from the PKN"
    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 13824 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 33 genes of interest in 146 gene sets, using a background of 1070 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 12324 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 84 genes of interest in 146 gene sets, using a background of 1019 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14025 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 1 genes of interest in 146 gene sets, using a background of 1102 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14015 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 9 genes of interest in 146 gene sets, using a background of 1094 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14025 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 2 genes of interest in 146 gene sets, using a background of 1101 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14021 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 3 genes of interest in 146 gene sets, using a background of 1100 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14025 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 2 genes of interest in 146 gene sets, using a background of 1101 non-interesting genes.

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 2 steps"
    ## [1] "COSMOS: 14025 from  14027 interactions are removed from the PKN"

    ## Warning in piano::runGSAhyper(genes = downstream_nodes, universe =
    ## background_nodes, : there are genes in gsc that are not in the universe, these
    ## will be removed before analysis

    ## Analyzing the overrepresentation of 2 genes of interest in 146 gene sets, using a background of 1101 non-interesting genes.

``` r
pathway_control_set <- do.call(rbind,pathway_control_set)
```

``` r
pathway_control_df <- reshape2::dcast(pathway_control_set, pathway~node_of_interest, value.var = "p-value")
row.names(pathway_control_df) <- pathway_control_df$pathway

pathway_control_df <- pathway_control_df[,-1]
pathway_control_df <- pathway_control_df[,apply(pathway_control_df, 2, function(x){min(x) < 0.1})]
```

``` r
threshold_pval <- 0.0000000000000001

pathway_control_df_top <- pathway_control_df[!grepl("CANCER",row.names(pathway_control_df)),]
pathway_control_df_top <- pathway_control_df_top[apply(pathway_control_df_top, 1, function(x){min(x) < threshold_pval}),apply(pathway_control_df_top, 2, function(x){min(x) < threshold_pval})]
pathway_control_df_top <- -log10(pathway_control_df_top)
# pathway_control_df_top[pathway_control_df_top < 3] <- NA
pathway_control_df_top[pathway_control_df_top >= 15] <- 15
pathway_control_df_top[pathway_control_df_top >= 8 & pathway_control_df_top < 15] <- 8
pathway_control_df_top[pathway_control_df_top >= 3 & pathway_control_df_top < 8] <- 3
pathway_control_df_top[pathway_control_df_top <3] <- 0

row.names(pathway_control_df_top) <- tolower(gsub("_"," ",gsub("KEGG","",row.names(pathway_control_df_top))))
names(pathway_control_df_top) <- gsub("Metab__","",gsub("_[a-z$]","",names(pathway_control_df_top)))
pheatmap::pheatmap(pathway_control_df_top, angle_col = 315, na_col = "grey", cluster_rows = T, cluster_cols = T, display_numbers = F, number_color = "black", color = colorRampPalette(c("white","red"))(100), treeheight_row = 0, treeheight_col = 0)
```

![](pathway_control_analysis_moon_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
pheatmap::pheatmap(pathway_control_df_top, angle_col = 315, na_col = "grey", cluster_rows = T, cluster_cols = T, display_numbers = F, number_color = "black", color = colorRampPalette(c("white","red"))(100), treeheight_row = 0, treeheight_col = 0, filename = "results/cosmos/moon/pathway_control_top.pdf", height = 3, width = 9)
```
