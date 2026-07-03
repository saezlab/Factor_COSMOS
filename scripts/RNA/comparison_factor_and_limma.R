library(limma)
library(readr)
library(reshape2)

source("scripts/support_functions_limma.R")


load("results/mofa/z_matrix.RData")
load("support/all_metadata.RData")
load("results/mofa/weights_RNA.RData")
load("data/RNA/mofa_ready_data_RNA.RData")

expression_mat <- dcast(mofa_ready_data_RNA, feature~sample, value.var = "value")
row.names(expression_mat) <- expression_mat$feature
expression_mat <- expression_mat[,-1]

targets <- Z_matrix[,4,drop = F]
targets$sample <- row.names(targets)
targets <- targets[,c(2,1)]
names(targets) <- c("sample","condition")

targets$condition <- ifelse(targets$condition < -0.8, "condition_A","control")

unique(targets$condition)

comparisons <- list("condition_A_vs_control" = c(2,-1))

limmaRes <- runLimma(expression_mat, targets, comparisons = comparisons)

ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                    comp_names = names(comparisons),
                                    number = length(expression_mat[,1]),
                                    adjust.method = "fdr")

ttop <- ttop_list[[1]]

RNA_weights_Factor <- weights_RNA[,4,drop = F]
RNA_weights_Factor$ID <- gsub("_RNA","",row.names(RNA_weights_Factor))

to_compare <- merge(RNA_weights_Factor, ttop[,c(1,4)], by = "ID")

plot(to_compare$Factor4, to_compare$t)
cor.test(to_compare$Factor4, to_compare$t)
