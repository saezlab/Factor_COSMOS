# We advise to instal from github to get the latest version of the tool.
# if (!requireNamespace("devtools", quietly = TRUE))
#   install.packages("devtools")
# 
# devtools::install_github("saezlab/cosmosR")

library(cosmosR)
library(reshape2)
library(readr)

data("meta_network")

meta_network <- meta_network_cleanup(meta_network)

load("data/cosmos/cosmos_inputs.RData")

names(cosmos_inputs)

cell_line <- "786-0"

sig_input <- cosmos_inputs[[cell_line]]$TF_scores
metab_input <- cosmos_inputs[[cell_line]]$metabolomic
RNA_input <- cosmos_inputs[[cell_line]]$RNA

#Choose which compartment to assign to the metabolic measurments
metab_input <- prepare_metab_inputs(metab_input, c("c","m"))

##Filter significant inputs
sig_input <- sig_input[abs(sig_input) > 2]
# metab_input <- metab_input[abs(metab_input) > 2]

#make sure the sig_input is a named vector

#Remove genes that are not expressed from the meta_network
meta_network <- cosmosR:::filter_pkn_expressed_genes(names(RNA_input), meta_pkn = meta_network)

#Filter inputs and prune the meta_network to only keep nodes that can be found downstream of the inputs
#The number of step is quite flexible, 7 steps already covers most of the network

n_steps <- 6

# in this step we prune the network to keep only the relevant part between upstream and downstream nodes
sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)
meta_network <- cosmosR:::keep_controllable_neighbours(meta_network, n_steps, names(sig_input))
metab_input <- cosmosR:::filter_input_nodes_not_in_pkn(metab_input, meta_network)
meta_network <- cosmosR:::keep_observable_neighbours(meta_network, n_steps, names(metab_input))
sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)


#compress the network
meta_network_compressed_list <- compress_same_children(meta_network, sig_input = sig_input, metab_input = metab_input)

meta_network_compressed <- meta_network_compressed_list$compressed_network

node_signatures <- meta_network_compressed_list$node_signatures

duplicated_parents <- meta_network_compressed_list$duplicated_signatures

meta_network_compressed <- meta_network_cleanup(meta_network_compressed)

load("support/dorothea_reg.RData")

meta_network_TF_to_metab <- meta_network_compressed

before <- 1
after <- 0
i <- 1
while (before != after & i < 10) {
  before <- length(meta_network_TF_to_metab[,1])
  moon_res <- moon(upstream_input = sig_input, 
                                                 downstream_input = metab_input, 
                                                 meta_network = meta_network_TF_to_metab, 
                                                 n_layers = n_steps, 
                                                 statistic = "ulm") 
  
  meta_network_TF_to_metab <- filter_incohrent_TF_target(moon_res, dorothea_reg, meta_network_TF_to_metab, RNA_input)
  after <- length(meta_network_TF_to_metab[,1])
  i <- i + 1
}

if(i < 10)
{
  print(paste("Converged after ",paste(i-1," iterations", sep = ""),sep = ""))
} else
{
  print(paste("Interupted after ",paste(i," iterations. Convergence uncertain.", sep = ""),sep = ""))
}

#####
write_csv(moon_res, file = paste("results/moon/",paste(cell_line, "_ATT_decouplerino_full.csv",sep = ""), sep = ""))

source("scripts/support_decompression.R")
moon_res <- decompress_moon_result(moon_res, meta_network_compressed_list, meta_network_TF_to_metab)

plot(density(moon_res$score))
abline(v = 1)
abline(v = -1)

solution_network <- reduce_solution_network(decoupleRnival_res = moon_res, 
                                            meta_network = meta_network,
                                            cutoff = 1, 
                                            upstream_input = sig_input, 
                                            RNA_input = RNA_input, 
                                            n_steps = n_steps)

SIF <- solution_network$SIF
names(SIF)[3] <- "sign"
ATT <- solution_network$ATT

data("HMDB_mapper_vec")

translated_res <- translate_res(SIF,ATT,HMDB_mapper_vec)

SIF <- translated_res[[1]]
ATT <- translated_res[[2]]

write_csv(SIF, file = paste("results/",paste(cell_line, "_dec_compressed_SIF.csv",sep = ""), sep = ""))
write_csv(ATT, file = paste("results/",paste(cell_line, "_dec_compressed_ATT.csv",sep = ""), sep = ""))


X7860_ATT_decouplerino_full <- as.data.frame(read_csv("results/moon/7860_ATT_decouplerino_full.csv"))

temp <- merge(moon_res, X7860_ATT_decouplerino_full, by = "source")

plot(temp$score.x, temp$score.y)

temp$diff <- abs(temp$score.y - temp$score.x)
temp <- temp[order(temp$diff, decreasing = T),]
