if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("saezlab/cosmosR")
devtools::install_github("saezlab/decoupleR")

library(cosmosR)
library(readr)
library(dplyr)
library(decoupleR)

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

#Remove genes that are not expressed from the meta_network
meta_network <- cosmosR:::filter_pkn_expressed_genes(names(RNA_input), meta_pkn = meta_network)

#Filter inputs and prune the meta_network to only keep nodes that can be found downstream of the inputs
#The number of step is quite flexible, 7 steps already covers most of the network

n_steps <- 10

# in this step we prune the network to keep only the relevant part between upstream and downstream nodes
sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)
meta_network <- cosmosR:::keep_controllable_neighbours(meta_network, n_steps, names(sig_input))
metab_input <- cosmosR:::filter_input_nodes_not_in_pkn(metab_input, meta_network)
meta_network <- cosmosR:::keep_observable_neighbours(meta_network, n_steps, names(metab_input))
sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)
#####
load("support/dorothea_reg.RData")

meta_network_TF_to_metab <- meta_network

before <- 1
after <- 0
i <- 1
while (before != after & i < 10) {
  before <- length(meta_network_TF_to_metab[,1])
  recursive_decoupleRnival_res <- decoupleRnival(upstream_input = sig_input, 
                                                              downstream_input = metab_input, 
                                                              meta_network = meta_network, 
                                                              n_layers = n_steps, 
                                                              n_perm = 100) # 1000 is better for definitive results
  
  meta_network_TF_to_metab <- filter_incohrent_TF_target(recursive_decoupleRnival_res, dorothea_reg, meta_network_TF_to_metab, RNA_input)
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
write_csv(recursive_decoupleRnival_res, file = paste("results/decoupleRnival/",paste(cell_line, "_ATT_decouplerino_full.csv",sep = ""), sep = ""))

plot(density(recursive_decoupleRnival_res$score))
abline(v = 1)
abline(v = -1)

solution_network <- reduce_solution_network(decoupleRnival_res = recursive_decoupleRnival_res, 
                                            meta_network = meta_network,
                                            cutoff = 1, 
                                            upstream_input = sig_input, 
                                            RNA_input = RNA_input, 
                                            n_steps = n_steps)

SIF <- solution_network$SIF
names(SIF)[3] <- "sign"
ATT <- solution_network$ATT


# You can export the table to csv to import them in cytoscape
write_csv(SIF, file = paste("results/decoupleRnival/",paste(cell_line, "_SIF.csv",sep = ""), sep = ""))
write_csv(ATT, file = paste("results/decoupleRnival/",paste(cell_line, "_ATT.csv",sep = ""), sep = ""))
