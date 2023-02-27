library(cosmosR)
library(readr)
library(dplyr)
library(decoupleR)

source("scripts/support_decoupleRnival.R")

data("meta_network")

meta_network <- meta_network[-which(meta_network$source == meta_network$target),]
meta_network <- unique(meta_network)
meta_network <- meta_network %>% group_by(source,target) %>% summarise_each(funs(mean(., na.rm = TRUE)))
meta_network <- as.data.frame(meta_network)
meta_network <- meta_network[meta_network$interaction %in% c(1,-1),]

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

sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)

meta_network <- cosmosR:::keep_controllable_neighbours(meta_network, 7, names(sig_input))
meta_network <- cosmosR:::filter_pkn_expressed_genes(names(RNA_input), meta_pkn = meta_network)

sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)
meta_network <- cosmosR:::keep_controllable_neighbours(meta_network, 7, names(sig_input))

#####
recursive_decoupleRnival_res <- decoupleRnival(upstream_input = sig_input, 
                                               downstream_input = metab_input, 
                                               meta_network = meta_network, 
                                               n_layers = 7)


# dorothea_reg <- get_dorothea(levels = c("A","B"))
load("support/dorothea_reg.RData")

meta_network <- filter_incohrent_TF_target(recursive_decoupleRnival_res, dorothea_reg, meta_network, RNA_input)
# temp <- filter_incohrent_TF_target(recursive_decoupleRnival_res, dorothea_reg, meta_network, RNA_input)

recursive_decoupleRnival_res <- decoupleRnival(upstream_input = sig_input, 
                                               downstream_input = metab_input, 
                                               meta_network = meta_network, 
                                               n_layers = 7)

meta_network <- filter_incohrent_TF_target(recursive_decoupleRnival_res, dorothea_reg, meta_network, RNA_input)

recursive_decoupleRnival_res <- decoupleRnival(upstream_input = sig_input, 
                                               downstream_input = metab_input, 
                                               meta_network = meta_network, 
                                               n_layers = 7)

meta_network <- filter_incohrent_TF_target(recursive_decoupleRnival_res, dorothea_reg, meta_network, RNA_input)

recursive_decoupleRnival_res <- decoupleRnival(upstream_input = sig_input, 
                                               downstream_input = metab_input, 
                                               meta_network = meta_network, 
                                               n_layers = 7)
#####
write_csv(recursive_decoupleRnival_res, file = "ATT_decouplerino_full.csv")

plot(density(recursive_decoupleRnival_res$score))
abline(v = 1)
abline(v = -1)

solution_network <- reduce_solution_network(decoupleRnival_res = recursive_decoupleRnival_res, 
                                            meta_network = meta_network,
                                            cutoff = 2, 
                                            sig_input = sig_input, 
                                            RNA_input = RNA_input)

SIF <- solution_network$SIF
ATT <- solution_network$ATT

write_csv(SIF, file = paste("results/decoupleRnival/",paste(cell_line, "_SIF.csv",sep = ""), sep = ""))
write_csv(ATT, file = paste("results/decoupleRnival/",paste(cell_line, "_ATT.csv",sep = ""), sep = ""))
