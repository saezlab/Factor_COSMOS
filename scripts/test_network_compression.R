library(cosmosR)
library(ocean)
library(reshape2)

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

df <- meta_network

nodes <- unique(c(df$source,df$target))

parents <- nodes[which(nodes %in% meta_network$source)]

children_signature <- sapply(parents, function(parent,df){
  return(paste("parent_of_",paste0(unlist(df[which(df[,1] == parent),2]), collapse = "_____"), sep = ""))
},df = df, USE.NAMES = T, simplify = F)

dubs <- children_signature[duplicated(children_signature)]

duplicated_parents <- unlist(children_signature[which(children_signature %in% dubs)])

df[,1] <- sapply(df[,1], function(node,duplicated_parents){
  if(node %in% names(duplicated_parents))
  {
    node <- duplicated_parents[node]
  }
  return(node)
},duplicated_parents = duplicated_parents, simplify = T)

df[,2] <- sapply(df[,2], function(node,duplicated_parents){
  if(node %in% names(duplicated_parents))
  {
    node <- duplicated_parents[node]
  }
  return(node)
},duplicated_parents = duplicated_parents, simplify = T)

df <- unique(df)

meta_network <- df

meta_network <- meta_network_cleanup(meta_network)

load("support/dorothea_reg.RData")

recursive_decoupleRnival_res <- decoupleRnival(upstream_input = sig_input, 
                                               downstream_input = metab_input, 
                                               meta_network = meta_network, 
                                               n_layers = n_steps, 
                                               n_perm = 100) # 1000 is better for definitive results
