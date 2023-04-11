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

meta_network_compressed <- df

meta_network_compressed <- meta_network_cleanup(meta_network_compressed)

load("support/dorothea_reg.RData")

meta_network_TF_to_metab <- meta_network_compressed

before <- 1
after <- 0
i <- 1
while (before != after & i < 10) {
  before <- length(meta_network_TF_to_metab[,1])
  recursive_decoupleRnival_res <- decoupleRnival(upstream_input = sig_input, 
                                                 downstream_input = metab_input, 
                                                 meta_network = meta_network_TF_to_metab, 
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
# write_csv(recursive_decoupleRnival_res, file = paste("results/decoupleRnival/",paste(cell_line, "_ATT_decouplerino_full.csv",sep = ""), sep = ""))

duplicated_parents_df <- data.frame(duplicated_parents)
duplicated_parents_df$source_original <- row.names(duplicated_parents_df)
names(duplicated_parents_df)[1] <- "source"

addons <- data.frame(names(children_signature)[-which(names(children_signature) %in% duplicated_parents_df$source_original)]) 
names(addons)[1] <- "source"
addons$source_original <- addons$source

mapping_table <- as.data.frame(rbind(duplicated_parents_df,addons))

recursive_decoupleRnival_res <- merge(recursive_decoupleRnival_res, mapping_table, by = "source")
recursive_decoupleRnival_res <- recursive_decoupleRnival_res[,c(3,2)]
names(recursive_decoupleRnival_res)[1] <- "source"

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




X7860_ATT_decouplerino_full <- as.data.frame(read_csv("results/decoupleRnival/7860_ATT_decouplerino_full.csv"))

temp <- merge(recursive_decoupleRnival_res, X7860_ATT_decouplerino_full, by = "source")

plot(temp$score.x, temp$score.y)

temp$diff <- abs(temp$score.y - temp$score.x)
temp <- temp[order(temp$diff, decreasing = T),]
