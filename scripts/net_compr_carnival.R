library(cosmosR)
library(ocean)
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
metab_input <- metab_input[abs(metab_input) > 2]

#Remove genes that are not expressed from the meta_network
meta_network <- cosmosR:::filter_pkn_expressed_genes(names(RNA_input), meta_pkn = meta_network)

#Filter inputs and prune the meta_network to only keep nodes that can be found downstream of the inputs
#The number of step is quite flexible, 7 steps already covers most of the network

n_steps <- 4

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



#In order to adapt options to users specification we can load them into a variable 
#that will then be passed to preprocess_COSMOS_signaling_to_metabolism CARNIVAL_options parameter
my_options <- default_CARNIVAL_options(solver = "cplex")

#Here the user should provide a path to its CPLEX executable (only cplex at the moment, other solvers will be documented soon !)
# my_options$solverPath <- "~/Documents/cplex" #or cbc solver executable
# my_options$solverPath <- "./cplex"
my_options$solverPath <- "cplex_macos/cplex" #or cbc solver executable
# my_options$solver <- "cplex" #or cbc
my_options$solver <- "cplex"
my_options$timelimit <- 3600/10
my_options$mipGAP <- 0.05
my_options$threads <- 6

metab_input <- cosmosR:::filter_input_nodes_not_in_pkn(metab_input, meta_network_compressed)
sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network_compressed)

test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = meta_network_compressed,
                                                      signaling_data = sig_input,
                                                      metabolic_data = metab_input,
                                                      diff_expression_data = RNA_input,
                                                      maximum_network_depth = n_steps,
                                                      remove_unexpressed_nodes = F,
                                                      filter_tf_gene_interaction_by_optimization = F,
                                                      CARNIVAL_options = my_options)

my_options$timelimit <- 3600/10

test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                      CARNIVAL_options = my_options)

formatted_res <- format_COSMOS_res(test_result_for)

SIF <- formatted_res[[1]]
ATT <- formatted_res[[2]]

SIF <- SIF[which(SIF$Weight != 0),]

RNA_input_df <- data.frame(Nodes = names(RNA_input), t = RNA_input)
ATT <- merge(ATT, RNA_input_df, all.x = T)
ATT <- ATT[ATT$AvgAct != 0,]


write_csv(SIF, file = paste("results/",paste(cell_line, "_compressed_SIF.csv",sep = ""), sep = ""))
write_csv(ATT, file = paste("results/",paste(cell_line, "_compressed_ATT.csv",sep = ""), sep = ""))
