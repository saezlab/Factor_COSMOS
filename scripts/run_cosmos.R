library(cosmosR)
library(readr)

data("meta_network")
load("data/cosmos/cosmos_inputs.RData")

names(cosmos_inputs)

cell_line <- "COLO 205"

sig_input <- cosmos_inputs[[cell_line]]$TF_scores
metab_input <- cosmos_inputs[[cell_line]]$metabolomic
RNA_input <- cosmos_inputs[[cell_line]]$RNA

compartment_code <- "_c"

names(metab_input) <- paste(names(metab_input),compartment_code, sep = "")
names(metab_input) <- paste("Metab__",names(metab_input),sep = "")
##Filter sugnificant inputs
sig_input <- sig_input[abs(sig_input) > 2]
metab_input <- metab_input[abs(metab_input) > 1.5]

#In order to adapt options to users specification we can load them into a variable 
#that will then be passed to preprocess_COSMOS_signaling_to_metabolism CARNIVAL_options parameter
my_options <- default_CARNIVAL_options()

#Here the user should provide a path to its CPLEX executable (only cplex at the moment, other solvers will be documented soon !)
# my_options$solverPath <- "~/Documents/cplex" #or cbc solver executable
my_options$solverPath <- "~/Documents/cbc-osx/cbc" #or cbc solver executable
# my_options$solver <- "cplex" #or cbc
my_options$solver <- "cbc"
my_options$timelimit <- 1800
my_options$mipGAP <- 0.05
my_options$threads <- 6

metab_input <- cosmosR:::filter_input_nodes_not_in_pkn(metab_input, meta_network)
sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)

test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = meta_network,
                                                      signaling_data = sig_input,
                                                      metabolic_data = metab_input,
                                                      diff_expression_data = RNA_input,
                                                      maximum_network_depth = 4,
                                                      remove_unexpressed_nodes = T,
                                                      filter_tf_gene_interaction_by_optimization = T,
                                                      CARNIVAL_options = my_options
                                                      
)

my_options$timelimit <- 7200

test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                      CARNIVAL_options = my_options)

formatted_res <- format_COSMOS_res(test_result_for)

SIF <- formatted_res[[1]]
ATT <- formatted_res[[2]]

RNA_input_df <- data.frame(Nodes = names(RNA_input), t = RNA_input)
ATT <- merge(ATT, RNA_input_df, all.x = T)

write_csv(SIF, file = paste("results/",paste(cell_line, "_SIF.csv",sep = ""), sep = ""))
write_csv(ATT, file = paste("results/",paste(cell_line, "_ATT.csv",sep = ""), sep = ""))

