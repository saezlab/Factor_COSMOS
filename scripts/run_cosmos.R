
# make sure that CARNIVAL and COSMOS are installed from github:
# if(!require("remotes")) install.packages("remotes")
remotes::install_github("saezlab/CARNIVAL")
remotes::install_github("saezlab/cosmosR", force = TRUE)

library(cosmosR)
library(readr)
library(dplyr)

data("meta_network")

#Seems to be an error in omnipath
meta_network <- meta_network[-which(meta_network$source == "PRKCA" & meta_network$target == "SRC"),]

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

#In order to adapt options to users specification we can load them into a variable 
#that will then be passed to preprocess_COSMOS_signaling_to_metabolism CARNIVAL_options parameter
my_options <- default_CARNIVAL_options(solver = "cbc")

#Here the user should provide a path to its CPLEX executable (only cplex at the moment, other solvers will be documented soon !)
# my_options$solverPath <- "~/Documents/cplex" #or cbc solver executable
# my_options$solverPath <- "./cplex"
my_options$solverPath <- "cbc/cbc-osx/cbc" #or cbc solver executable
# my_options$solver <- "cplex" #or cbc
my_options$solver <- "cbc"
my_options$timelimit <- 3600*0.5
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
                                                      CARNIVAL_options = my_options)

my_options$timelimit <- 3600*2

test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                      CARNIVAL_options = my_options)

formatted_res <- format_COSMOS_res(test_result_for)

SIF <- formatted_res[[1]]
ATT <- formatted_res[[2]]

SIF <- SIF[which(SIF$Weight != 0),]

RNA_input_df <- data.frame(Nodes = names(RNA_input), t = RNA_input)
ATT <- merge(ATT, RNA_input_df, all.x = T)
ATT <- ATT[ATT$AvgAct != 0,]


write_csv(SIF, file = paste("results/",paste(cell_line, "_SIF.csv",sep = ""), sep = ""))
write_csv(ATT, file = paste("results/",paste(cell_line, "_ATT.csv",sep = ""), sep = ""))

my_options$timelimit <- 3600*0.5


meta_network <- meta_network[-which(meta_network$source == meta_network$target),]


meta_network <- meta_network[-which(meta_network$source == meta_network$target),]

test_back <- preprocess_COSMOS_metabolism_to_signaling(meta_network = meta_network,
                                                      signaling_data = sig_input,
                                                      metabolic_data = metab_input,
                                                      diff_expression_data = RNA_input,
                                                      maximum_network_depth = 4,
                                                      remove_unexpressed_nodes = T,
                                                      filter_tf_gene_interaction_by_optimization = T,
                                                      CARNIVAL_options = my_options)

my_options$timelimit <- 3600*2

test_result_back <- run_COSMOS_metabolism_to_signaling(data = test_back,
                                                      CARNIVAL_options = my_options)

formatted_res_back <- format_COSMOS_res(test_result_back)

SIF_back <- formatted_res_back[[1]]
ATT_back <- formatted_res_back[[2]]

SIF_back <- SIF_back[which(SIF_back$Weight != 0),]

ATT_back <- merge(ATT_back, RNA_input_df, all.x = T)
ATT_back <- ATT_back[ATT_back$AvgAct != 0,]

write_csv(SIF_back, file = paste("results/",paste(cell_line, "_SIF_back.csv",sep = ""), sep = ""))
write_csv(ATT_back, file = paste("results/",paste(cell_line, "_ATT_back.csv",sep = ""), sep = ""))

SIF_full <- as.data.frame(rbind(SIF,SIF_back))
SIF_full <- unique(SIF_full)

ATT_full <- as.data.frame(rbind(ATT,ATT_back))
ATT_full <- unique(ATT_full)

S_nodes <- ATT_full[ATT_full$NodeType == "S","Nodes"]
T_nodes <- ATT_full[ATT_full$NodeType == "T","Nodes"]
C_nodes <- union(S_nodes,T_nodes)
# C_nodes <- C_nodes[which(C_nodes %in% c(SIF$Node1,SIF$Node2) & C_nodes %in% c(SIF_back$Node1,SIF_back$Node2))]

ATT_full <- ATT_full[,-6]

ATT_full <- ATT_full %>% group_by(Nodes) %>% summarise_each(funs(mean(., na.rm = TRUE)))
ATT_full <- as.data.frame(ATT_full)

ATT_full$NodeType <- ifelse(ATT_full$Nodes %in% C_nodes,"C",ifelse(ATT_full$Nodes %in% S_nodes,"S",ifelse(ATT_full$Nodes %in% T_nodes,"T","")))

write_csv(SIF_full, file = paste("results/",paste(cell_line, "_SIF_full.csv",sep = ""), sep = ""))
write_csv(ATT_full, file = paste("results/",paste(cell_line, "_ATT_full.csv",sep = ""), sep = ""))

