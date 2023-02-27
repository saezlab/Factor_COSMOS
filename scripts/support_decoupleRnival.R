decoupleRnival <- function(upstream_input, downstream_input, meta_network, n_layers, n_perm = 1000){
  
  sig_input <- upstream_input
  metab_input <- downstream_input
  
  regulons <- meta_network
  
  names(regulons)[3] <- "mor"
  regulons <- regulons[!(regulons$source %in% names(metab_input)),]
  
  n_plus_one <- run_wmean(mat = data.frame(metab_input), network = regulons, times = n_perm, minsize = 1)
  n_plus_one <- n_plus_one[n_plus_one$statistic == "norm_wmean",c(2,4)]
  # regulons <- regulons[!(regulons$source %in% n_plus_one$source),]
  
  res_list <- list()
  res_list[[1]] <- as.data.frame(n_plus_one)
  i <- 2
  while(length(regulons[,1]) > 1 & sum(regulons$target %in% res_list[[i - 1]]$source) > 1 & i <= n_layers)
  {
    print(i)
    regulons <- regulons[!(regulons$source %in% res_list[[i - 1]]$source),]
    previous_n_plu_one <- res_list[[i - 1]]
    row.names(previous_n_plu_one) <- previous_n_plu_one$source
    previous_n_plu_one <- previous_n_plu_one[,-1,drop = F]
    n_plus_one <- run_wmean(mat = as.matrix(previous_n_plu_one), network = regulons, times = n_perm, minsize = 1)
    n_plus_one <- n_plus_one[n_plus_one$statistic == "norm_wmean",c(2,4)]
    regulons <- regulons[!(regulons$source %in% n_plus_one$source),]
    res_list[[i]] <- as.data.frame(n_plus_one)
    i <- i +1
  }
  
  recursive_decoupleRnival_res <- as.data.frame(do.call(rbind,res_list))
  
  
  
  metabs <- as.data.frame(metab_input)
  metabs$source <- row.names(metabs)
  names(metabs)[1] <- "score"
  metabs <- metabs[abs(metabs$score) > 2,]
  
  recursive_decoupleRnival_res <- as.data.frame(rbind(recursive_decoupleRnival_res,metabs))
  
  sig_input_df <- as.data.frame(sig_input)
  sig_input_df$source <- row.names(sig_input_df)
  names(sig_input_df)[1] <- "score"
  
  sig_input_df <- merge(sig_input_df, recursive_decoupleRnival_res, by = "source")
  sig_input_df$filterout <- sign(sig_input_df$score.x) != sign(sig_input_df$score.y)
  
  recursive_decoupleRnival_res <- recursive_decoupleRnival_res[!(recursive_decoupleRnival_res$source %in% sig_input_df[sig_input_df$filterout,"source"]),]
  
  return(return(recursive_decoupleRnival_res))
}

filter_incohrent_TF_target <- function(decouplRnival_res, TF_reg_net, meta_network, RNA_input){
  recursive_decoupleRnival_res <- decouplRnival_res
  dorothea_reg <- TF_reg_net
  
  RNA_df <- data.frame(RNA_input)
  RNA_df$node <- row.names(RNA_df)
  
  reg_meta <- recursive_decoupleRnival_res[recursive_decoupleRnival_res$source %in% dorothea_reg$source,]
  reg_meta <- merge(reg_meta,dorothea_reg, by.x = "source", by.y = "source")
  names(reg_meta)[2] <- "TF_score"
  # print(names(RNA_df))
  reg_meta <- merge(reg_meta, RNA_df, by.x = "target", by.y = "node")
  reg_meta$incoherent <- sign(reg_meta$TF_score * reg_meta$RNA_input * reg_meta$mor) < 0
  # View(reg_meta)
  reg_meta$edgeID <- paste(reg_meta$source, reg_meta$target, sep = "_")
  incoherent_edges <- reg_meta[reg_meta$incoherent, "edgeID"]
  # View(data.frame(incoherent_edges))
  # View(reg_meta)
  meta_network$edgeID <- paste(meta_network$source, meta_network$target, sep = "_")
  # View(meta_network)
  meta_network <- meta_network[!(meta_network$edgeID %in% incoherent_edges),]
  
  return(meta_network[,names(meta_network) != "edgeID"]) 
}

reduce_solution_network <- function(decoupleRnival_res, meta_network, cutoff, sig_input, RNA_input = NULL)
{
  recursive_decoupleRnival_res <- decoupleRnival_res
  
  recursive_decoupleRnival_res <- recursive_decoupleRnival_res[abs(recursive_decoupleRnival_res$score) > cutoff,]
  consistency_vec <- recursive_decoupleRnival_res$score
  names(consistency_vec) <- recursive_decoupleRnival_res$source
  
  res_network <- meta_network[meta_network$source %in% recursive_decoupleRnival_res$source & meta_network$target %in% recursive_decoupleRnival_res$source,]
  res_network$consistency <- sign(consistency_vec[res_network$source] * consistency_vec[res_network$target]) == res_network$interaction
  res_network <- res_network[res_network$consistency,]
  
  names(recursive_decoupleRnival_res)[1] <- "nodes"
  
  names(res_network)[3] <- "interaction"
  
  sig_input_df <- as.data.frame(sig_input)
  sig_input_df$source <- row.names(sig_input_df)
  names(sig_input_df) <- c("score","nodes")
  
  sig_input_df <- merge(sig_input_df, recursive_decoupleRnival_res, by = "nodes")
  sig_input_df$filterout <- sign(sig_input_df$score.x) != sign(sig_input_df$score.y)
  
  upstream_nodes <- sig_input_df[!(sig_input_df$filterout), "nodes"]
  upstream_nodes <- upstream_nodes[upstream_nodes %in% res_network$source | upstream_nodes %in% res_network$target]
  
  res_network <- cosmosR:::keep_controllable_neighbours(res_network, 7, upstream_nodes)
  
  SIF <- res_network
  ATT <- recursive_decoupleRnival_res[recursive_decoupleRnival_res$nodes %in% SIF$source | recursive_decoupleRnival_res$nodes %in% SIF$target,]
  
  if(!is.null(RNA_input))
  {
    RNA_df <- data.frame(RNA_input)
    RNA_df$nodes <- row.names(RNA_df)
    
    ATT <- merge(ATT, RNA_df, all.x = T)
  } else
  {
    ATT$RNA_input <- NA
  }
  return(list("SIF" = SIF, "ATT" = ATT))
}