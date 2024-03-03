decompress_moon_result <- function(moon_res, meta_network_compressed_list, meta_network) {
  # Extract node_signatures and duplicated_parents from the list
  node_signatures <- meta_network_compressed_list$node_signatures
  duplicated_parents <- meta_network_compressed_list$duplicated_signatures
  
  # Create a dataframe for duplicated parents
  duplicated_parents_df <- data.frame(duplicated_parents)
  duplicated_parents_df$source_original <- row.names(duplicated_parents_df)
  names(duplicated_parents_df)[1] <- "source"
  
  # Create a dataframe for addons
  addons <- data.frame(names(node_signatures)[-which(names(node_signatures) %in% duplicated_parents_df$source_original)]) 
  names(addons)[1] <- "source"
  addons$source_original <- addons$source
  
  # Get final leaves
  final_leaves <- meta_network[!(meta_network$target %in% meta_network$source),"target"]
  final_leaves <- as.data.frame(cbind(final_leaves,final_leaves))
  names(final_leaves) <- names(addons)
  
  # Combine addons and final leaves
  addons <- as.data.frame(rbind(addons,final_leaves))
  
  # Create mapping table by combining duplicated parents and addons
  mapping_table <- as.data.frame(rbind(duplicated_parents_df,addons))
  
  mapping_table <- unique(mapping_table)
  # Merge the moon_res data frame with the mapping table
  moon_res <- merge(moon_res, mapping_table, by = "source")
  
  # Return the merged data frame
  return(moon_res)
}
