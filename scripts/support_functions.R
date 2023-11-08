format_LR_ressource <- function(ligrec_ressource)
{
  ligrec_geneset <- ligrec_ressource[,c("source_genesymbol","target_genesymbol")]
  ligrec_geneset$set <- paste(ligrec_geneset$source_genesymbol, ligrec_geneset$target_genesymbol, sep = "___")
  ligrec_geneset <- reshape2::melt(ligrec_geneset, id.vars = "set")[,c(3,1)]
  names(ligrec_geneset)[1] <- "gene"
  ligrec_geneset$mor <- 1
  ligrec_geneset$likelihood <- 1
  ligrec_geneset <- distinct(ligrec_geneset)
  
  return(ligrec_geneset)
}

wide_ulm_res <- function(ulm_result)
{
  ulm_result_df <- reshape2::dcast(ulm_result, formula = source~condition, value.var = "score")
  row.names(ulm_result_df) <- ulm_result_df$source
  ulm_result_df <- ulm_result_df[,-1]
  
  return(ulm_result_df)
}

make_heatmap_color_palette <- function(my_matrix)
{
  t <- as.vector(t(my_matrix))
  palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
  palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
  palette <- c(palette1, palette2)
}

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
  
  # Merge the moon_res data frame with the mapping table
  moon_res <- merge(moon_res, mapping_table, by = "source")
  
  # Return the merged data frame
  return(moon_res)
}

translate_column_HMDB <- function(my_column, HMDB_mapper_vec)
{
  return(sapply(my_column, function(x, HMDB_mapper_vec) {
    x <- gsub("Metab__", "", x)
    x <- gsub("^Gene", "Enzyme", x)
    suffixe <- stringr::str_extract(x, "_[a-z]$")
    x <- gsub("_[a-z]$", "", x)
    if (x %in% names(HMDB_mapper_vec)) {
      x <- HMDB_mapper_vec[x]
      x <- paste("Metab__", x, sep = "")
    }
    if (!is.na(suffixe)) {
      x <- paste(x, suffixe, sep = "")
    }
    return(x)
  }, HMDB_mapper_vec = HMDB_mapper_vec))
}