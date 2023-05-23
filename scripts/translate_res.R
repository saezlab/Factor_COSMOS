translate_res <- function(SIF,ATT,metab_mapping)
{
  if (is.null(metab_mapping)) {
    data("HMDB_mapper_vec", package = "cosmosR", envir = environment())
  }
  colnames(ATT)[1] <- "Nodes"
  for (i in c(1, 2)) {
    SIF[, i] <- sapply(SIF[, i], function(x, HMDB_mapper_vec) {
      x <- gsub("Metab__", "", x)
      x <- gsub("^Gene", "Enzyme", x)
      suffixe <- stringr::str_extract(x, "_[a-z]$")
      x <- gsub("_[a-z]$", "", x)
      if (x %in% names(HMDB_mapper_vec)) {
        x <- HMDB_mapper_vec[x]
        x <- paste("Metab__", paste(x, suffixe, sep = ""), 
                   sep = "")
      }
      return(x)
    }, HMDB_mapper_vec = HMDB_mapper_vec)
  }
  ATT[, 1] <- sapply(ATT[, 1], function(x, HMDB_mapper_vec) {
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
  }, HMDB_mapper_vec = HMDB_mapper_vec)
  return(list(SIF, ATT))
}