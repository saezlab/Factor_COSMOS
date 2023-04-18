#' Compress Network by Merging Nodes with Identical Children
#'
#' This function compresses a network by merging nodes that have the same children.
#' The input network is represented as a data frame with three columns: source, target, and sign of interaction.
#' The function returns a list containing the compressed network, node signatures, and duplicated signatures.
#'
#' @param df A data frame representing the network with three columns: source, target, and sign of interaction.
#' @param sig_input A list of input node signatures to be considered for the merging process.
#' @param metab_input A list of input metabolic signatures to be considered for the merging process.
#'
#' @return A list containing the following elements:
#'   \item{compressed_network}{A data frame representing the compressed network.}
#'   \item{node_signatures}{A list of signatures of nodes in the network after the merging process.}
#'   \item{duplicated_signatures}{A list of duplicated signatures in the network after the merging process.}
#'
#' @examples
#' # Create a sample network
#' df <- data.frame(source = c("A", "A", "B", "B"),
#'                  target = c("C", "D", "C", "D"),
#'                  sign_of_interaction = c(1, 1, 1, 1))
#'
#' # Define input node and metabolic signatures
#' sig_input <- list()
#' metab_input <- list()
#'
#' # Compress the network
#' result <- compress_same_children(df, sig_input, metab_input)
#' compressed_network <- result$compressed_network
#'
#' @export
compress_same_children <- function(df, sig_input, metab_input)
{
  nodes <- unique(c(df$source,df$target))
  
  parents <- nodes[which(nodes %in% df$source)]
  
  df_signature <- df
  df_signature[,2] <- paste(df_signature[,2],df_signature[,3],sep = "")
  
  children_signature <- sapply(parents, function(parent,df_signature){
    
    return(paste("parent_of_",paste0(unlist(df_signature[which(df_signature[,1] == parent),2]), collapse = "_____"), sep = ""))
  },df_signature = df_signature, USE.NAMES = T, simplify = F)
  
  dubs <- children_signature[duplicated(children_signature) & 
                               !(names(children_signature) %in% names(metab_input) | 
                                   names(children_signature) %in% names(sig_input))]
  
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
  
  return(list("compressed_network" = df, "node_signatures" = children_signature, "duplicated_signatures" = dubs))
}

