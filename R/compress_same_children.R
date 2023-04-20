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
  
  node_signatures <- sapply(parents, function(parent,df_signature){
    
    return(paste("parent_of_",paste0(unlist(df_signature[which(df_signature[,1] == parent),2]), collapse = "_____"), sep = ""))
  },df_signature = df_signature, USE.NAMES = T, simplify = F)
  
  dubs <- node_signatures[duplicated(node_signatures) & 
                            !(names(node_signatures) %in% names(metab_input) | 
                                names(node_signatures) %in% names(sig_input))]
  
  duplicated_parents <- unlist(node_signatures[which(node_signatures %in% dubs)])
  
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
  
  return(list("compressed_network" = df, "node_signatures" = node_signatures, "duplicated_signatures" = duplicated_parents))
}

#' Decompress Solution Network
#'
#' This function decompresses a solution network by mapping node signatures back to their original identifiers.
#' The input is a formatted solution network, a meta network, node signatures, and duplicated parents.
#' The function returns a list containing the decompressed solution network and attribute table.
#'
#' @param formatted_res A list containing the solution network and attribute table.
#' @param meta_network A data frame representing the meta network.
#' @param node_signatures A list of node signatures.
#' @param duplicated_parents A list of duplicated parents from the compression process.
#'
#' @return A list containing the following elements:
#'   \item{SIF}{A data frame representing the decompressed solution network.}
#'   \item{ATT}{A data frame containing the attributes of the decompressed solution network.}
#'
#' @examples
#' # Create a sample formatted_res
#' formatted_res <- list(
#'   SIF = data.frame(source = c("parent_of_D1", "D"),
#'                    target = c("D", "F"),
#'                    interaction = c(1, 1),
#'                    Weight = c(1, 1)),
#'   ATT = data.frame(Nodes = c("parent_of_D1", "D", "F"),
#'                    NodeType = c("","",""),
#'                    ZeroAct = c(0,0,0),
#'                    UpAct = c(1,1,1),
#'                    DownAct = c(0,0,0),
#'                    AvgAct = c(1,1,1),
#'                    measured = c(0,0,0),
#'                    Activity = c(1,1,1))
#' )
#'
#' # Create a sample meta_network
#' meta_network <- data.frame(source = c("A", "B", "D"),
#'                            target = c("D", "D", "F"),
#'                            interaction_type = c(1, 1, 1))
#'
#' # Define node_signatures and duplicated_parents
#' node_signatures <- list("A" = "parent_of_D1","B" = "parent_of_D1","D" = "parent_F1")
#' duplicated_parents <- c("A" = "parent_of_D1","B" = "parent_of_D1")
#'
#' # Decompress the solution network
#' result <- decompress_solution_network(formatted_res, meta_network, node_signatures, duplicated_parents)
#' decompressed_network <- result[[1]]
#' attribute_table <- result[[2]]
#'
#' @export
decompress_solution_network <- function(formatted_res, meta_network, node_signatures, duplicated_parents)
{
  SIF <- formatted_res[[1]]
  ATT <- formatted_res[[2]]
  
  duplicated_parents_df <- data.frame(duplicated_parents)
  print(duplicated_parents_df)
  duplicated_parents_df$source_original <- row.names(duplicated_parents_df)
  names(duplicated_parents_df)[1] <- "Nodes"
  
  
  addons <- data.frame(names(node_signatures)[-which(names(node_signatures) %in% duplicated_parents_df$source_original)]) 
  names(addons)[1] <- "Nodes"
  addons$source_original <- addons$Nodes
  
  mapping_table <- as.data.frame(rbind(duplicated_parents_df,addons))
  
  data("HMDB_mapper_vec")
  
  mapping_table[, 1] <- sapply(mapping_table[, 1], function(x, HMDB_mapper_vec) {
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
  
  mapping_table[, 2] <- sapply(mapping_table[, 2], function(x, HMDB_mapper_vec) {
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
  
  # View(ATT[!(ATT$Nodes %in% mapping_table$Nodes),])
  
  
  
  ATT <- merge(ATT, mapping_table, by = "Nodes")
  ATT <- ATT[,c(9,2:8)]
  names(ATT)[1] <- "Nodes"
  
  SIF <- meta_network
  
  SIF[, 1] <- sapply(SIF[, 1], function(x, HMDB_mapper_vec) {
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
  
  SIF[, 2] <- sapply(SIF[, 2], function(x, HMDB_mapper_vec) {
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
  
  SIF <- SIF[SIF$source %in% ATT$Nodes & 
               SIF$target %in% ATT$Nodes,]
  
  print(SIF)
  print(ATT)
  
  SIF$Weight <- apply(SIF, 1, function(x, ATT){
    source_act <- ATT[ATT$Nodes == x[1],"Activity"]
    target_act <- ATT[ATT$Nodes == x[2],"Activity"]
    coherence <- source_act * target_act 
    weight <- ifelse(sign(coherence) == sign(as.numeric(x[3])), min(abs(c(source_act, target_act))), 0)
    return(weight)
  }, ATT = ATT)
  
  formatted_res <- list(SIF,ATT)
  
  return(formatted_res)
}