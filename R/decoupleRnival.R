#' moon
#'
#' Iteratively propagate downstream input activity through a signed directed network
#' using the weighted mean enrichment score from decoupleR package
#'
#' @param upstream_input A named vector with up_stream nodes and their corresponding activity.
#' @param downstream_input A named vector with down_stream nodes and their corresponding activity.
#' @param meta_network A network data frame containing signed directed prior knowledge of molecular interactions.
#' @param n_layers The number of layers that will be propagated upstream.
#' @param n_perm The number of permutations to use in decoupleR's algorithm.
#' @param downstream_cutoff If downstream measurments should be included above a given threshold
#' @param statistic the decoupleR stat to consider: "wmean", "norm_wmean", or "ulm"
#' @param return_level true or false, if true the layers that the protein belongs to will be returned alongside the scores
#'
#' @return A data frame containing the score of the nodes upstream of the 
#' downstream input based on the iterative propagation
#'
#' @export
#' 
#' @examples
#' # Example input data
#' upstream_input <- c("A" = 1, "B" = -1, "C" = 0.5)
#' downstream_input <- c("D" = 2, "E" = -1.5)
#' meta_network <- data.frame(
#'   source = c("A", "A", "B", "C", "C", "D", "E"),
#'   target = c("B", "C", "D", "E", "D", "B", "A"),
#'   sign = c(1, -1, -1, 1, -1, -1, 1)
#' )
#' 
#' # Run the function with the example input data
#' result <- moon(upstream_input, downstream_input, meta_network, n_layers = 2, n_perm = 100)
#' 
#' # View the results
#' print(result)
moon <- function(upstream_input = NULL, downstream_input, meta_network, n_layers, n_perm = 1000, downstream_cutoff = 0, statistic = "ulm", return_levels = F){
  
  
  regulons <- meta_network
  
  names(regulons)[names(regulons) == "sign" | names(regulons) == "interaction"] <- "mor"
  regulons <- regulons[!(regulons$source %in% names(downstream_input)),]
  
  switch(statistic,
         "norm_wmean" = {
           n_plus_one <- run_wmean(mat = as.matrix(data.frame(downstream_input)), network = regulons, times = n_perm, minsize = 1)
         },
         "wmean" = {
           n_plus_one <- run_wmean(mat = as.matrix(data.frame(downstream_input)), network = regulons, times = 2, minsize = 1)
         },
         "ulm" = {
           n_plus_one <- run_ulm(mat = as.matrix(data.frame(downstream_input)), network = regulons, minsize = 1)
         })
  
  n_plus_one <- n_plus_one[n_plus_one$statistic == statistic,c(2,4)]
  n_plus_one$level <- 1
  # regulons <- regulons[!(regulons$source %in% n_plus_one$source),]
  
  res_list <- list()
  res_list[[1]] <- as.data.frame(n_plus_one)
  i <- 2
  while(length(regulons[,1]) > 1 & sum(regulons$target %in% res_list[[i - 1]]$source) > 1 & i <= n_layers)
  {
    print(i)
    regulons <- regulons[!(regulons$source %in% res_list[[i - 1]]$source),]
    previous_n_plu_one <- res_list[[i - 1]][,-3,drop = F] #remove the layer indice
    row.names(previous_n_plu_one) <- previous_n_plu_one$source
    previous_n_plu_one <- previous_n_plu_one[,-1,drop = F]
    switch(statistic,
           "norm_wmean" = {
             n_plus_one <- run_wmean(mat = as.matrix(previous_n_plu_one), network = regulons, times = n_perm, minsize = 1)
           },
           "wmean" = {
             n_plus_one <- run_wmean(mat = as.matrix(previous_n_plu_one), network = regulons, times = 2, minsize = 1)
           },
           "ulm" = {
             n_plus_one <- run_ulm(mat = as.matrix(previous_n_plu_one), network = regulons, minsize = 1)
           })
    n_plus_one <- n_plus_one[n_plus_one$statistic == statistic,c(2,4)]
    regulons <- regulons[!(regulons$source %in% n_plus_one$source),]
    n_plus_one$level <- i
    res_list[[i]] <- as.data.frame(n_plus_one)
    i <- i +1
  }
  
  recursive_decoupleRnival_res <- as.data.frame(do.call(rbind,res_list))
  
  
  
  downstream_names <- as.data.frame(downstream_input)
  downstream_names$source <- row.names(downstream_names)
  names(downstream_names)[1] <- "score"
  downstream_names <- downstream_names[abs(downstream_names$score) > downstream_cutoff,]
  downstream_names$level <- 0
  
  recursive_decoupleRnival_res <- as.data.frame(rbind(recursive_decoupleRnival_res,downstream_names))
  
  if(!is.null(upstream_input))
  {
    upstream_input_df <- as.data.frame(upstream_input)
    upstream_input_df$source <- row.names(upstream_input_df)
    names(upstream_input_df)[1] <- "score"
    
    upstream_input_df <- merge(upstream_input_df, recursive_decoupleRnival_res, by = "source")
    upstream_input_df$filterout <- sign(upstream_input_df$score.x) != sign(upstream_input_df$score.y)
    
    recursive_decoupleRnival_res <- recursive_decoupleRnival_res[!(recursive_decoupleRnival_res$source %in% upstream_input_df[upstream_input_df$filterout,"source"]),]
  }
  
  return(return(recursive_decoupleRnival_res))
}


#' DecoupleRnival
#'
#' Iteratively propagate downstream input activity through a signed directed network
#' using the weighted mean enrichment score from decoupleR package
#'
#' @param upstream_input A named vector with up_stream nodes and their corresponding activity.
#' @param downstream_input A named vector with down_stream nodes and their corresponding activity.
#' @param meta_network A network data frame containing signed directed prior knowledge of molecular interactions.
#' @param n_layers The number of layers that will be propagated upstream.
#' @param n_perm The number of permutations to use in decoupleR's algorithm.
#' @param downstream_cutoff If downstream measurments should be included above a given threshold
#' @param statistic the decoupleR stat to consider: "wmean", "norm_wmean", or "ulm"
#'
#' @return A data frame containing the score of the nodes upstream of the 
#' downstream input based on the iterative propagation
#'
#' @export
#' 
#' @examples
#' # Example input data
#' upstream_input <- c("A" = 1, "B" = -1, "C" = 0.5)
#' downstream_input <- c("D" = 2, "E" = -1.5)
#' meta_network <- data.frame(
#'   source = c("A", "A", "B", "C", "C", "D", "E"),
#'   target = c("B", "C", "D", "E", "D", "B", "A"),
#'   sign = c(1, -1, -1, 1, -1, -1, 1)
#' )
#' 
#' # Run the function with the example input data
#' result <- decoupleRnival(upstream_input, downstream_input, meta_network, n_layers = 2, n_perm = 100)
#' 
#' # View the results
#' print(result)
decoupleRnival <- function(upstream_input = NULL, downstream_input, meta_network, n_layers, n_perm = 1000, downstream_cutoff = 0, statistic = "norm_wmean"){
  
  print("Warning, this function is deprecated and will no longer receive futur support. Please use the 'moon' function instead")
  regulons <- meta_network
  
  names(regulons)[names(regulons) == "sign" | names(regulons) == "interaction"] <- "mor"
  regulons <- regulons[!(regulons$source %in% names(downstream_input)),]
  
  switch(statistic,
         "norm_wmean" = {
           n_plus_one <- run_wmean(mat = as.matrix(data.frame(downstream_input)), network = regulons, times = n_perm, minsize = 1)
         },
         "wmean" = {
           n_plus_one <- run_wmean(mat = as.matrix(data.frame(downstream_input)), network = regulons, times = 2, minsize = 1)
         },
         "ulm" = {
           n_plus_one <- run_ulm(mat = as.matrix(data.frame(downstream_input)), network = regulons, minsize = 1)
         })
  
  n_plus_one <- n_plus_one[n_plus_one$statistic == statistic,c(2,4)]
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
    switch(statistic,
           "norm_wmean" = {
             n_plus_one <- run_wmean(mat = as.matrix(previous_n_plu_one), network = regulons, times = n_perm, minsize = 1)
           },
           "wmean" = {
             n_plus_one <- run_wmean(mat = as.matrix(previous_n_plu_one), network = regulons, times = 2, minsize = 1)
           },
           "ulm" = {
             n_plus_one <- run_ulm(mat = as.matrix(previous_n_plu_one), network = regulons, minsize = 1)
           })
    n_plus_one <- n_plus_one[n_plus_one$statistic == statistic,c(2,4)]
    regulons <- regulons[!(regulons$source %in% n_plus_one$source),]
    res_list[[i]] <- as.data.frame(n_plus_one)
    i <- i +1
  }
  
  recursive_decoupleRnival_res <- as.data.frame(do.call(rbind,res_list))
  
  
  
  downstream_names <- as.data.frame(downstream_input)
  downstream_names$source <- row.names(downstream_names)
  names(downstream_names)[1] <- "score"
  downstream_names <- downstream_names[abs(downstream_names$score) > downstream_cutoff,]
  
  recursive_decoupleRnival_res <- as.data.frame(rbind(recursive_decoupleRnival_res,downstream_names))
  
  if(!is.null(upstream_input))
  {
    upstream_input_df <- as.data.frame(upstream_input)
    upstream_input_df$source <- row.names(upstream_input_df)
    names(upstream_input_df)[1] <- "score"
    
    upstream_input_df <- merge(upstream_input_df, recursive_decoupleRnival_res, by = "source")
    upstream_input_df$filterout <- sign(upstream_input_df$score.x) != sign(upstream_input_df$score.y)
    
    recursive_decoupleRnival_res <- recursive_decoupleRnival_res[!(recursive_decoupleRnival_res$source %in% upstream_input_df[upstream_input_df$filterout,"source"]),]
  }
  
  return(return(recursive_decoupleRnival_res))
}

#' filter_incohrent_TF_target
#'
#' Filters incoherent target genes from a regulatory network based on a decoupling analysis
#' of upstream and downstream gene expression.
#'
#' @param decouplRnival_res A data frame resulting from the decoupleRnival function.
#' @param TF_reg_net A data frame containing prior knowledge of transcription factor (TF) regulatory interactions.
#' @param meta_network A network data frame containing signed directed prior knowledge of molecular interactions.
#' @param RNA_input A named vector containing differential gene expression data.
#'
#' @return A network data frame containing the genes that are not incoherently regulated by TFs.
#'
#' @import decoupleR
#' @import dplyr
#' @export
#' 
#' @examples
#' # Example input data
#' upstream_input <- c("A" = 1, "B" = -1, "C" = 0.5)
#' downstream_input <- c("D" = 2, "E" = -1.5)
#' meta_network <- data.frame(
#'   source = c("A", "A", "B", "C", "C", "D", "E"),
#'   target = c("B", "D", "D", "E", "D", "B", "A"),
#'   interaction = c(-1, 1, -1, 1, -1, -1, 1)
#' )
#' RNA_input <- c("A" = 1, "B" = -1, "C" = 5, "D" = -0.7, "E" = -0.3)
#' 
#' TF_reg_net <- data.frame(
#' source = c("B"),
#' target = c("D"),
#' mor = c(-1)
#' )
#'
#' # Run the decoupleRnival function to get the upstream influence scores
#' upstream_scores <- decoupleRnival(upstream_input, downstream_input, meta_network, n_layers = 2, n_perm = 100)
#'
#'filtered_network <- filter_incohrent_TF_target(upstream_scores, TF_reg_net, meta_network, RNA_input)
#'
#'print(filtered_network)
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

#' reduce_solution_network
#'
#' Reduces a solution network based on a decoupling analysis of upstream and downstream gene expression,
#' by filtering out edges that do not meet a consistency threshold, and limiting the network to a
#' certain number of steps from upstream input nodes.
#'
#' @param decoupleRnival_res A data frame resulting from the decoupleRnival function.
#' @param meta_network A network data frame containing signed directed prior knowledge of molecular interactions.
#' @param cutoff The consistency threshold for filtering edges from the solution network.
#' @param upstream_input A named vector with up_stream nodes and their corresponding activity.
#' @param RNA_input A named vector containing differential gene expression data.
#' @param n_steps The maximum number of steps from upstream input nodes to include in the solution network.
#'
#' @return A list containing the solution network (SIF) and an attribute table (ATT) with gene expression data.
#'
#' @import decoupleR
#'
#' @export
#' 
#' @examples
#' # Example input data
#' upstream_input <- c("A" = 1, "B" = -1, "C" = 0.5)
#' downstream_input <- c("D" = 2, "E" = -1.5)
#' meta_network <- data.frame(
#'   source = c("A", "A", "B", "C", "C", "D", "E"),
#'   target = c("B", "D", "D", "E", "D", "B", "A"),
#'   interaction = c(-1, 1, -1, 1, -1, -1, 1)
#' )
#' RNA_input <- c("A" = 1, "B" = -1, "C" = 5, "D" = 0.7, "E" = -0.3)
#'
#' # Run the decoupleRnival function to get the upstream influence scores
#' upstream_scores <- decoupleRnival(upstream_input, downstream_input, meta_network, n_layers = 2, n_perm = 100)
#'
#' # Reduce the solution network based on the upstream influence scores
#' reduced_network <- reduce_solution_network(upstream_scores, meta_network, 0.4, upstream_input, RNA_input, 3)
#'
#' # View the resulting solution network and attribute table
#' print(reduced_network$SIF)
#' print(reduced_network$ATT)
reduce_solution_network <- function(decoupleRnival_res, meta_network, cutoff, upstream_input, RNA_input = NULL, n_steps = 10)
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
  
  upstream_input_df <- as.data.frame(upstream_input)
  upstream_input_df$source <- row.names(upstream_input_df)
  names(upstream_input_df) <- c("score","nodes")
  
  upstream_input_df <- merge(upstream_input_df, recursive_decoupleRnival_res, by = "nodes")
  upstream_input_df$filterout <- sign(upstream_input_df$score.x) != sign(upstream_input_df$score.y)
  
  upstream_nodes <- upstream_input_df[!(upstream_input_df$filterout), "nodes"]
  upstream_nodes <- upstream_nodes[upstream_nodes %in% res_network$source | upstream_nodes %in% res_network$target]
  
  res_network <- cosmosR:::keep_controllable_neighbours(res_network, n_steps, upstream_nodes)
  
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

#' meta_network_cleanup
#'
#' This function cleans up a meta network data frame by removing self-interactions,
#' calculating the mean interaction values for duplicated source-target pairs, and
#' keeping only interactions with values of 1 or -1.
#'
#' @param meta_network A data frame with columns 'source', 'interaction', and 'target'.
#' @import dplyr
#' @return A cleaned up meta network data frame.
#' 
#' @export
#' 
#' @examples
#' # Create a meta network data frame
#' example_meta_network <- data.frame(
#' source = c("A", "B", "C", "D", "A", "B", "C", "D", "A"),
#' interaction = c(1, 1, 1, -1, 1, -1, 1, -1, 1),
#' target = c("B", "C", "D", "A", "C", "D", "A", "B", "B")
#' )
#'
#' # Clean up the example meta network
#' cleaned_meta_network <- meta_network_cleanup(example_meta_network)
#' print(cleaned_meta_network)
#'
meta_network_cleanup <- function(meta_network)
{
  if(sum(meta_network$source == meta_network$target) != 0)
  {
    meta_network <- meta_network[-which(meta_network$source == meta_network$target),]
  }  
  meta_network <- unique(meta_network)
  meta_network <- meta_network %>% group_by(source,target) %>% summarise_each(funs(mean(., na.rm = TRUE)))
  meta_network <- as.data.frame(meta_network)
  meta_network <- meta_network[meta_network$interaction %in% c(1,-1),]
  return(meta_network)
}

#' translate_res
#'
#' formats the network with readable names
#'
#' @param SIF  result SIF of decoupleRnival pipeline
#' @param ATT  result ATT of decoupleRnival pipeline
#' @param HMDB_mapper_vec a named vector with HMDB Ids as names and desired metabolite names as values.
#' @return list with network and attribute tables.
#' @importFrom stringr str_extract
#' @export
#' 
#' @examples
#' # Create a meta network data frame
#' example_SIF <- data.frame(
#' source = c("GPX1", "Gene863__GPX1"),
#' target = c("Gene863__GPX1", "Metab__HMDB0003337_c"),
#' sign = c(1, 1)
#' )
#' 
#' example_ATT <- data.frame(
#' Nodes = c("GPX1", "Gene863__GPX1","Metab__HMDB0003337_c"),
#' sign = c(1, 1, 1)
#' )
#' 
#' example_SIF
#' 
#' data("HMDB_mapper_vec")
#' 
#' translated_res <- translate_res(example_SIF,example_ATT,HMDB_mapper_vec)
#' 
#' translated_res$SIF
translate_res <- function(SIF,ATT,HMDB_mapper_vec = NULL)
{
  if (is.null(HMDB_mapper_vec)) {
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
  return(list("SIF" = SIF, "ATT" = ATT))
}