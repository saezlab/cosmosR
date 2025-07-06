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
#' @param return_levels true or false, if true the layers that the protein belongs to will be returned alongside the scores
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
#' result <- moon(upstream_input, downstream_input, meta_network, n_layers = 2, statistic = "wmean")
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
#' Extracts a subnetwork from a decoupleRnival result and a prior knowledge network,
#' enforcing score thresholds, neighbourhood distance, and structural constraints:
#' - Pure children (no outgoing edges) must be level 0 in decoupleRnival_res.
#' - Pure parents (no incoming edges) must be among the provided upstream_input nodes.
#' - All nodes must have |score| > cutoff.
#'
#' @param decoupleRnival_res A data.frame with columns `source`, `score`, and `level` from decoupleRnival().
#' @param meta_network A data.frame with columns `source`, `target`, `interaction` for signed directed edges.
#' @param cutoff Numeric. Absolute score threshold for node filtering.
#' @param upstream_input A named numeric vector of upstream seed nodes and their activity values.
#' @param RNA_input Optional named numeric vector of differential expression values; merged into ATT.
#' @param n_steps Integer. Maximum number of steps away from any upstream_input node.
#'
#' @return A list with:
#'   - SIF: data.frame of filtered edges (`source`, `target`, `interaction`, `consistency`).
#'   - ATT: data.frame of node attributes (`nodes`, `score`, `level`, `RNA_input`).
#'
#' @importFrom igraph graph_from_data_frame distances induced_subgraph 
#' @importFrom igraph delete_vertices degree V get.data.frame
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
#' # Run decoupleRnival to generate scores
#' dec_res <- moon(upstream_input, downstream_input, meta_network,
#'                          n_layers = 2, n_perm = 100)
#' # Extract solution network
#' sol_net <- reduce_solution_network(dec_res, meta_network,
#'                                   cutoff = 0.4, upstream_input,
#'                                   RNA_input, n_steps = 3)
#' # View outputs
#' print(sol_net$SIF)
#' print(sol_net$ATT)
reduce_solution_network <- function(decoupleRnival_res, meta_network, cutoff,
                                    upstream_input, RNA_input = NULL, n_steps = 10) {
  # Ensure required columns
  stopifnot(all(c("source", "score", "level") %in% colnames(decoupleRnival_res)))
  stopifnot(all(c("source", "target", "interaction") %in% colnames(meta_network)))
  
  # 1. Filter nodes by absolute score
  nodes_df <- decoupleRnival_res[abs(decoupleRnival_res$score) > cutoff, ]
  score_map <- setNames(nodes_df$score, nodes_df$source)
  valid_nodes <- nodes_df$source
  
  # 2. Build initial edge set and enforce consistency
  net <- meta_network[meta_network$source %in% valid_nodes &
                        meta_network$target %in% valid_nodes, ]
  net$consistency <- sign(score_map[net$source] * score_map[net$target]) == net$interaction
  net <- net[net$consistency, ]
  
  # 3. Create igraph
  g <- igraph::graph_from_data_frame(net[, c("source", "target")], directed = TRUE,
                                     vertices = data.frame(name = valid_nodes, stringsAsFactors = FALSE))
  
  # 4. Restrict to n_steps from upstream seeds
  seeds <- intersect(names(upstream_input), igraph::V(g)$name)
  if (length(seeds) == 0) {
    warning("No upstream_input nodes found in network; returning empty network.")
    return(list(SIF = data.frame(), ATT = data.frame()))
  }
  dists <- igraph::distances(g, v = seeds, to = igraph::V(g), mode = "out")
  keep_nodes <- names(which(apply(dists, 2, min) <= n_steps))
  g <- igraph::induced_subgraph(g, vids = keep_nodes)
  
  # 5. Prune pure children/parents by level/upstream rules
  repeat {
    deg_out <- igraph::degree(g, mode = "out")
    deg_in  <- igraph::degree(g, mode = "in")
    # identify pure children not at level 0
    pure_children <- names(deg_out)[deg_out == 0]
    bad_children <- pure_children[ nodes_df$level[match(pure_children, nodes_df$source)] != 0 ]
    # identify pure parents not in upstream_input
    pure_parents <- names(deg_in)[deg_in == 0]
    bad_parents <- setdiff(pure_parents, names(upstream_input))
    to_drop <- unique(c(bad_children, bad_parents))
    if (length(to_drop) == 0) break
    g <- igraph::delete_vertices(g, to_drop)
  }
  
  # 6. Prepare outputs
  # SIF
  sif <- igraph::get.data.frame(g, what = "edges")
  # rename for consistency
  colnames(sif)[colnames(sif) == "from"] <- "source"
  colnames(sif)[colnames(sif) == "to"]   <- "target"
  # merge interaction and consistency
  sif <- merge(sif, meta_network, by = c("source", "target"))
  sif$consistency <- sign(score_map[sif$source] * score_map[sif$target]) == sif$interaction
  sif <- sif[, c("source", "target", "interaction", "consistency")]
  
  # ATT
  final_nodes <- igraph::V(g)$name
  att <- nodes_df[nodes_df$source %in% final_nodes, c("source", "score", "level")]
  names(att)[names(att) == "source"] <- "nodes"
  if (!is.null(RNA_input)) {
    rna_df <- data.frame(nodes = names(RNA_input), RNA_input = as.numeric(RNA_input),
                         stringsAsFactors = FALSE)
    att <- merge(att, rna_df, by = "nodes", all.x = TRUE)
  } else {
    att$RNA_input <- NA_real_
  }
  
  return(list(SIF = sif, ATT = att))
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

#' Decompress Moon Result
#'
#' This function decompresses the results obtained from moon analysis by incorporating 
#' node signatures and handling duplicated parents. It merges these details with the 
#' provided meta network data and returns a comprehensive data frame.
#'
#' @param moon_res A data frame containing the results of a moon analysis.
#' @param meta_network_compressed_list A list containing compressed meta network details, 
#'        including node signatures and duplicated parents.
#' @param meta_network A data frame representing the original meta network.
#'
#' @return A data frame which merges the moon analysis results with the meta network data,
#'         including additional details about node signatures and handling of duplicated parents.
#'
#' @examples
#' # Example usage (requires appropriate data structures for moon_res, 
#' # meta_network_compressed_list, and meta_network)
#' # decompressed_result <- decompress_moon_result(moon_res, meta_network_compressed_list, meta_network)
#'
#' @export
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

#' Get Moon Scoring Network
#'
#' This function analyzes a given meta network based on moon scores and an upstream node. 
#' It filters and processes the network by controlling and observing neighbours 
#' according to specified parameters. The function returns a list containing a filtered 
#' network and updated moon scores.
#'
#' @param upstream_node The node from which the network analysis starts.
#' @param meta_network The complete network data.
#' @param moon_scores Scores associated with each node in the network.
#' @param keep_upstream_node_peers Logical; whether to keep peers of the upstream node. Default is FALSE.
#' 
#' @return A list with two elements: 
#'   - `SIF`: A data frame representing the filtered meta network.
#'   - `ATT`: A data frame representing the updated moon scores.
#' 
#' @examples
#' # Example usage (requires appropriate data structures for meta_network and moon_scores)
#' # result <- get_moon_scoring_network(upstream_node, meta_network, moon_scores)
#' 
#' @export
get_moon_scoring_network <- function(upstream_node,
                                     meta_network,
                                     moon_scores,
                                     keep_upstream_node_peers = F){
  
  # Determine the number of steps from the upstream node based on moon score level.
  n_steps <- moon_scores[moon_scores$source == upstream_node,"level"]
  
  # If level peers of the upstream node are not to be kept, filter out these nodes.
  if(!keep_upstream_node_peers)
  {
    moon_scores <- moon_scores[!(moon_scores$level == n_steps & moon_scores$source != upstream_node),]
  }
  
  # Filter the meta network to keep only controllable neighbours of the upstream node.
  meta_network_filtered <- cosmosR:::keep_controllable_neighbours(network = meta_network, n_steps = n_steps,input_nodes = upstream_node)
  
  # Identify downstream inputs from the moon scores.
  downstream_inputs <- moon_scores[which(moon_scores$level == 0 & moon_scores$source %in% meta_network_filtered$target),"source"]
  
  # Further filter the network to keep only observable neighbours.
  meta_network_filtered <- cosmosR:::keep_observable_neighbours(network = meta_network_filtered, n_steps = n_steps,observed_nodes = downstream_inputs)
  
  # Update moon scores to include only those present in the filtered network.
  moon_scores <- moon_scores[moon_scores$source %in% meta_network_filtered$source |
                               moon_scores$source %in% meta_network_filtered$target,]
  
  # Filter the network to include only those connections present in moon scores.
  meta_network_filtered <- meta_network_filtered[meta_network_filtered$source %in% moon_scores$source &
                                                   meta_network_filtered$target %in% moon_scores$source,]
  
  if(n_steps > 1 & !keep_upstream_node_peers)
  {
    remaining_level <- n_steps
    while (remaining_level >= 0) {
      top_nodes <- moon_scores[moon_scores$level == remaining_level,"source"]
      child_nodes <- meta_network_filtered[meta_network_filtered$source %in% top_nodes,"target"]
      
      moon_scores <- moon_scores[moon_scores$source %in% child_nodes | moon_scores$level != (remaining_level - 1),]
      meta_network_filtered <- meta_network_filtered[meta_network_filtered$source %in% moon_scores$source &
                                                       meta_network_filtered$target %in% moon_scores$source,]
      
      remaining_level <- remaining_level - 1
    }
  }
  
  # Return a list containing the filtered network (SIF) and the updated moon scores (ATT).
  return(list("SIF" = meta_network_filtered, "ATT" = moon_scores))
}

#' reduce_solution_network_double_thresh
#'
#' Extracts a subnetwork using a two-tier absolute-score threshold, then restricts
#' to only paths connecting upstream_input seeds to level 0 nodes, with recursive
#' consistency filtering to ensure edge signs match node activities.
#'
#' @param decoupleRnival_res A data.frame with columns `source`, numeric `score`, and integer `level`.
#' @param meta_network A data.frame with columns `source`, `target`, and `interaction`.
#' @param primary_thresh Numeric. Absolute score cutoff for primary node selection.
#' @param secondary_thresh Numeric. Absolute score cutoff for secondary node restriction.
#' @param upstream_input A named numeric vector of upstream seed nodes.
#' @param RNA_input Optional named numeric vector of differential expression values; merged into ATT.
#'
#' @return A list with:
#'   - SIF: data.frame of filtered edges (`source`,`target`,`interaction`).
#'   - ATT: data.frame of node attributes (`source`,`score`,`level`,`type`,`RNA_input`).
#'
#' @examples
#' dec_res <- data.frame(
#'   source = paste0("G", 1:6),
#'   score  = c(2.5, 1.2, 0.8, -2.2, 1.5, -0.5),
#'   level  = c(0, 0, 1, 0, 1, 1)
#' )
#' meta_net <- data.frame(
#'   source      = c("G1","G1","G2","G3","G4","G5"),
#'   target      = c("G2","G3","G4","G5","G2","G6"),
#'   interaction = c(1, -1, 1, 1, -1, 1)
#' )
#' upstream_input <- c(G1 = 1, G4 = -1)
#' RNA_input <- c(G1 = 0.5, G2 = -0.2, G4 = 1.1)
#' dbl_net <- reduce_solution_network_double_thresh(
#'   decoupleRnival_res = dec_res,
#'   meta_network       = meta_net,
#'   primary_thresh     = 2,
#'   secondary_thresh   = 1,
#'   upstream_input     = upstream_input,
#'   RNA_input          = RNA_input
#' )
#' print(dbl_net$SIF)
#' print(dbl_net$ATT)
#'
#' @importFrom igraph graph_from_data_frame neighborhood V vcount
#' @export
reduce_solution_network_double_thresh <- function(
    decoupleRnival_res,
    meta_network,
    primary_thresh,
    secondary_thresh,
    upstream_input,
    RNA_input = NULL
) {
  stopifnot(all(c("source","score","level") %in% names(decoupleRnival_res)))
  stopifnot(all(c("source","target","interaction") %in% names(meta_network)))
  
  # 1. Primary nodes
  prim_nodes <- decoupleRnival_res$source[abs(decoupleRnival_res$score) > primary_thresh]
  
  # 2. Sub-network touching primary nodes
  sub_net <- subset(meta_network,
                    source %in% prim_nodes | target %in% prim_nodes)
  # keep edges where both ends are primary or serve as intermediary
  cond_pp  <- sub_net$source %in% prim_nodes & sub_net$target %in% prim_nodes
  cond_int <- (sub_net$source %in% prim_nodes & sub_net$target %in% prim_nodes)
  sub_net <- sub_net[cond_pp | cond_int, ]
  
  # 3. Secondary filter
  sec_nodes <- decoupleRnival_res$source[abs(decoupleRnival_res$score) > secondary_thresh]
  sub_net <- subset(sub_net, source %in% sec_nodes & target %in% sec_nodes)
  
  # 4. Recursive consistency pruning
  score_map <- setNames(decoupleRnival_res$score, decoupleRnival_res$source)
  repeat {
    valid_edges <- with(sub_net,
                        sign(score_map[source] * score_map[target]) == interaction
    )
    if (all(valid_edges)) break
    sub_net <- sub_net[valid_edges, ]
    # rebuild graph and enforce seed->level0 connectivity
    g <- igraph::graph_from_data_frame(sub_net[,c("source","target")],
                                       directed = TRUE,
                                       vertices = unique(c(sub_net$source, sub_net$target)))
    seeds  <- intersect(names(upstream_input), igraph::V(g)$name)
    level0 <- intersect(decoupleRnival_res$source[decoupleRnival_res$level == 0], igraph::V(g)$name)
    from   <- unlist(igraph::neighborhood(g, order = igraph::vcount(g), nodes = seeds, mode = "out"))
    to     <- unlist(igraph::neighborhood(g, order = igraph::vcount(g), nodes = level0, mode = "in"))
    keep   <- intersect(igraph::V(g)$name[from], igraph::V(g)$name[to])
    sub_net <- subset(sub_net, source %in% keep & target %in% keep)
  }
  
  # 5. Final path restriction
  final_edges <- sub_net
  
  # 6. Build ATT
  final_nodes <- unique(c(final_edges$source, final_edges$target))
  att <- subset(decoupleRnival_res, source %in% final_nodes,
                select = c(source, score, level))
  att$type <- ifelse(att$source %in% names(upstream_input),
                     "upstream_input",
                     ifelse(att$level == 0, "level0", "other"))
  if (!is.null(RNA_input)) {
    rna_df <- data.frame(source = names(RNA_input),
                         RNA_input = as.numeric(RNA_input),
                         stringsAsFactors = FALSE)
    att <- merge(att, rna_df, by = "source", all.x = TRUE)
  } else {
    att$RNA_input <- NA_real_
  }
  
  return(list(SIF = final_edges, ATT = att))
}

