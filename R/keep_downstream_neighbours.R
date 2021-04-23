#' keep controllable neighbors
#' 
#' keeps the nodes in network that are no more then n_steps away from the starting nodes.
#' @param network network in 3 column data.frame format
#'  (source, interaction, target)
#' @param n_steps largest distance t consider
#' @param input_nodes names of the input nodes in the network
#' @noRd
keep_controllable_neighbours <- function(network, n_steps, input_nodes)
{
    stopifnot(all(c("source","target","interaction") %in% colnames(network)))
    stopifnot(all(input_nodes %in% c(network$source,network$target)))
    
    print(paste("COSMOS: removing nodes that are not reachable from inputs within",n_steps,"steps"))
    meta_g <- igraph::graph_from_data_frame(network[,c("source",'target',"interaction")],directed = TRUE) 
    dn_nbours <- igraph::ego(graph = meta_g, order = n_steps,nodes = input_nodes, mode = "out")
    
    sub_nodes <- c(unique(names(unlist(dn_nbours))), input_nodes)
    
    to_keep = network$source %in% sub_nodes & network$target %in% sub_nodes
    
    print(paste("COSMOS:", sum(!to_keep), "from ", length(to_keep), 
                "interactions are removed from the PKN"))
    
    network <- network[to_keep,]
    
    
    return(network)
}

#' keep observable neighbors
#' 
#' keeps the nodes in network that are no more then n_steps upstreams from the 
#' measured nodes 
#' @param network network in 3 column data.frame format
#'  (source, interaction, target)
#' @param n_steps largest distance t consider
#' @param observed_nodes names of the measured nodes in the network
#' @noRd
keep_observable_neighbours <- function(network, n_steps, observed_nodes)
{
    stopifnot(all(c("source","target","interaction") %in% colnames(network)))
    stopifnot(all(observed_nodes %in% c(network$source,network$target)))
    
    print(paste("COSMOS: removing nodes that are not observable by measurements within",n_steps,"steps"))
    meta_g <- igraph::graph_from_data_frame(network[,c("source",'target',"interaction")],directed = TRUE) 
    up_nbours <- igraph::ego(graph = meta_g, order = n_steps, nodes = observed_nodes, mode = "in")
    
    sub_nodes <- c(unique(names(unlist(up_nbours))), observed_nodes)
    
    to_keep = network$source %in% sub_nodes & network$target %in% sub_nodes
    
    print(paste("COSMOS:", sum(!to_keep), "from ", length(to_keep), 
                "interactions are removed from the PKN"))
    
    network <- network[to_keep,]
    
    
    return(network)
}