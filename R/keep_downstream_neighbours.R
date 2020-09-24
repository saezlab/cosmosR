#' keep downstream neighbors
#' 
#' keeps the nodes in network that are no more then n_steps away from the starting nodes.
#' @param network network in 3 column data.frame format
#'  (source, interaction, target)
#' @param n_steps largest distance t consider
#' @param starting_nodes names of nodes in the network
#' @import igraph
keep_downstream_neighbours <- function(network, n_steps, starting_nodes)
{
    stopifnot(all(c("source","target","interaction") %in% colnames(network)))
    stopifnot(all(starting_nodes %in% c(network$source,network$target)))
    
    meta_g <- igraph::graph_from_data_frame(network[,c("source",'target',"interaction")],directed = TRUE) 
    dn_nbours <- igraph::ego(graph = meta_g, order = n_steps,nodes = starting_nodes, mode = "out")
    
    sub_nodes <- c(unique(names(unlist(dn_nbours))), starting_nodes)
    
    to_keep = network$source %in% sub_nodes & network$target %in% sub_nodes
    
    print(paste("COSMOS:", sum(!to_keep), "from ", length(to_keep), 
                "interactions are removed from the PKN because nodes are not",
                "reachable from input nodes"))
    
    network <- network[to_keep,]
    
    return(network)
}
