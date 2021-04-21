
#' process_CARNIVAL_results
#' 
#' formats the raw CARNIVAL results to a more appealing list of networks. 
#' @param  CARNIVAL_results list of matrices received from 
#' CRANIVAL::run_CARNIVAL()
#' @return list with the following elements:
#' - `aggregated_network` the averaged networks found by optimization in a 
#' format of a Simple Interaction network, i.e. each row codes an edge
#' - `N_networks`: number of solutions found by the optimization
#' - `aggregated_network_node_attributes`: estimated node properties
#' - `individual_networks`: list of optimial networks found
#' - `individual_networks_node_attributes`: node activity in each network.
#' @importFrom rlang .data
#' @importFrom purrr map
process_CARNIVAL_results <- function(CARNIVAL_results){
    
    network_output = list()
    
    network_output$weightedSIF <- as_tibble(CARNIVAL_results$weightedSIF) %>%
        rename(Node1 = "Node1",
               Node2 = "Node2",
               Sign = "Sign",
               Weight = "Weight") %>%
        mutate(Sign = as.numeric(.data$Sign),
               Weight = as.numeric(.data$Weight)/100)
    
    network_output$N_networks = length(CARNIVAL_results$sifAll)
    
    network_output$nodesAttributes = 
        as_tibble(CARNIVAL_results$nodesAttributes) %>%
        mutate(ZeroAct = as.numeric(.data$ZeroAct)/100,
               UpAct = as.numeric(.data$UpAct)/100,
               DownAct = as.numeric(.data$DownAct)/100,
               AvgAct = as.numeric(.data$AvgAct)/100,
               ) 
    
    network_output$individual_networks = 
        purrr::map(CARNIVAL_results$sifAll,function(Net){
            as_tibble(Net) %>%
                rename(source = "Node1",
                       target = "Node2",
                       interaction = "Sign") %>%
                mutate(interaction = as.numeric(.data$interaction))
            
        } )
    
    
    network_output$individual_networks_node_attributes = 
        purrr::map(CARNIVAL_results$attributesAll,function(Net){
            tibble::as_tibble(Net) %>%
                rename(node = "Nodes",
                       activity = "Activity") %>%
                mutate(activity = as.numeric(.data$activity))
            
        } )
    
    return(network_output)
    
}