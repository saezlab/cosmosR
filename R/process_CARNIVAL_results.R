
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
#' 
process_CARNIVAL_results <- function(CARNIVAL_results){
    
    
    network_output = list()
    
    network_output$aggregated_network <- as_tibble(CARNIVAL_results$weightedSIF) %>%
        rename(source = "Node1",
               target = "Node2",
               interaction = "Sign",
               weight = "Weight") %>%
        mutate(interaction = as.numeric(interaction),
               weight = as.numeric(weight)/100)
    
    network_output$N_networks = length(CARNIVAL_results$sifAll)
    
    network_output$aggregated_network_node_attributes = 
        as_tibble(CARNIVAL_results$nodesAttributes) %>%
        mutate(ZeroAct = as.numeric(ZeroAct)/100,
               UpAct = as.numeric(UpAct)/100,
               DownAct = as.numeric(DownAct)/100,
               AvgAct = as.numeric(AvgAct)/100,
               ) 
    
    network_output$individual_networks = 
        map(CARNIVAL_results$sifAll,function(Net){
            as_tibble(Net) %>%
                rename(source = "Node1",
                       target = "Node2",
                       interaction = "Sign") %>%
                mutate(interaction = as.numeric(interaction))
            
        } )
    
    
    network_output$individual_networks_node_attributes = 
        map(CARNIVAL_results$attributesAll,function(Net){
            as_tibble(Net) %>%
                rename(node = "Nodes",
                       activity = "Activity") %>%
                mutate(activity = as.numeric(activity))
            
        } )
    
    return(network_output)
    
}