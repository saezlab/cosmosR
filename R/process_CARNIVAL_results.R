
#' process_CARNIVAL_results
#' 
#' formats the raw CARNIVAL results to a more appealing list of networks. 
#' @param  CARNIVAL_results list of matrices received from 
#' CRANIVAL::run_CARNIVAL()
#' @return list with the following elements:
#'   \describe{
#'     \item{\code{aggregated_network}}{The averaged networks found by 
#'     optimization in a format of a Simple Interaction network, i.e. each row 
#'     codes an edge}
#'     \item{\code{N_networks}}{Number of solutions found by the 
#'     optimization}
#'     \item{\code{aggregated_network_node_attributes}}{Estimated node properties}
#'     \item{\code{individual_networks}}{List of optimial networks found}
#'     \item{\code{individual_networks_node_attributes}}{Node activity in each
#'     network}
#'   }
#' @importFrom rlang .data

process_CARNIVAL_results <- function(CARNIVAL_results){
    
    network_output = list()
    
    network_output$weightedSIF <- tibble::as_tibble(CARNIVAL_results$weightedSIF) %>%
        dplyr::rename(Node1 = "Node1",
               Node2 = "Node2",
               Sign = "Sign",
               Weight = "Weight") %>%
        dplyr::mutate(Sign = as.numeric(.data$Sign),
               Weight = as.numeric(.data$Weight)/100)
    
    network_output$N_networks = length(CARNIVAL_results$sifAll)
    
    network_output$nodesAttributes = 
        tibble::as_tibble(CARNIVAL_results$nodesAttributes) %>%
        dplyr::mutate(ZeroAct = as.numeric(.data$ZeroAct)/100,
               UpAct = as.numeric(.data$UpAct)/100,
               DownAct = as.numeric(.data$DownAct)/100,
               AvgAct = as.numeric(.data$AvgAct)/100,
               ) 
    
    network_output$individual_networks = 
        purrr::map(CARNIVAL_results$sifAll,function(Net){
            tibble::as_tibble(Net) %>%
                dplyr::rename(source = "Node1",
                       target = "Node2",
                       interaction = "Sign") %>%
                dplyr::mutate(interaction = as.numeric(.data$interaction))
            
        } )
    
    
    network_output$individual_networks_node_attributes = 
        purrr::map(CARNIVAL_results$attributesAll,function(Net){
            tibble::as_tibble(Net) %>%
                dplyr::rename(node = "Nodes",
                       activity = "Activity") %>%
                dplyr::mutate(activity = as.numeric(.data$activity))
            
        } )
    
    return(network_output)
    
}