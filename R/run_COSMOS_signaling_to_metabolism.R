#' run COSMOS signaling to metabolism
#' 
#' runs COSMOS from signaling to metabolism. This function finds a subset of the
#' prior knowledge network based on optimisation that (1) includes the most 
#' measured and input nodes and (2) which is in agreement with the data. 
#' Use \code{\link{preprocess_COSMOS_signaling_to_metabolism}} to prepare 
#' inputs, measurements and prior knowledge network.
#' @param data \code{\link{cosmos_data}} object. Use the \code{\link{preprocess_COSMOS_signaling_to_metabolism}}
#' function to create one. 
#' @param CARNIVAL_options list that controls the options of CARNIVAL. See details 
#'  in \code{\link{default_CARNIVAL_options}}. 
#' @export
#' @import dplyr
#' @return list with the following elements:
#' - `aggregated_network` the averaged networks found by optimization in a 
#' format of a Simple Interaction network, i.e. each row codes an edge
#' - `N_networks`: number of solutions found by the optimization
#' - `aggregated_network_node_attributes`: estimated node properties
#' - `individual_networks`: list of optimial networks found
#' - `individual_networks_node_attributes`: node activity in each network.

run_COSMOS_signaling_to_metabolism <- function(data,
                                               CARNIVAL_options = default_CARNIVAL_options()){
    
    validate_cosmos_data_signaling_to_metabolism(data)
    
    check_CARNIVAL_options(CARNIVAL_options)
    
    disc_signaling_data <- discretize_input(data$signaling_data)
    
    CARNIVAL_results = runCARNIVAL_wrapper(network = data$meta_network,
                                           input_data = disc_signaling_data,
                                           measured_data = data$metabolic_data,
                                           options = CARNIVAL_options)
    
    network_results <- process_CARNIVAL_results(CARNIVAL_results)
    
    
}


validate_cosmos_data_signaling_to_metabolism <- function(data){
    
    validate_cosmos_data(data)
    
    if(!all(c("signaling_data_bin", "diff_expression_data_bin") %in% names(data)))
        stop("missing inputs detected. Input data should be obtained by running preprocess_cosmos_signaling_to_metabolism.")
    
}

