#' run COSMOS signaling to metabolism
#' 
#' Runs COSMOS from signaling to metabolism.  This function uses CARNIVAL to find
#' a subset of the prior knowledge network based on optimisation that (1) 
#' includes the most measured and input nodes and (2) which is in agreement with
#' the data.  Use \code{\link{preprocess_COSMOS_signaling_to_metabolism}} to
#' prepare inputs, measurements and prior knowledge network.
#' 
#' @param data \code{\link{cosmos_data}} object.  Use the
#'   \code{\link{preprocess_COSMOS_signaling_to_metabolism}} function to create
#'   an instance. 
#' @param CARNIVAL_options List that controls the options of CARNIVAL. See 
#'   details in \code{\link{default_CARNIVAL_options}}. 
#' @export
#' @return List with the following elements:
#'   \describe{
#'     \item{\code{weightedSIF}}{The averaged networks found by 
#'     optimization in a format of a Simple Interaction network, i.e. each row 
#'     codes an edge}
#'     \item{\code{N_networks}}{Number of solutions found by the 
#'     optimization}
#'     \item{\code{nodesAttributes}}{Estimated node properties}
#'     \item{\code{individual_networks}}{List of optimial networks found}
#'     \item{\code{individual_networks_node_attributes}}{Node activity in each
#'     network}
#'   }
#' @seealso \code{\link{preprocess_COSMOS_metabolism_to_signaling}},  
#'   \code{\link[CARNIVAL]{runCARNIVAL}}, \code{\link{cosmos_data}}
#' @examples 
#' data(toy_network)
#' data(toy_signaling_input)
#' data(toy_metabolic_input)
#' data(toy_RNA)
#' test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = toy_network,
#'                              signaling_data = toy_signaling_input,
#'                              metabolic_data = toy_metabolic_input,
#'                              diff_expression_data = toy_RNA,
#'                              maximum_network_depth = 15,
#'                              remove_unexpressed_nodes = TRUE,
#'                              CARNIVAL_options = default_CARNIVAL_options("lpSolve"))
#' 
#' test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
#'                              CARNIVAL_options = default_CARNIVAL_options("lpSolve"))
run_COSMOS_signaling_to_metabolism <- function(data,
                                               CARNIVAL_options = default_CARNIVAL_options("lpSolve")){
    
    ## Checking COSMOS input format
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
    
    if(!all(c("signaling_data_bin") %in% names(data)))
        stop("missing inputs detected. Input data should be obtained by running preprocess_cosmos_signaling_to_metabolism.")
    
}

