#' run COSMOS metabolism to signaling 
#' 
#' Runs COSMOS from metabolism to signaling.  This function uses CARNIVAL to find
#' a subset of the prior knowledge network based on optimization that (1) 
#' includes the most measured and input nodes and (2) which is in agreement with
#' the data.  Use \code{\link{preprocess_COSMOS_metabolism_to_signaling}} to 
#' prepare the the inputs, measurements and the prior knowledge network.
#'  
#' @param data \code{\link{cosmos_data}} object.  Use the 
#'   \code{\link{preprocess_COSMOS_metabolism_to_signaling}} function to create 
#'   an instance.   
#' @param CARNIVAL_options List that controls the options of CARNIVAL. See details 
#'   in \code{\link{default_CARNIVAL_options}}. 
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

run_COSMOS_metabolism_to_signaling <- function(data,
                                               CARNIVAL_options = default_CARNIVAL_options()){
    
    ## Checking COSMOS input format
    validate_cosmos_data_metabolism_to_signaling(data)
    
    check_CARNIVAL_options(CARNIVAL_options)
    
    disc_metabolic_data <- discretize_input(data$metabolic_data)
    
    CARNIVAL_results = runCARNIVAL_wrapper(network = data$meta_network,
                                           input_data = disc_metabolic_data,
                                           measured_data = data$signaling_data,
                                           options = CARNIVAL_options)

    network_results <- process_CARNIVAL_results(CARNIVAL_results)
    
}


validate_cosmos_data_metabolism_to_signaling <- function(data){
    
    validate_cosmos_data(data)
    
    if(!all(c("signaling_data_bin", "diff_expression_data_bin") %in% names(data)))
        stop("missing inputs detected. Input data should be obtained by running preprocess_cosmos_metabolism_to_signaling.")
    
}