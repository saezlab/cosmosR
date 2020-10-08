#' run_COSMOS_signaling2metabolism
#' 
#' runs COSMOS from signaling to metabolism. This function finds a subset of the
#' prior knowledge network based on optimisation that (1) includes the most 
#' measured and input nodes and (2) which is in agreement with the data. 
#' Use \code{\link{preprocess_COSMOS}} to prepare the prior knowledge network or 
#' load the one in the toolbox.  
#'  
#' @param meta_network prior knowledge network. 
#' @param signaling_data numerical vector, where names are signaling nodes 
#' in the PKN and values are from \{1, 0, -1\}. Continuous data will be 
#' discretized using the \code{\link{sign}} function.  
#' @param metabolic_data numerical vector, where names are metabolic nodes 
#' in the PKN and values are continuous values that represents log2 fold change 
#' or t values from a differential analysis. These values are compared to 
#' the simulation results (simulated nodes can take value -1, 0 or 1)
#' @param CARNIVAL_options list that controls the options of CARNIVAL. See details 
#'  in \code{\link{default_CARNIVAL_options()}}. 
#' @export
#' @import dplyr
#' @return named list with the following members: 
#'  - `meta_network`  filtered PKN
#'  - `tf_regulon`  TF - target regulatory network
#'  - `signaling_data_bin` binarised signaling data 
#'  - `metabolic_data`  metabolomics data
#'  - `optimized_network` initial optimized network if filter_tf_gene_interaction_by_optimization is TRUE. 


run_COSMOS_signaling_to_metabolism <- function(meta_network,
                                               signaling_data,
                                               metabolic_data,
                                               CARNIVAL_options = default_CARNIVAL_options(),
                                               test_run = FALSE){
    
    ## Checking COSMOS input format

    
    check_COSMOS_inputs(meta_network = meta_network,
                        signaling_data = signaling_data,
                        metabolic_data = metabolic_data)
    
    check_gene_names(signaling_data)
    
    # Check overlap among node names in the inputs
    check_network_data_coverage(meta_network = meta_network,
                                signaling_data = signaling_data,
                                metabolic_data = metabolic_data)
    
    if(!test_run){
        check_CARNIVAL_options(CARNIVAL_options)
        
        disc_signaling_data <- discretize_input(signaling_data)
        
        CARNIVAL_results = runCARNIVAL_wrapper(network = meta_network,
                                               input_data = disc_signaling_data,
                                               measured_data = metabolic_data,
                                               options = CARNIVAL_options)
        
    }else{
        CARNIVAL_results <- CARNIVAL_results_up_down
    }
    
    network_results <- process_CARNIVAL_results(CARNIVAL_results)
    
    
}


