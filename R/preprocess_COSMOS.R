#' preprocess_COSMOS
#' 
#' runs checks on the input data and simplifies the prior knowledge network.
#' Simplification includes the removal of (1) nodes that are not reachable from 
#' signaling nodes and  (2) interactions between transcription factors and target
#' genes if the target gene does not respond or the response is contradictory
#' with the change in the transcription factor activity. 
#' Optionally, further TF activities are estimated via network optimization and 
#' the interactions between TF and genes are filtered again. 
#' 
#' @param meta_network prior knowledge network. By default COSMOS use a PKN 
#' derived from Omnipath, STITCHdb and Recon3D. See details on the function 
#' \cite{\code{load_meta_pkn}}
#' @param signaling_data numerical vector, where names are signaling nodes 
#' in the PKN and values are from \{1, 0, -1\}. Continuous data will be 
#' discretized using the \link\code{sign} function.  
#' @param metabolic_data numerical vector, where names are metabolic nodes 
#' in the PKN and values are continuous values. These values are compared to 
#' with the the simulation [? range of values]
#' @param expression_data (optional) numerical vector, where names are gene
#' names and values are from \{1,-1\}
#' @param filter_tf_gene_interaction_by_optimization (default:TRUE), if TRUE then runs 
#' a network optimization that estimates TF activity not included in the inputs
#' and checks the consistency between the estimated activity and change in gene 
#' expression. Removes interactions where TF and gene expression are inconsistent 
#' @param solver_path argument passed to \link\code{CARNIVAL::runCARNIVAL}. 
#' used if filter_tf_gene_interaction_by_optimization is TRUE
#' @param solver argument passed to \link\code{CARNIVAL::runCARNIVAL}
#' used if filter_tf_gene_interaction_by_optimization is TRUE
#' @param time_limit argument passed to \link\code{CARNIVAL::runCARNIVAL}
#' used if filter_tf_gene_interaction_by_optimization is TRUE
#' @export
#' @import dorothea biomaRt igraph dplyr
#' @return named list with the following members: 
#'  - `meta_network`  filtered PKN
#'  - `tf_regulon`  TF - target regulatory network
#'  - `signaling_data_bin` binarised signaling data 
#'  - `metabolic_data`  metabolomics data
#'  - `expression_data_bin`  binarised gene expression data 
#'  - `optimized_network` initial optimized network if filter_tf_gene_interaction_by_optimization is TRUE. 
#'
preprocess_COSMOS <- function(meta_network = load_meta_pkn(),
                              tf_regulon = load_tf_regulon_dorothea(),
                              signaling_data,
                              metabolic_data,
                              expression_data,
                              filter_tf_gene_interaction_by_optimization = TRUE,
                              solver_path = NULL, 
                              solver = "cplex",
                              time_limit = 3600){
    
    ## Checking COSMOS input format
    check_COSMOS_inputs(meta_network,
                        tf_regulon,
                        signaling_data,
                        metabolic_data,
                        expression_data)
    
    # Check overlap among node names in the inputs
    check_network_data_coverage(meta_network,
                                tf_regulon,
                                signaling_data,
                                metabolic_data,
                                expression_data)
    
    # Preprocess PKN
    # - cut unreachable nodes from inputs
    meta_network <- keep_downstream_neighbours(
        network = meta_network,
        n_steps =  8,
        starting_nodes = names(signaling_data))
    
    
    # Filter TF -> target interaction from PKN if target expression not changing
    gene_expression_binarized <- binarize_with_sign(expression_data,
                                                    threshold = 1)
    
    filtered_meta_network <- filter_transcriptional_regulations(
        network = meta_network[,1:3], 
        gene_expression_binarized = gene_expression_binarized,
        signaling_data  = signaling_data,
        tf_regulon=tf_regulon[,c("tf","sign","target")])
    
    # run CARNIVAL from signaling to metabolism,
    # this may estimate the activity of other TF-s.
    if(filter_tf_gene_interaction_by_optimization){
        
        CARNIVAL_results = runCARNIVAL_wrapper(network = meta_network,
                                               input_data = sign(signaling_data),
                                               measured_data = metabolic_data,
                                               solver_path = solver_path,
                                               solver = solver,
                                               timelimit = timelimit,
                                               mipGAP = 0.2)
        
        # get the estimated activity of TFs from CARNIVAL results
        estimated_TF_activity <- get_TF_activity_from_CARNIVAL(CARNIVAL_results, tf_regulon$tf)
        
        filtered_meta_network <- filter_transcriptional_regulations(
            network = filtered_meta_network, 
            gene_expression_binarized = binarize_with_sign(expression_data,
                                                           threshold = 1),
            signaling_data  = estimated_TF_activity,
            tf_regulon=tf_regulon[,c("tf","sign","target")])
        
        
    }else{
        CARNIVAL_results = list()
    }
    
    results <- list(meta_network = filtered_meta_network,
                    tf_regulon = tf_regulon,
                    signaling_data_bin = sign(signaling_data),
                    metabolic_data = metabolic_data,
                    expression_data = gene_expression_binarized,
                    optimized_network = CARNIVAL_results)
}
