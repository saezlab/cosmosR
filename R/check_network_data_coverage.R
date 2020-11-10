
#' checks coverage between data and prior knowledge
#' 
#' checks the following
#' - signaling data nodes are in PKN
#' - metabolic data nodes are in PKN
#' - expression data genes overlap with TF targets
#' Stops if a check fails
#'
#' @param meta_network contains the PKN
#' @param tf_regulon (optional) tf-target network with EntrezID
#' @param signaling_data numerical vector, where names are signaling nodes 
#' in the PKN and values are from \{1,-1\}
#' @param metabolic_data numerical vector, where names are metabolic nodes 
#' in the PKN and values are continuous values. These values are compared to 
#' with the the simulation [? range of values]
#' @param expression_data (optional) numerical vector, where names are gene names  
#' and values are from \{1,-1\}
#' @param expand_metabolic_data format metabolic data to match meta_network nodes?
#' It will add "XMetab__" before the pubchem IDs and all the possible compartments 
#' after it. (e.g, "1150" will become "Xmetab__1150___c____" and "Xmetab__1150___e____")
#' @param verbose (default: TRUE) reports coverage
check_network_data_coverage <- function(meta_network,
                                        tf_regulon = NULL,
                                        signaling_data,
                                        metabolic_data,
                                        expression_data = NULL,
                                        verbose = TRUE){
    
    # signaling should be in PKN
    signaling_nodes = names(signaling_data)
    missing_nodes <- signaling_nodes[!signaling_nodes %in% c(meta_network$source,
                                                             meta_network$target)]
    if(length(missing_nodes)>0){
        
        stop(paste("The following signaling nodes are not found in the PKN:",
                   limit_string_vec(missing_nodes)))    
        
        
    }else{
        if(verbose) print(paste("COSMOS: all", length(signaling_nodes),
                                "signaling nodes from data were found in the meta PKN"))
    }
    
    # metabolic nodes should be in PKN
    metabolic_nodes = names(metabolic_data)
    missing_nodes <- metabolic_nodes[!metabolic_nodes %in% c(meta_network$source,
                                                             meta_network$target)]
    if(length(missing_nodes)>0){
        
        stop(paste("The following metabolic nodes are not found in the PKN:",
                   limit_string_vec(missing_nodes)))    
        
        
        
        
    }else{
        if(verbose) print(paste("COSMOS: all",length(metabolic_nodes) ,
                                "metabolic nodes from data were found in the meta PKN"))
        
    }
    
    # the expression data should overlap with the TF targets
    if(!is.null(tf_regulon)){
        genes = names(expression_data)
        genes_as_tf_target <- sum(genes %in% tf_regulon$target)
        
        if(verbose) print(paste0("COSMOS: ", genes_as_tf_target,
                                 " of the ",length(genes),
                                 " genes in expression data were found as transcription factor target"))
        if(verbose) print(paste0("COSMOS: ", sum( unique(tf_regulon$target) %in% genes),
                                 " of the ",length(unique(tf_regulon$target)),
                                 " transcription factor targets were found in expression data"))
        if(genes_as_tf_target==0){
            stop("Expression data contains no gene that appear as transcription factor target.
             The expression_data must be a named vector using EntrezIDs.")
        }
    }
    
    
}



