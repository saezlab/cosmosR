
# Constructor ------------------------------------------------------------------
# for internal use
#' Create New Cosmos Data
#' 
#' Constructor. Should not be exported. The object is not validated.
#' @param meta_network Prior knowledge network (PKN).  By default COSMOS use a 
#'   PKN derived from Omnipath, STITCHdb and Recon3D. See details on the data 
#'   \code{\link{meta_network}}.
#' @param tf_regulon Collection of transcription factor - target interactions.
#'   A default collection from dorothea can be obtained by the 
#'   \code{\link{load_tf_regulon_dorothea}} function.
#' @param signaling_data Numerical vector, where names are signaling nodes 
#'   in the PKN and values are from \{1, 0, -1\}.  Continuous data will be 
#'   discretized using the \code{\link{sign}} function.  
#' @param metabolic_data Numerical vector, where names are metabolic nodes 
#'   in the PKN and values are continuous values that represents log2 fold change 
#'   or t-values from a differential analysis.  These values are compared to 
#'   the simulation results (simulated nodes can take value -1, 0 or 1).
#' @param expression_data Numerical vector that represents the results of a 
#'   differential gene expression analysis. Names are gene names using EntrezID 
#'   starting with an X and values are log fold change or t-values. Genes with 
#'   NA values are considered none expressed and they will be removed from the
#'   TF-gene expression interactions. 
#' @noRd
new_cosmos_data <- function(meta_network = data.frame(),
                            tf_regulon = data.frame(),
                            signaling_data = double(),
                            metabolic_data = double(),
                            expression_data = double()){
    
    stopifnot(is.data.frame(meta_network))
    stopifnot(is.data.frame(tf_regulon))
    stopifnot(is.double(signaling_data))
    stopifnot(is.double(metabolic_data))
    stopifnot(is.double(expression_data))
    
    
    structure(list(
        meta_network = meta_network,
        tf_regulon = tf_regulon,
        signaling_data = signaling_data,
        metabolic_data = metabolic_data,
        expression_data = expression_data,
        history = c()
    ), class = "cosmos_data")
    
}

# validator  -------------------------------------------------------------------
#' Validate cosmos data
#' 
#' Checks the cosmos object.
#' 
#' @param x \code{\link{cosmos_data}} object.  Use the 
#'   \code{\link{preprocess_COSMOS_metabolism_to_signaling}} function to create 
#'   one.   
#' @noRd
validate_cosmos_data <- function(x){

    if(!all(c("meta_network", "tf_regulon", "signaling_data", "metabolic_data", "expression_data") %in% names(x))){
        stop(
            "All fields should present in the cosmos object: meta_network, tf_regulon, signaling_data, metabolic_data, expression_data",
            call. = FALSE
        )
    }
    
    meta_network = x$meta_network
    tf_regulon = x$tf_regulon
    signaling_data = x$signaling_data
    metabolic_data = x$metabolic_data
    expression_data = x$expression_data
    

    check_COSMOS_inputs(meta_network = meta_network,
                        tf_regulon = tf_regulon,
                        expression_data = expression_data,
                        signaling_data = signaling_data,
                        metabolic_data = metabolic_data)
    
    check_gene_names(signaling_data)
    
    # Check overlap among node names in the inputs
    check_network_data_coverage(meta_network = meta_network,
                                signaling_data = signaling_data,
                                metabolic_data = metabolic_data,
                                verbose = FALSE)
    
    
}





# helper -----------------------------------------------------------------------
#' Create Cosmos Data
#' 
#' An S3 class that combines the required data into a comprehensive list.  Use 
#' the \code{\link{preprocess_COSMOS_signaling_to_metabolism}} or 
#' \code{\link{preprocess_COSMOS_metabolism_to_signaling}} to create an instance.
#' 
#' @param meta_network Prior knowledge network (PKN).  By default COSMOS use a
#'   PKN derived from Omnipath, STITCHdb and Recon3D.  See details on the data 
#'   \code{\link{meta_network}}.
#' @param tf_regulon Collection of transcription factor - target interactions.  
#'   A default collection from dorothea can be obtained by the 
#'   \code{\link{load_tf_regulon_dorothea}} function.
#' @param signaling_data Numerical vector, where names are signaling nodes 
#'   in the PKN and values are from \{1, 0, -1\}.  Continuous data will be 
#'   discretized using the \code{\link{sign}} function.  
#' @param metabolic_data Numerical vector, where names are metabolic nodes 
#'   in the PKN and values are continuous values that represents log2 fold change 
#'   or t-values from a differential analysis.  These values are compared to 
#'   the simulation results (simulated nodes can take value -1, 0 or 1).
#' @param expression_data Numerical vector that represents the results of a 
#'   differential gene expression analysis. Names are gene names using EntrezID 
#'   starting with an X and values are log fold change or t-values. Genes with 
#'   NA values are considered none expressed and they will be removed from the
#'   TF-gene expression interactions.
#' @param verbose (default: TRUE) Reports details about the
#'   \code{\link{cosmos_data}} object. 
#' @return \code{cosmos data} class instance.  
#'   

cosmos_data <- function(meta_network,
                        tf_regulon = NULL,
                        signaling_data,
                        metabolic_data,
                        expression_data,
                        verbose = TRUE){
    
    check_gene_names(signaling_data)
    
    # signaling should be in PKN
    signaling_nodes = names(signaling_data)
    missing_nodes <- signaling_nodes[!signaling_nodes %in% c(meta_network$source,
                                                             meta_network$target)]
    if(length(missing_nodes)>0){
        message_missing_nodes <- paste("The following signaling nodes are not found in the PKN and will be removed from input:",
                                       limit_string_vec(missing_nodes))
        message("COSMOS: ", message_missing_nodes)
        signaling_data <- signaling_data[signaling_nodes != missing_nodes]
        
        
    }else{
        if(verbose) print(paste("COSMOS: all", length(signaling_nodes),
                                "signaling nodes from data were found in the meta PKN"))
    }
    
    # metabolic nodes should be in PKN
    metabolic_nodes = names(metabolic_data)
    missing_nodes <- metabolic_nodes[!metabolic_nodes %in% c(meta_network$source,
                                                             meta_network$target)]
    if(length(missing_nodes)>0){
        
        message_missing_nodes <- paste("The following metabolic nodes are not found in the PKN and will be removed from input:",
                                       limit_string_vec(missing_nodes))
        message("COSMOS: ", message_missing_nodes)
        # remove metabolic inputs not matching into the meta_network
        metabolic_data <- metabolic_data[metabolic_nodes != missing_nodes]
        
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
    
    new_cosmos_data(meta_network = meta_network,
                    tf_regulon = tf_regulon,
                    signaling_data = signaling_data,
                    metabolic_data = metabolic_data,
                    expression_data = expression_data)
    
    
    
}


### Generic functions ------------------------------------------------------------------
#' Print Cosmos Data Summary

#' Print a summary of cosmos data.

#' @param x \code{\link{cosmos_data}} object.  Use the 
#'   \code{\link{preprocess_COSMOS_metabolism_to_signaling}} function to 
#'   create one.
#' @param ... Further print arguments passed to or from other methods.
#' 
#' @seealso \code{\link[base]{print}}, \code{\link{cosmos_data}}
#'   
#' 
#' @return input (invisible)
#' @export
print.cosmos_data <- function(x, ...) {
    cat("cosmos_data contains the following data: \n")
    cat("Meta network:   ", nrow(x$meta_network), "interactions and", length(unique(c(x$meta_network$source,x$meta_network$target))), "unique nodes\n")
    cat("TF regulon:     ", nrow(x$tf_regulon), "interactions,", length(unique(c(x$tf_regulon$tf))),"TFs and", length(unique(c(x$tf_regulon$target))), "targets\n")
    cat("Signaling data: ", length(x$signaling_data), "measured nodes\n")
    cat("Metabolic data: ", length(x$metabolic_data), "measured nodes\n")
    cat("Expression data:", length(x$expression_data), "measured nodes\n")
    if(!is.null(x$history)) print(x$history)
    
    invisible(x)
}
