#' check_COSMOS_inputs
#' 
#' checks the format of the main data inputs. Checks only the data that was 
#' passed (not NULL). 
#' 
check_COSMOS_inputs <- function(meta_network = NULL,
                                tf_regulon = NULL,
                                signaling_data = NULL,
                                metabolic_data = NULL,
                                expression_data = NULL){
    # checking the data
    if(!is.null(signaling_data))  stopifnot(is.vector(signaling_data))
    if(!is.null(metabolic_data))  stopifnot(is.vector(metabolic_data))
    if(!is.null(expression_data))  stopifnot(is.vector(expression_data))
    
    # checking the networks
    if(!is.null(meta_network)){
        stopifnot(is.data.frame(meta_network))
        stopifnot(all(c("source","interaction","target" ) %in% names(meta_network)))
        stopifnot(ncol(meta_network)==3)
    }
    
    if(!is.null(tf_regulon)){
        stopifnot(is.data.frame(tf_regulon))
        stopifnot(all(c("tf","sign","target" ) %in% names(tf_regulon)))
        stopifnot(ncol(meta_network)==3)
    }
    return(TRUE)
}


