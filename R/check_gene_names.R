#' check_gene_names
#' 
#' checks if gene names looks like an EntrezID with a prefix of X. 
#' @param ... takes named vectors
#' @noRd
check_gene_names <- function(...){
    
    res = lapply(list(...),function(x){
        if(any(duplicated(names(x)))){
            stop("Duplicated gene names detected. Gene names in inputs should be unique.")
        }
        TRUE
    }
    )
    
    return(unlist(res))
    
}
