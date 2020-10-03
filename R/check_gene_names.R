#' check_gene_names
#' 
#' checks if gene names looks like an EntrezID with a prefix of X. 
#' @param ... takes named vectors
#' 
check_gene_names <- function(...){
    
    res = lapply(list(...),function(x){
        if(!all(grepl("^X",names(x)))){
            stop("gene names in inputs should be entrezID with a prefix character X. eg: X5062")
        }
        TRUE
    }
    )
    
    res = lapply(list(...),function(x){
        if(any(duplicated(x))){
            stop("Duplicated gene names detected. Gene names in inputs should be unique.")
        }
        TRUE
    }
    )
    
    return(unlist(res))
    
}
