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
    return(unlist(res))
    
}
