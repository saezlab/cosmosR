
filter_input_nodes_not_in_pkn <- function(data,pkn){
    
    new_data = data[names(data) %in% c(pkn$source,pkn$target)]
    
    if(length(data) != length(new_data)){
        removed_nodes <- names(data)[!names(data) %in% names(new_data)]
        
        print(paste("COSMOS:",length(removed_nodes), 
                    "input/measured nodes are not in PKN any more:",
                    limit_string_vec(removed_nodes)))
    }
    return(new_data)
}
