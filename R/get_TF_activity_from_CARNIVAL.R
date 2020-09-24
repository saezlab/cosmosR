#' get_TF_activity_from_CARNIVAL
#' 
#' screens the CARNIVAL results and obtains the activity of transcription factors
#' 
#' @param carnival_result list obtained from runCARNIVAL
#' @param TFs character vector of transcription factors using EntrezIDs 
#' @return named numerical vector with TF activity found in CARNIVAL results. 

get_TF_activity_from_CARNIVAL <- function(carnival_result, TFs)
{
    
    if(!validate_CARNIVAL_results(carnival_result)) warning("we failed to validate CARNIVAL results.")
    
    estimated_activity <- as_tibble(carnival_result$nodesAttributes) %>%
        dplyr::select(Node,AvgAct) %>%
        dplyr::filter(AvgAct!=0) %>%
        dplyr::filter(Node %in% TFs) 
    
    activity = as.numeric(estimated_activity$AvgAct)/100
    names(activity) = estimated_activity$Node
    
    return(activity)
}

#' validate_CARNIVAL_results
#' 
#' implement here basic checks to see if we got back from CARNIVAL that we
#' expected. Subject to change if CARNIVAL results changes
validate_CARNIVAL_results <- function(CR){
 
    list_names = c("weightedSIF","nodesAttributes","sifAll","attributesAll" ) %in% names(CR)
    if(!all(list_names)) return(FALSE)
    
    if(is.null(CR$weightedSIF)) return(FALSE)
    if(is.null(CR$nodesAttributes)) return(FALSE)
    if(is.null(CR$sifAll)) return(FALSE)
    if(is.null(CR$attributesAll)) return(FALSE)
    
    return(TRUE)
}
