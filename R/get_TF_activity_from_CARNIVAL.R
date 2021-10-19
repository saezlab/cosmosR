#' get_TF_activity_from_CARNIVAL
#' 
#' screens the CARNIVAL results and obtains the activity of transcription factors
#' 
#' @param carnival_result list obtained from runCARNIVAL
#' @param TFs character vector of transcription factors using EntrezIDs
#' @return named numerical vector with TF activity found in CARNIVAL results.
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @noRd
get_TF_activity_from_CARNIVAL <- function(carnival_result, TFs)
{
    
    if(!validate_CARNIVAL_results(carnival_result)) warning("we failed to validate CARNIVAL results.")
    
    estimated_activity <- tibble::as_tibble(carnival_result$nodesAttributes) %>%
        dplyr::select(.data$Node,.data$AvgAct) %>%
        dplyr::filter(.data$AvgAct!=0) %>%
        dplyr::filter(.data$Node %in% TFs) 
    
    activity = as.numeric(estimated_activity$AvgAct)/100
    names(activity) = estimated_activity$Node
    
    return(activity)
}

#' Validate CARNIVAL Results
#' 
#' Implement here basic checks to see if we got back from CARNIVAL what we
#' expected.  Subject to change if CARNIVAL results changes.
#' 
#' @param CR CARNIVAL result object
#' @noRd
validate_CARNIVAL_results <- function(CR){
 
    list_names = c("weightedSIF","nodesAttributes","sifAll","attributesAll" ) %in% names(CR)
    if(!all(list_names)) return(FALSE)
    
    if(is.null(CR$weightedSIF)) return(FALSE)
    if(is.null(CR$nodesAttributes)) return(FALSE)
    if(is.null(CR$sifAll)) return(FALSE)
    if(is.null(CR$attributesAll)) return(FALSE)
    
    return(TRUE)
}
