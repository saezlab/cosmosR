#' format COSMOS signaling input
#'
#' @param signaling_data A named numeric vector, containing the values to be used for the signaling layer in COSMOS. 
#' @param pkn Prior knowledge network created with \code{\link{load_meta_pkn()}}.
#' @param format_names Add "X" to the vector names?
#'
#' @return A new vector ready to be used as COSMOS input.
#' 
#' @import dplyr 
#' 
prepare_signaling_data <- function(signaling_data, pkn, format_names = TRUE) {
    
    # if required, format ENTREZ IDs to match COSMOS nodes
    if(format_names) {
        
        names(signaling_data) <- paste0("X", names(signaling_data))
        
    }
    
    # collect all nodes in the PKN
    allPknNodes <- unique( c(dplyr::pull(pkn, source), dplyr::pull(pkn, target)) )
    
    # filter inputs to PKN and message info
    inPkn <- names(signaling_data) %in% allPknNodes
    message("Original number of inputs: ", length(signaling_data))
    message("Number of inputs in the network: ", length(inPkn))
    
    # return filtered Data
    filteredData <- signaling_data[inPkn]
    
    return(signaling_data)
    
}

#' format COSMOS metabolic input
#'
#' @param metab_data A named numeric vector, containing the values to be used for the metabolic layer in COSMOS. 
#' @param pkn Prior knowledge network created with \code{\link{load_meta_pkn()}}.
#' @param format_names Add "XMetab" to the vector names?
#'
#' @return A new vector ready to be used as COSMOS input.
#' 
#' @import dplyr 
#' @importFrom stringr str_extract
#'
prepare_metab_data <- function(metab_data, pkn, format_names = TRUE) {
    
    # if required, format PUBCHEM IDs to match COSMOS nodes
    if(format_names) {
        
        names(metab_data) <- paste0("XMetab__", names(metab_data))
        
    }
    
    # collect all nodes in the PKN
    allPknNodes <- unique( c(dplyr::pull(pkn, source),dplyr::pull(pkn, target)) ) %>%
        .[grepl("XMetab__", .)]
    
    # extract all compartments in the PKN
    compartments <- stringr::str_extract(string = allPknNodes, pattern = "___.____") %>%
        unique() %>%
        .[ !is.na(.) ]
    
    # create data frame with all the possible combinations and restrict to IDs in the PKN
    nameToId <- expand.grid(names(metab_data), compartments, stringsAsFactors = FALSE) %>%
        as.data.frame() %>%
        dplyr::rename(name = Var1, compartment = Var2) %>%
        dplyr::mutate(nodeId = paste0(name, compartment)) %>%
        subset(nodeId %in% allPknNodes)
    
    # create the expanded version of the inputs
    expandedInputs <- metab_data[nameToId$name]
    names(expandedInputs) <- nameToId$nodeId
    
    # message mapping info
    message("Original number of inputs: ", length(metab_data))
    message("Number of inputs in the network: ", length(unique(nameToId$name)))
    message("Number of inputs in the expanded version across compartments: ", length(expandedInputs))
    
    return(expandedInputs)
    
}
