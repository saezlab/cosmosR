#' format COSMOS signaling input
#'
#' @param signalingData A named numeric vector, containing the values to be used for the signaling layer in COSMOS. 
#' @param pkn Prior knowledge network created with \code{\link{load_meta_pkn()}}.
#' @param formatNames Add "X" to the vector names?
#'
#' @return A new vector ready to be used as COSMOS input.
#' 
#' @import dplyr 
#' 
prepareSignalingData <- function(signalingData, pkn, formatNames = TRUE) {
    
    # if required, format ENTREZ IDs to match COSMOS nodes
    if(formatNames) {
        
        names(signalingData) <- paste0("X", names(signalingData))
        
    }
    
    # collect all nodes in the PKN
    allPknNodes <- unique( c(dplyr::pull(pkn, source),dplyr::pull(pkn, target)) )
    
    # filter inputs to PKN and message info
    inPkn <- names(signalingData) %in% allPknNodes
    message("Original number of inputs: ", length(signalingData))
    message("Number of inputs in the network: ", length(inPkn))
    
    # return filtered Data
    filteredData <- signalingData[inPkn]
    
    return(signalingData)
    
}

#' format COSMOS metabolic input
#'
#' @param metabData A named numeric vector, containing the values to be used for the metabolic layer in COSMOS. 
#' @param pkn Prior knowledge network created with \code{\link{load_meta_pkn()}}.
#' @param formatNames Add "X" to the vector names?
#'
#' @return A new vector ready to be used as COSMOS input.
#' 
#' @import dplyr 
#' @import stringr
#'
prepareMetabData <- function(metabData, pkn, formatNames = TRUE) {
    
    # if required, format PUBCHEM IDs to match COSMOS nodes
    if(formatNames) {
        
        names(metabData) <- paste0("XMetab__", names(metabData))
        
    }
    
    # collect all nodes in the PKN
    allPknNodes <- unique( c(dplyr::pull(pkn, source),dplyr::pull(pkn, target)) ) %>%
        .[grepl("XMetab__", .)]
    
    # extract all compartments in the PKN
    compartments <- stringr::str_extract(string = allPknNodes, pattern = "___.____") %>%
        unique() %>%
        .[ !is.na(.) ]
    
    # create data frame with all the possible combinations and restrict to IDs in the PKN
    nameToId <- expand.grid(names(metabData), compartments, stringsAsFactors = FALSE) %>%
        as.data.frame() %>%
        dplyr::rename(name = Var1, compartment = Var2) %>%
        dplyr::mutate(nodeId = paste0(name, compartment)) %>%
        subset(nodeId %in% allPknNodes)
    
    # create the expanded version of the inputs
    expandedInputs <- metabData[nameToId$name]
    names(expandedInputs) <- nameToId$nodeId
    
    # message mapping info
    message("Original number of inputs: ", length(metabData))
    message("Number of inputs in the network: ", length(unique(nameToId$name)))
    message("Number of inputs in the expanded version across compartments: ", length(expandedInputs))
    
    return(expandedInputs)
    
}
