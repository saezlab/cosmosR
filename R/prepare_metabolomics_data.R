#' format COSMOS metabolic input
#'
#' This function prepares the metabolic data to be used in the COSMOS 
#' optimization steps. It takes as input a vector with the metabolic
#'  data (e.g, limma t values) named with PUBCHEM IDs and expand it to the 
#'  multi-compartment COSMOS format. 
#'  It also messages the number of final inputs in the meta network.
#'
#' @param metabolic_data A named numeric vector, containing the values to be
#'  used for the metabolic layer in COSMOS. The names of the vector 
#'  should be PUBCHEM IDs.
#' @param meta_network Prior knowledge network
#'  created with \code{\link{load_meta_pkn()}}.
#'
#' @return A new vector ready to be used as COSMOS input.
#' 
#' @import dplyr 
#' @importFrom stringr str_extract
#' @export
prepare_metabolomics_data <- function(metabolic_data, meta_network) {
    
    # if required, format PUBCHEM IDs to match COSMOS nodes
    if( ! all( grepl("XMetab__", names(metabolic_data) ))) {
        
        message("COSMOS: Adding `XMetab__` label to metabolic data names", length(metabolic_data))
        names(metabolic_data) <- paste0("XMetab__", names(metabolic_data))
        
    }
    
    # collect all nodes in the PKN
    allPknNodes <- unique( c(dplyr::pull(meta_network, source),dplyr::pull(meta_network, target)) ) %>%
        .[grepl("XMetab__", .)]
    
    # extract all compartments in the PKN
    compartments <- stringr::str_extract(string = allPknNodes, pattern = "___.____") %>%
        unique() %>%
        .[ !is.na(.) ]
    
    # create data frame with all the possible combinations and restrict to IDs in the PKN
    nameToId <- expand.grid(names(metabolic_data), compartments, stringsAsFactors = FALSE) %>%
        as.data.frame() %>%
        dplyr::rename(name = Var1, compartment = Var2) %>%
        dplyr::mutate(nodeId = paste0(name, compartment)) 
    
    # create the expanded version of the inputs
    expandedInputs <- metabolic_data[nameToId$name]
    names(expandedInputs) <- nameToId$nodeId
    
    # message mapping info
    message("COSMOS: Original number of metabolic inputs = ", length(metabolic_data) )
    message("COSMOS: Resulting number of expanded metabolic inputs: ", length(expandedInputs) )
    
    return(expandedInputs)
    
}
