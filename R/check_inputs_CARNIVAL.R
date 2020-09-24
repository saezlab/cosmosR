#' check_inputs_for_CARNIVAL
#' 
#' checks the format of the main data inputs for CARNIVL. Checks the data format
#' and coverage of nodes in the PKN and data. All nodes in input_data and 
#' measured data must apper in the PKN 
check_inputs_for_CARNIVAL <- function(meta_network,
                                      input_data,
                                      measured_data){
    # checking the data
    stopifnot(is.vector(input_data))
    stopifnot(is.vector(measured_data))
    
    # checking the networks
    
    stopifnot(is.data.frame(meta_network))
    stopifnot(all(c("source","interaction","target" ) %in% names(meta_network)))
    stopifnot(ncol(meta_network)==3)
    
    # check inputs and measurements are in the network
    stopifnot(all(names(input_data) %in% c(meta_network$source, meta_network$target)))
    stopifnot(all(names(measured_data) %in% c(meta_network$source, meta_network$target)))
    
    return(TRUE)
}