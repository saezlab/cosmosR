#' Check Inputs For CARNIVAL
#' 
#' Checks the format of the main data inputs for CARNIVL.  Checks the data format
#' and coverage of nodes in the PKN and data.  Ensures that all nodes in 
#' input_data and measured_data apper in the PKN.
#' 
#' @param meta_network Prior knowledge network (PKN).  \dQuote{\code{data.frame}} 
#'   object with source, sign and target columns.  By default COSMOS uses a PKN 
#'   derived from Omnipath, STITCHdb and Recon3D.  See details on the data 
#'   \code{\link{meta_network}}.
#' 
#' @param input_data Numerical vector, where names are input nodes in the PKN 
#'   and values are from \{1, 0, -1\}.
#' 
#' @param measured_data Numerical vector, where names are measured nodes in the
#'   PKN and values are continuous values.  These values are compared to with
#'   the simulation.
#' 
#' @seealso \code{\link{default_CARNIVAL_options}}
#' 
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