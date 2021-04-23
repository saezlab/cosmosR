#' Run CARNIVAL Wrapper
#' 
#' Checks and formats the COSMOS data to CARNIVAL inputs and runs CARNIVAL. 
#' 
#' @param network Prior knowledge network (PKN).  \dQuote{\code{data.frame}} 
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
#' @param options An object of type \dQuote{\code{list}} defining the run 
#'   parameters for CARNIVAL.  Use the \code{\link{default_CARNIVAL_options}}
#'   function to create a list with default parameter settings. If cplex or cbc 
#'   are chosen as the solver, the parameter solverPath needs to be supplied 
#'   (not automatically added by \code{default_CARNIVAL_options()}).
#' @noRd
runCARNIVAL_wrapper <- function(network, 
                                input_data,
                                measured_data,
                                options
                                ){
    
    check_CARNIVAL_options(options)
    check_inputs_for_CARNIVAL(meta_network = network,
                              input_data = input_data,
                              measured_data = measured_data)
    
    
    inputObj = dplyr::bind_rows(input_data)
    measObj  = dplyr::bind_rows(measured_data)
    # netObj = meta_network
    netObj = network
    
    CARNIVAL_Result <- CARNIVAL::runCARNIVAL(inputObj = inputObj,
                                             measObj = measObj,
                                             netObj = netObj,
                                             poolrelGAP = options$poolrelGAP,
                                             limitPop = options$limitPop,
                                             poolCap = options$poolCap,
                                             poolIntensity = options$poolIntensity,
                                             alphaWeight = options$alphaWeight,
                                             betaWeight = options$betaWeight,
                                             poolReplace = options$poolReplace,
                                             threads = options$threads,
                                             solverPath = options$solverPath,
                                             solver = options$solver,
                                             timelimit = options$timelimit,
                                             mipGAP = options$mipGAP,
                                             dir_name = options$dir_name)
    
    if(!validate_CARNIVAL_results(CARNIVAL_Result)) warning("We failed to validate CARNIVAL results.")
    
    return(CARNIVAL_Result)
}