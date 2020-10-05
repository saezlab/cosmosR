#' runCARNIVAL_wrapper
#' 
#' checks and formats the COSMOS data to CARNIVAL inputs and runs CARNIVAL. 
#' 
#' @param network prior knowledge network. data.frame object with source, sign and 
#' target columns
#' @param input_data numerical vector, where names are input nodes 
#' in the PKN and values are from \{1, 0, -1\}.
#' @param measured_data numerical vector, where names are measured nodes 
#' in the PKN and values are continuous values. These values are compared to 
#' with the the simulation
#' @param solver_path argument passed to \code{\link{CARNIVAL::runCARNIVAL}}
#' @param solver argument passed to \code{\link{CARNIVAL::runCARNIVAL}}
#' @param time_limit argument passed to \code{\link{CARNIVAL::runCARNIVAL}}
#' @import CARNIVAL
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
    netObj = meta_network
    
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
    
    if(!validate_CARNIVAL_results(CARNIVAL_Result)) warning("we failed to validate CARNIVAL results.")
    
    return(CARNIVAL_Result)
}