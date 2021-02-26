#' default_CARNIVAL_options
#' 
#' generates default CARNIVAL options.  
#' 
#' @return returns a list with all possible options implemented in CARNIVAL.
#' see the documentation on \code{\link[CARNIVAL]{runCARNIVAL}}.
#' @export
#' 
default_CARNIVAL_options = function(){
    
    opts <- list(solverPath=NULL,
         solver=c("cplex"), 
         timelimit=3600, 
         mipGAP=0.05,
         poolrelGAP=0.0001,
         limitPop=500, 
         poolCap=100,
         poolIntensity=4,
         poolReplace=2,
         alphaWeight=1, 
         betaWeight=0.2,
         threads = 1,
         dir_name=NULL
    )
    return(opts)
}



#' check_CARNIVAL_options
#' 
#' checks options for CARNIVAL
#' 
check_CARNIVAL_options <- function(opts){
    
    
    if(!is.list(opts)) stop("CARNIVAL options should be a list")
    req_names <- c(
        "solverPath",
        "solver", 
        "timelimit",
        "mipGAP",
        "poolrelGAP",
        "limitPop", 
        "poolCap",
        "poolIntensity",
        "poolReplace",
        "alphaWeight", 
        "betaWeight",
        "threads",
        "dir_name")
    
    if(!all(req_names %in% names(opts))){
        stop("CARNIVAL options should contain all options. 
             Start by calling default_CARNIVAL_options() and replace entries. ")
    }
    
    
    if(!(opts$solver %in% c("cplex","cbc"))) stop("current version supports only the CPLEX or cbc solver")
    if(opts$solver == "cbc") print("COSMOS wasn't tested thoroughly with the cbc solver. We recommend the users to use CPLEX if possible, and use cbc as a backup solution.")
    if(is.null(opts$solverPath)) stop("path to CPLEX solver must be provided")
    
    
    
}
