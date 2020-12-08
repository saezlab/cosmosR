#' default_CARNIVAL_options
#' 
#' returns the default CARNIVAL options as a list. You can modify the elements of
#' the list and then use it as an argument in `run_COSMOS()`
#' 
#' @return returns a list with all possible options implemented in CARNIVAL.
#' see the documentation on \code{\link[CARNIVAL]{runCARNIVAL}}.
#' @examples 
#' # load and change default options: 
#' my_options = default_CARNIVAL_options()
#' 
#' my_options$solverPath = "/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex"
#' my_options$threads = 2
#' my_options$timelimit = 3600*15
#' # pass the options to COSMOS:
#' \dontrun{
#' run_COSMOS_metabolism_to_signaling([other inputs], CARNIVAL_options = my_options)
#' }
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
    
    
    if(opts$solver!="cplex") stop("current version supports only the CPLEX solver")
    if(is.null(opts$solverPath)) stop("path to CPLEX solver must be provided")
    
    
    
}
