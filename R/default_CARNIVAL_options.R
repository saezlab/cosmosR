#' Setting Default CARNIVAL Options
#' 
#' Returns the default CARNIVAL options as a list.  You can modify the elements
#' of the list and then use it as an argument in \dQuote{\code{run_COSMOS()}}.
#' 
#' @param solver one of `cplex` (recommended, but require 3rd party tool), `cbc` (also require 3rd party tool) or `lpSolve` (default, only for small networks) 
#' @param verbose (default: TRUE)  notify the user to include path for CPLEX or CBC path
#' @return returns a list with all possible options implemented in CARNIVAL.
#' see the documentation on \code{\link[CARNIVAL]{runCARNIVAL}}.
#' @examples 
#' # load and change default options: 
#' my_options = default_CARNIVAL_options(solver = "cplex")
#' 
#' my_options$solverPath = "/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex"
#' my_options$threads = 2
#' my_options$timelimit = 3600*15
#' @export
#' 
default_CARNIVAL_options = function(solver = c("lpSolve","cplex","cbc"),verbose = TRUE){
    
    solver <- match.arg(solver)
    
    if(solver == "lpSolve"){
        opts = CARNIVAL::defaultLpSolveCarnivalOptions()
    }else if(solver == "cplex"){
        opts = CARNIVAL::defaultCplexCarnivalOptions()
        if(verbose){
            warning("you selected CPLEX solver. Make sure to set the solverPath field to cplex executable.")
        }
    }else if(solver == "cbc"){
        opts = CARNIVAL::defaultCbcSolveCarnivalOptions()
        if(verbose){
            warning("you selected CBC solver. Make sure to set the solverPath field to CBC executable.")
        }
    }
    opts$keepLPFiles = FALSE
    
    return(opts)
}



#' Check CARNIVAL Options
#' 
#' Checks the list of input options for CARNIVAL for completeness.
#' 
#' @param opts An object of type \dQuote{\code{list}} defining the run parameters
#'   for CARNIVAL.  Use the \code{\link{default_CARNIVAL_options}} function to 
#'   create a list with default parameter settings. If cplex or cbc are chosen 
#'   as the solver, the parameter solverPath needs to be supplied (not 
#'   automatically added by \code{default_CARNIVAL_options()}).
#'   
#' @seealso \code{\link{default_CARNIVAL_options}},  
#'   \code{\link[CARNIVAL]{runCARNIVAL}}
#' @noRd
check_CARNIVAL_options <- function(opts){
    
    
    if(!is.list(opts)) stop("CARNIVAL options should be a list")
    req_names <- c(
        # "solverPath",
        "solver" 
        #"timelimit"
        # "mipGAP",
        # "poolrelGAP",
        # "limitPop", 
        # "poolCap",
        # "poolIntensity",
        # "poolReplace",
        # "alphaWeight", 
        # "betaWeight",
        # "threads",
        # "dir_name"
        )
    
    # if(opts$solver == "lpSolve")
    # {
    #     req_names <- req_names[-which(req_names == "solverPath")]
    # }
    # 
    if(!all(req_names %in% names(opts))){
        stop("CARNIVAL options should contain all options. 
             Start by calling CARNIVAL::defaultCbcCarnivalOptions() for CPLEX  or
             CARNIVAL::defaultCplexCarnivalOptions() for CBC solver, or
             CARNIVAL::defaultLpSolveCarnivalOptions() for LpSolve and replace entries. ")
    }
    
    
    if(!(opts$solver %in% c("cplex","cbc","lpSolve"))) stop("current version supports only the CPLEX or cbc solver. lpSolve is also available for test runs.")
    if(opts$solver == "cbc") print("COSMOS wasn't tested thoroughly with the cbc solver. We recommend the users to use CPLEX if possible, and use cbc as a backup solution.")
    if(opts$solver == "lpSolve") print("lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers.")
    if(is.null(opts$solverPath) & opts$solver != "lpSolve") stop("path to CPLEX or cbc solver must be provided")
    
    return('List of CARNIVAL options is complete.')
    
    
}
