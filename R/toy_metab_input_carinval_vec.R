#' Toy Metabolic Input Data
#' 
#' This metabolic data are a subset from the metabolic measurements used as an
#' input in the case study of the COSMOS paper.  The subset contains a random 
#' selection of metabolites present in the toy network.
#'
#' @docType data
#'
#' @usage toy_metab_input_carnival_vec
#'
#' @format An object of class \dQuote{\code{numeric}} containing the t-values of
#'   10 metabolites, which are named with metabolite PubChem CIDs matching the
#'   toy network.
#'
#' @source
#'   \url{https://github.com/saezlab/COSMOS_MSB/blob/main/data/metab_input_COSMOS.csv}
#' 
#' @references {
#'   Dugourd, A., Kuppe, C. and Sciacovelli, M. et. al. (2021) \emph{Molecular 
#'   Systems Biology}. \bold{17}, e9730.
#' }
#' 
#' @examples
#' toy_metabolic_input = toy_metab_input_carnival_vec
#' 
#' @keywords datasets
"toy_metab_input_carnival_vec"