#' Toy Input Network
#' 
#' This signaling network is the reduced COSMOS network solution obtained in the
#' case study of the COSMOS paper.  Here, this network solution is reused as an
#' exemplary input prior knowledge network (PKN).
#'
#' @docType data
#'
#' @usage toy_network
#'
#' @format An object of class \dQuote{\code{data.frame}} with 720 rows
#'   (interactions) and three variables:
#'   \describe{
#'     \item{\code{source}}{Source node, either metabolite or protein}
#'     \item{\code{interaction}}{Type of interaction, 1 = Activation, -1 = Inhibition}
#'     \item{\code{target}}{Target node, either metabolite or protein}
#'   A detailed description of the identifier formatting can be found under 
#'   \url{https://metapkn.omnipathdb.org/00__README.txt}.
#'   }
#'
#' @source  The network is available on github:
#'   \url{https://github.com/saezlab/COSMOS_MSB/tree/main/results/COSMOS_result}
#'   
#' @references {
#'   Dugourd, A., Kuppe, C. and Sciacovelli, M. et. al. (2021) \emph{Molecular 
#'   Systems Biology}. \bold{17}, e9730.
#' }
#' 
#' @examples
#' toy_pkn = toy_network
#' 
#' @keywords datasets
"toy_network"