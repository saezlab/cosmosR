#' Meta Prior Knowledge Network
#' 
#' Comprehensive Prior Knowledge Network (PKN), which combines signaling and
#' metabolic interaction networks.  The network was constructed using the Recon3D 
#' and STITCH metabolic networks as well as the signaling network from
#' OmniPath.
#' 
#' @docType data
#'
#' @usage meta_network
#'
#' @format An object of class \dQuote{\code{tibble}} with 117065 rows
#'   (interactions) and three variables:
#'   \describe{
#'     \item{\code{source}}{Source node, either metabolite or protein}
#'     \item{\code{interaction}}{Type of interaction, 1 = Activation, -1 = Inhibition}
#'     \item{\code{target}}{Target node, either metabolite or protein}
#'   A detailed description of the identifier formatting can be found under 
#'   \url{https://metapkn.omnipathdb.org/00__README.txt}.
#'   }
#'
#' @source The network is available in Omnipath:
#'   \url{https://metapkn.omnipathdb.org/metapkn__20200122.txt}, the scripts 
#'   used for the build of the network are available under
#'   \url{https://github.com/saezlab/Meta_PKN}.
#' 
#' @references {
#'   Dugourd, A., Kuppe, C. and Sciacovelli, M. et. al. (2021) \emph{Molecular 
#'   Systems Biology}. \bold{17}, e9730.
#' }
#' 
#' @examples
#' meta_pkn = meta_network
#' 
#' @keywords datasets
"meta_network"