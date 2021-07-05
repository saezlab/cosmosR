#' OmniPath PTMs 
#' 
#' Collection of annotated enzyme-substrate post translational modifications
#' obtained from OmniPath.
#'
#' @docType data
#'
#' @usage data(omnipath_ptm)
#'
#' @format An object of class \dQuote{\code{data.frame}} with 39201 rows (PTMs) 
#'   and 12 variables:
#'   \describe{
#'     \item{\code{enzyme}}{}
#'     \item{\code{substrate}}{}
#'     \item{\code{enzyme_genesymbol}}{}
#'     \item{\code{substrate_genesymbol}}{}
#'     \item{\code{residue_type}}{}
#'     \item{\code{residue_offset}}{}
#'     \item{\code{modification}}{}
#'     \item{\code{sources}}{}
#'     \item{\code{references}}{}
#'     \item{\code{curation_effort}}{}
#'     \item{\code{n_references}}{}
#'     \item{\code{n_resources}}{}
#'   }
#'   
#' @source Default resource collection of OmniPath:
#'   \url{http://omnipathdb.org/ptms?fields=sources,references&genesymbols=1},
#'   version of Feb 5th, 2020.
#'
#' @references {
#'   Turei, D., Korcsmaros, T. and Saez-Rodriguesz, J. (2016) \emph{Nature Methods}.
#'   \bold{13}(12), 966--967.
#' }
#' 
#' @examples
#' data(omnipath_ptm)
#' 
#' @keywords datasets
"omnipath_ptm"