#' Metabolite-PubChem CID Mapping
#' 
#' Mapping between metabolite names and PubChem CIDs obtained from the recon3D
#' metabolic model.  Combined table version from the recon3D matlab object.
#'
#' @docType data
#'
#' @usage data(metabolite_to_pubchem)
#'
#' @format An object of class \dQuote{\code{data.frame}} with 1131 rows (metabolites) 
#'   and two variables:
#'   \describe{
#'     \item{\code{name}}{Metabolite name synonym}
#'     \item{\code{pubchem}}{Pubchem CID}
#'   }
#' 
#' @source \url{https://www.vmh.life/#downloadview}, downloaded on Feb 19th, 2018.
#' 
#' @references {
#'   Brunk, E. et al. (2018) \emph{Nature Biotechnology}. \bold{36}(3), 272--281.
#' }
#' 
#' @examples
#' metab_to_pub = metabolite_to_pubchem
#' 
#' @keywords datasets
"metabolite_to_pubchem"