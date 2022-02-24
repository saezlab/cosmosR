#' Toy Signaling Input
#' 
#' This signaling data are a subset of the footprint-based signaling activity
#' estimates of transcription factors of the 786-O cell line from the NCI60 dataset.
#'
#' @docType data
#'
#' @usage data(toy_signaling_input)
#'
#' @format  An object of class \dQuote{\code{data.frame}} containing the normalised
#'   enrichment scores (NES) of 2 signaling proteins, which are named with their
#'   respective gene Entrez ID matching the toy network.
#'
#' @source Subset of: 
#'   \url{https://github.com/saezlab/COSMOS_MSB/blob/main/data/signaling_input_COSMOS.csv}
#'   
#' @references {
#'   Dugourd, A., Kuppe, C. and Sciacovelli, M. et. al. (2021) \emph{Molecular 
#'   Systems Biology}. \bold{17}, e9730.
#' }
#' 
#' @examples
#' data(toy_signaling_input)
#' 
#' @keywords datasets
"toy_signaling_input"