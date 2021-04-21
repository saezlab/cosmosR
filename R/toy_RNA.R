#' Toy Input Transcription Data Set
#' 
#' This transcription data are a subset from the differential expression results
#' analysed in the case study of the COSMOS paper.  The subset consists of all 
#' genes present in the toy network.
#'
#' @docType data
#'
#' @usage toy_RNA
#'
#' @format An object of class \dQuote{\code{numeric}} containing the t-values of
#'   385 genes, which are named with gene Entrez IDs matching the toy network.
#'
#' @source 
#'   \url{https://github.com/saezlab/COSMOS_MSB/blob/main/data/RNA_ttop_tumorvshealthy.csv}
#' 
#' @references {
#'   Dugourd, A., Kuppe, C. and Sciacovelli, M. et. al. (2021) \emph{Molecular 
#'   Systems Biology}. \bold{17}, e9730.
#' }
#' 
#' @examples
#' toy_transcription_input = toy_RNA
#' 
#' @keywords datasets
"toy_RNA"