#' load transcription factor regulon
#' 
#' load the transcription factors from \code{DOROTHEA} package and converts 
#' gene symbols to EntrezID using org.Hs.eg.db
#' 
#' @param confidence strong vector (by default: c("A","B","C")). Subset of \{A, B,
#'  C, D, E\}. See the `dorothea` for the meaning of confidence levels. 
#' package for further details. 
#' @return returns a PKN of a form of a data table. Each row is an interaction.
#' Columns names are:
#' 
#' - `tf` transcription factor
#' - `confidence` class of confidence
#' - `target` target gene
#' - `sign` indicates if interaction is up (1) or down-regulation (-1). 
#' 
#' @importFrom magrittr %>%
#' @export
#' @examples 
#' load_tf_regulon_dorothea()
load_tf_regulon_dorothea <- function(confidence = c("A","B","C")){
    . <- NULL
    
    # load regulon from dorothea:
    regulon = dorothea::dorothea_hs
    regulon <- regulon %>% dplyr::rename(sign = "mor")
    
    conf = confidence
    regulon <- regulon %>% dplyr::filter(.data$confidence %in% conf)
    regulon <- regulon[,-which(names(regulon) == "confidence")]
    
    return(regulon)
}



