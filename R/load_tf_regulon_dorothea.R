#' load transcription factor regulon
#' 
#' load the transcription factors from \code{DOROTHEA} package and converts 
#' gene symbols to EntrezID using org.Hs.eg.db
#' 
#' @param toEntrez if TRUE (default), converts gene symbols to EntrezID
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
load_tf_regulon_dorothea <- function(toEntrez = TRUE, confidence = c("A","B","C")){
    . <- NULL
    
    # load regulon from dorothea:
    regulon = dorothea::dorothea_hs
    regulon <- regulon %>% dplyr::rename(sign = "mor")
    
    conf = confidence
    regulon <- regulon %>% dplyr::filter(.data$confidence %in% conf)
    
    # transform names to EntrezID
    if(toEntrez){
        
        tf_table <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, regulon$tf, 'ENTREZID', 'SYMBOL')
        regulon$tf_entrez = tf_table
        
        target_table <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, regulon$target, 'ENTREZID', 'SYMBOL')
        regulon$target_entrez = target_table
       
        regulon <- regulon %>% dplyr::filter(stats::complete.cases(.)) %>%
            dplyr::mutate(tf = paste0("X",.data$tf_entrez),
                   target = paste0("X",.data$target_entrez)) %>%
            dplyr::select(.data$tf,.data$sign,.data$target)

    }
    
    
    return(regulon)
}



