

#' load transcription factor regulon
#' 
#' load the transcription factors from \code{DOROTHEA} package and converts 
#' gene symbols to EntrezID using biomart
#' 
#' @param toEntrez if TRUE (default), converts gene symbols to EntrezID
#' @import dorothea
#' @importFrom biomaRt getBM useEnsembl
load_tf_regulon_dorothea <- function(toEntrez = TRUE, confidence = c("A","B","C")){
    
    # load regulon from dorothea:
    regulon = dorothea::dorothea_hs
    regulon <- regulon %>% rename(sign = "mor")
    
    # transform names to EntrezID
    if(toEntrez){
        
        
        if(!dir.exists("./cosmos_cache")) dir.create("./cosmos_cache")

        
        Biobase::cache(
            ensembl <- biomaRt::useEnsembl(biomart="ensembl",
                                           dataset="hsapiens_gene_ensembl"),
            dir = "./cosmos_cache")
        Biobase::cache(
            tf_entrez <- biomaRt::getBM(filters = "hgnc_symbol",
                                        attributes = c('hgnc_symbol',
                                                       'entrezgene_id'),
                                        values = unique(regulon$tf)[1:4],
                                        mart = ensembl),
            dir = "./cosmos_cache")
        
        
        
        regulon <- regulon %>% dplyr::left_join(tf_entrez, by=c("tf"="hgnc_symbol" )) %>%
            mutate(tf = paste0("X",entrezgene_id)) %>% 
            dplyr::select(-entrezgene_id)
        
        Biobase::cache( target_entrez <- biomaRt::getBM(filters = "hgnc_symbol",
                                        attributes = c('hgnc_symbol','entrezgene_id'),
                                        values = unique(regulon$target), mart = ensembl),
                        dir = "./cosmos_cache")
        
        regulon <- regulon %>% dplyr::left_join(target_entrez, by=c("target"="hgnc_symbol" )) %>%
            mutate(target = paste0("X",entrezgene_id)) %>% 
            dplyr::select(-entrezgene_id)
        
        
        
        regulon <- regulon %>% dplyr::filter(complete.cases(.))
    }
    conf = confidence
    
    regulon <- regulon %>% dplyr::filter(confidence %in% conf)
    
    return(regulon)
}



