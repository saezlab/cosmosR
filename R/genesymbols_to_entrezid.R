#' convert gene symbols to entrez id
#'
#' @param symbols vector of genesymbols 
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @export
genesymbols_to_entrezid <- function(symbols){
    
    require(org.Hs.eg.db)
    require(AnnotationDbi)
    
    map_table <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
    
    stopifnot(length(map_table) == length(symbols))
    stopifnot(all(names(map_table) == symbols))
    
    map_table
}