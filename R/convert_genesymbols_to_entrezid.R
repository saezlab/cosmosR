#' convert gene symbols to entrez id
#'
#' @param symbols vector of genesymbols 
#' @export
#' @seealso [convert_ensembl_to_entrezid()]
convert_genesymbols_to_entrezid <- function(symbols){

    map_table <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
    
    stopifnot(length(map_table) == length(symbols))
    stopifnot(all(names(map_table) == symbols))
    
    map_table
}