#' convert gene symbols to entrez id
#'
#' @param symbols vector of genesymbols
#' @return data.frame with human gene ENTREZID and SYMBOL mapping
#' @export
#' @seealso \code{\link{convert_ensembl_to_entrezid}}
#' @examples 
#' symbols <- c("MDH1", "PARP1", "IL6")
#' symbol_entrez_map <- convert_genesymbols_to_entrezid(symbols)
convert_genesymbols_to_entrezid <- function(symbols){

    map_table <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
    
    stopifnot(length(map_table) == length(symbols))
    stopifnot(all(names(map_table) == symbols))
    
    map_table
}