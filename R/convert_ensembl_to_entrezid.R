#' convert gene ensembl to entrez id
#'
#' @param ensembl vector of genes with ensembl id
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @export
#' @seealso [convert_genesymbols_to_entrezid()]
convert_ensembl_to_entrezid <- function(ensembl){
    
    require(org.Hs.eg.db)
    require(AnnotationDbi)
    
    map_table <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, ensembl, 'ENTREZID', 'ENSEMBL')
    
    stopifnot(length(map_table) == length(ensembl))
    stopifnot(all(names(map_table) ==ensembl))
    
    map_table
}