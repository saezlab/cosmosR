#' convert gene ensembl to entrez id
#'
#' @param ensembl vector of genes with ensembl id
#'
#' @return named vector, where names are the old ensemblIDs and values are the 
#' entrezIDs
#' @export
#' @seealso [convert_genesymbols_to_entrezid()]
#' @examples 
#' ensembl <- c("ENSG00000100601", "ENSG00000178826", "ENSG00000138231")
#' entrez_map <- convert_ensembl_to_entrezid(ensembl)
convert_ensembl_to_entrezid <- function(ensembl){
    
    map_table <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, ensembl, 'ENTREZID', 'ENSEMBL')
    
    stopifnot(length(map_table) == length(ensembl))
    stopifnot(all(names(map_table) == ensembl))
    if(any(is.na(map_table))) message("WARNING: some ensemble IDs were not found.")
    
    map_table
}