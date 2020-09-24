#' load_tf_regulon_dorothea_omnipath
#'
#' downloads the TF-regulons from Omnipaht. Different from DOROTHEA regulon because
#' sign is handled differently: if both stimulation and inhibition was reported then
#' it is removed. 
#'
load_tf_regulon_dorothea_omnipath <- function()
{
    
    url = 'http://omnipathdb.org/interactions?datasets=tfregulons&tfregulons_levels=A,B,C&genesymbols=1&fields=sources,tfregulons_level'
    dorothea <- read.table(url, sep = '\t', header = TRUE)
    
    dorothea <- dorothea[,c(4,3,6,7)]
    dorothea$sign <- dorothea$is_stimulation - dorothea$is_inhibition
    dorothea <- dorothea[dorothea$sign != 0,]
    dorothea <- dorothea[,c(1,2,5)]
    
    ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    
    G_list <- biomaRt::getBM(filters = "hgnc_symbol",
                             attributes = c('hgnc_symbol','entrezgene_id', "description"),
                             values = unique(dorothea$source_genesymbol), mart = ensembl)
    
    gene_mapping <- as.character(G_list[,2])
    names(gene_mapping) <- G_list[,1]
    
    for(i in 1:length(dorothea[,1]))
    {
        dorothea[i,2] <- gene_mapping[dorothea[i,2]]
    }
    
    dorothea$source_genesymbol <- paste0("X",dorothea$source_genesymbol)
    
    return(dorothea)
}


