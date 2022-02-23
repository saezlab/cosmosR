#' filter_pkn_expressed_genes
#'
#' filters the non-expressed genes from the prior knowledge network
#'
#' @param expressed_genes_entrez  EntrezID of expressed genes
#' @param meta_pkn  COSMOS prior knowledge network
#' @noRd
filter_pkn_expressed_genes <- function(expressed_genes_entrez,meta_pkn){
    
    print(paste("COSMOS: removing unexpressed nodes from PKN..."))    

    is_expressed <- function(x)
    {
        if(!grepl("Metab",x) & !grepl("orphanReac",x))
        {
            if(x %in% expressed_genes_entrez)
            {
                return(x)
            } else
            {
                if(grepl("Gene[0-9]+__[A-Z0-9_]+$",x))
                {
                    genes <- gsub("Gene[0-9]+__","",x)
                    genes <- strsplit(genes,"_")[[1]]
                    if(sum(genes %in% expressed_genes_entrez) != length(genes))
                    {
                        return(NA)
                    } else
                    {
                        return(x)
                    }
                } else
                {
                    if(grepl("Gene[0-9]+__[^_][a-z]",x))
                    {
                        print(x)
                        return(x)
                    } else
                    {
                        if(grepl("Gene[0-9]+__[A-Z0-9_]+reverse",x))
                        {
                            genes <- gsub("Gene[0-9]+__","",x)
                            genes <- gsub("_reverse","",genes)
                            genes <- strsplit(genes,"_")[[1]]
                            if(sum(genes %in% expressed_genes_entrez) != length(genes))
                            {
                                return(NA)
                            } else
                            {
                                return(x)
                            }
                        } else
                        {
                            return(NA)
                        }
                    }
                }
            }
        } else
        {
            return(x)
        }
    }
    
    # is_expressed("XGene3004__124975_91227")
    
    meta_pkn$source <- sapply(meta_pkn$source,is_expressed)
    n_removed = sum(!stats::complete.cases(meta_pkn))
    meta_pkn <- meta_pkn[stats::complete.cases(meta_pkn),]
    meta_pkn$target <- sapply(meta_pkn$target,is_expressed)
    meta_pkn <- meta_pkn[stats::complete.cases(meta_pkn),]
    n_removed = n_removed + sum(!stats::complete.cases(meta_pkn))
    print(paste("COSMOS:", n_removed, "interactions removed"))
    return(meta_pkn)
}
