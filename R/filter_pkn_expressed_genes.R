#' filter_pkn_expressed_genes
#'
#' filters the non-expressed genes from the prior knowledge network
#'
#' @param expressed_genes_entrez  Gene_symboles of expressed genes (mut be same ID as used in meta_PKN)
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

#' filter_pkn_expressed_genes_fast
#'
#' filters the non-expressed genes from the prior knowledge network 
#' USING A FASTER OPTIMISED ALGORITHM
#'
#' @param expressed_genes_entrez  Gene_symboles of expressed genes (mut be same ID as used in meta_PKN)
#' @param meta_pkn  COSMOS prior knowledge network
#' @noRd
filter_pkn_expressed_genes_fast <- function(expressed_genes_entrez, meta_pkn) {
  message("COSMOS: removing unexpressed nodes from PKN...")
  
  # Create an environment to act as a hash table for fast membership checks.
  expressed_env <- new.env(hash = TRUE, parent = emptyenv())
  for (gene in expressed_genes_entrez) {
    expressed_env[[gene]] <- TRUE
  }
  
  vectorized_is_expressed <- function(x, env) {
    out <- rep(NA_character_, length(x))
    
    # 1. Nodes that contain "Metab" or "orphanReac" are kept as is.
    keep_idx <- grepl("Metab|orphanReac", x, perl = TRUE)
    out[keep_idx] <- x[keep_idx]
    
    # 2. Direct membership check: if the node exists in the environment.
    not_checked <- which(is.na(out))
    if (length(not_checked) > 0) {
      direct <- sapply(x[not_checked], function(s) exists(s, envir = env, inherits = FALSE))
      out[not_checked[direct]] <- x[not_checked[direct]]
    }
    
    # 3. For nodes that follow the pattern "Gene...":
    not_checked <- which(is.na(out))
    if (length(not_checked) > 0) {
      # Pattern for concatenated genes (e.g. "Gene1234__ID1_ID2")
      pattern_concat <- grepl("^Gene[0-9]+__[A-Z0-9_]+$", x[not_checked], perl = TRUE)
      if (any(pattern_concat)) {
        idx <- not_checked[pattern_concat]
        out[idx] <- sapply(x[idx], function(s) {
          genes_str <- sub("^Gene[0-9]+__", "", s)
          genes <- unlist(strsplit(genes_str, "_"))
          if (all(sapply(genes, exists, envir = env, inherits = FALSE))) s else NA_character_
        })
      }
      
      # Pattern for reverse-tagged genes (e.g. "Gene1234__ID1_ID2_reverse")
      not_checked <- which(is.na(out))
      pattern_reverse <- grepl("^Gene[0-9]+__[A-Z0-9_]+_reverse$", x[not_checked], perl = TRUE)
      if (any(pattern_reverse)) {
        idx <- not_checked[pattern_reverse]
        out[idx] <- sapply(x[idx], function(s) {
          genes_str <- sub("^Gene[0-9]+__", "", s)
          genes_str <- sub("_reverse$", "", genes_str)
          genes <- unlist(strsplit(genes_str, "_"))
          if (all(sapply(genes, exists, envir = env, inherits = FALSE))) s else NA_character_
        })
      }
      
      # Special pattern with a lower-case letter in the tag:
      not_checked <- which(is.na(out))
      pattern_lower <- grepl("^Gene[0-9]+__[^_][a-z]", x[not_checked], perl = TRUE)
      if (any(pattern_lower)) {
        idx <- not_checked[pattern_lower]
        # Optionally, print messages for these cases.
        for (s in x[idx]) message("Note: ", s)
        out[idx] <- x[idx]
      }
    }
    
    return(out)
  }
  
  meta_pkn$source <- vectorized_is_expressed(meta_pkn$source, expressed_env)
  meta_pkn$target <- vectorized_is_expressed(meta_pkn$target, expressed_env)
  
  initial_n <- nrow(meta_pkn)
  meta_pkn <- meta_pkn[stats::complete.cases(meta_pkn), ]
  removed <- initial_n - nrow(meta_pkn)
  message("COSMOS: ", removed, " interactions removed")
  
  return(meta_pkn)
}

