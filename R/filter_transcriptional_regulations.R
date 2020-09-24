
#' filter_transcriptional_regulations
#' 
#' Filters interactions between TFs and genes in the PKN by the following points
#' - finds the TF target genes that change strongly (threshold)
#' - finds TFs from the inputs (known TF activity)
#' - remove interactions if the target gene expression doesn't change (regulation 
#' should not go through this node )
#' - remove interactions if the TF activity is not in agreement with the change
#'  in the target gene expression (other factors influence these genes)
#'  
#'  @param network meta PKN
#'  @param gene_expression_binarized named vector of {-1,0,1} difnining which genes changed
#'  @param signaling_data named vector containing known activity of signaling nodes (kinases, TFs)
#'  @param tf_regulon transcription factor regulon, e.g. from dorothea
#'  @return a filtered version of the network 
filter_transcriptional_regulations <- function(network, 
                                               gene_expression_binarized,
                                               signaling_data,
                                               tf_regulon)
{
    check_COSMOS_inputs(meta_network = network,
                        tf_regulon = tf_regulon,
                        signaling_data = signaling_data,
                        expression_data = gene_expression_binarized)
    
    # map the TF values and gene expression levels on the PKN
    annotated_network <- network %>% 
        # add the gene expression data: (non-genes will be filled with NA)
        left_join(
            enframe(gene_expression_binarized, name = "gene", value="target_sign"),
            by=c(target="gene")) %>%
        # add the TF data (non-TFs will be filled with NAs)
        left_join(
            enframe(signaling_data, name = "TF", value="TF_sign") %>%
                filter(TF %in% tf_regulon$tf),
            by=c(source="TF")) %>%
        mutate(source_is_TF = source %in% tf_regulon$tf)
    
    # find interactions where TF regulates a gene, but target 
    # - does not change 
    # - target gene is not measured
    # - target changes inconsistently
    annotated_network <- annotated_network %>% 
        # genes didn't change
        mutate(target_gene_unchanged = target_sign==0) %>%
        # gene didnt change and the source is a TF (it is a transcriptional regulation)
        mutate(TF_target_unchanged = source_is_TF &
                   (target_gene_unchanged | is.na(target_gene_unchanged))) %>%
        # gene changed, but not consistently with TF activity
        mutate(inconsistent_TF_gene_sign = 
                   sign(TF_sign) != interaction * sign(target_sign) ) 
    
    # TODO: option to return these interactions
    removed_interactions <- annotated_network %>%
        filter(TF_target_unchanged | inconsistent_TF_gene_sign)
    
    print(paste("COSMOS: ", nrow(removed_interactions), 
                "interactions are removed from the PKN based on",
                "consistency check between TF activity and gene expression"))
    
    
    kept_interactions <-  annotated_network %>%
        filter(!TF_target_unchanged | is.na(TF_target_unchanged)) %>%
        filter(!inconsistent_TF_gene_sign | is.na(inconsistent_TF_gene_sign)) 
    
    out_pkn <- kept_interactions %>% select(source,interaction,target)    
    return(out_pkn)
}


