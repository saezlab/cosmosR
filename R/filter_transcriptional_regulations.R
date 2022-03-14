#' Filter Based On Transcriptional Regulations
#' 
#' Filters interactions between TFs and genes in the PKN by the following points:
#' - finds the TF target genes that change strongly (threshold)
#' - finds TFs from the inputs (known TF activity)
#' - remove interactions if the target gene expression doesn't change (regulation 
#' should not go through this node )
#' - remove interactions if the TF activity is not in agreement with the change
#'  in the target gene expression (other factors influence these genes)
#'  
#' @param network Prior knowledge network (PKN).  By default COSMOS use a PKN 
#'   derived from Omnipath, STITCHdb and Recon3D. See details on the data 
#'   \code{\link{meta_network}}.
#' @param gene_expression_binarized Named vector of {-1,0,1} difnining which 
#'   genes changed.
#' @param signaling_data Named vector containing known activity of signaling
#'   nodes (kinases, TFs).
#' @param tf_regulon Collection of transcription factor - target interactions.
#'  A default collection from dorothea can be obtained by the 
#'  \code{\link{load_tf_regulon_dorothea}} function.
#' @return A filtered version of the network.
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @noRd
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
    
    gene_exp_df = data.frame(gene =names(gene_expression_binarized),
                             target_sign = gene_expression_binarized)

    signaling_df = data.frame(TF =names(signaling_data),
                              TF_sign = signaling_data) %>%
        dplyr::filter(.data$TF %in% tf_regulon$tf)
    
    annotated_network <- network %>% 
        # add the gene expression data: (non-genes will be filled with NA)
        dplyr::left_join(gene_exp_df, by=c(target="gene")) %>%
        # add the TF data (non-TFs will be filled with NAs)
        dplyr::left_join( signaling_df, by=c(source="TF")) %>%
        dplyr::mutate(source_is_TF = .data$source %in% tf_regulon$tf)
    
    # find interactions where TF regulates a gene, but target 
    # - does not change 
    # - target gene is not measured
    # - target changes inconsistently
    annotated_network <- annotated_network %>% 
        # genes didn't change
        dplyr::mutate(target_gene_unchanged = .data$target_sign==0) %>%
        # gene didnt change and the source is a TF (it is a transcriptional regulation)
        dplyr::mutate(TF_target_unchanged = .data$source_is_TF &
                   (.data$target_gene_unchanged | is.na(.data$target_gene_unchanged))) %>%
        # gene changed, but not consistently with TF activity
        dplyr::mutate(inconsistent_TF_gene_sign = 
                   sign(.data$TF_sign) != interaction * sign(.data$target_sign) ) 
    
    # TODO: option to return these interactions
    removed_interactions <- annotated_network %>%
        dplyr::filter(.data$TF_target_unchanged | .data$inconsistent_TF_gene_sign)
    
    print(paste("COSMOS: ", nrow(removed_interactions), 
                "interactions are removed from the PKN based on",
                "consistency check between TF activity and gene expression"))
    
    
    kept_interactions <-  annotated_network %>%
        dplyr::filter(!.data$TF_target_unchanged | is.na(.data$TF_target_unchanged)) %>%
        dplyr::filter(!.data$inconsistent_TF_gene_sign | is.na(.data$inconsistent_TF_gene_sign)) 
    
    out_pkn <- kept_interactions %>% dplyr::select(.data$source,.data$interaction,.data$target)    
    return(out_pkn)
}


