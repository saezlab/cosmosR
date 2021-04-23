#' Preprocess COSMOS Inputs For Signaling to Metabolism
#' 
#' Runs checks on the input data and simplifies the prior knowledge network.
#' Simplification includes the removal of (1) nodes that are not reachable from 
#' signaling nodes and  (2) interactions between transcription factors and target
#' genes if the target gene does not respond or the response is contradictory
#' with the change in the transcription factor activity. 
#' Optionally, further TF activities are estimated via network optimization via
#' CARNIVAL and the interactions between TF and genes are filtered again. 
#' 
#' @param meta_network prior knowledge network. By default COSMOS use a PKN 
#' derived from Omnipath, STITCHdb and Recon3D. See details on the data 
#' \code{\link{meta_network}}.
#' @param tf_regulon collection of transcription factor - target interactions.
#' A default collection from dorothea can be obtained by the 
#' \code{\link{load_tf_regulon_dorothea}} function.
#' @param signaling_data numerical vector, where names are signaling nodes 
#' in the PKN and values are from \{1, 0, -1\}. Continuous data will be 
#' discretized using the \code{\link{sign}} function.  
#' @param metabolic_data numerical vector, where names are metabolic nodes 
#' in the PKN and values are continuous values that represents log2 fold change 
#' or t-values from a differential analysis. These values are compared to 
#' the simulation results (simulated nodes can take value -1, 0 or 1)
#' @param diff_expression_data (optional) numerical vector that represents the 
#' results of a differential gene expression analysis. Names are gene
#' names using EntrezID starting with an X and values are log fold change or
#'  t-values.  \code{\link{convert_genesymbols_to_entrezid}} can be used for
#'  conversion.  We use the \dQuote{\code{diff_exp_threshold} parameter to decide
#'  which genes changed significantly.  Genes with NA values are considered none
#'  expressed and they will be removed from the TF-gene expression interactions. 
#' @param diff_exp_threshold threshold parameter (default 1) used to binarize
#'  the values of \dQuote{\code{diff_expression_data}. 
#' @param maximum_network_depth integer > 0 (default: 8). Nodes that are further 
#' than \dQuote{\code{maximum_network_depth} steps from the signaling nodes on 
#' the directed graph of the PKN are considered non-reachable and are removed. 
#' @param remove_unexpressed_nodes if TRUE (default) removes nodes from the PKN 
#' that are not expressed, see input \dQuote{\code{expressed_genes}.
#' @param expressed_genes character vector. Names of nodes that are expressed. By 
#' default we consider all the nodes that appear in \code{diff_expression_data} with
#' a numeric value (i.e. nodes with NA are removed) 
#' @param filter_tf_gene_interaction_by_optimization (default:TRUE), if TRUE then runs 
#' a network optimization that estimates TF activity not included in the inputs
#' and checks the consistency between the estimated activity and change in gene 
#' expression. Removes interactions where TF and gene expression are inconsistent 
#' @param CARNIVAL_options list that controls the options of CARNIVAL. See details 
#'  in \code{\link{default_CARNIVAL_options}}. 
#' @export
#' @return cosmos_data object with the following fields:
#'   \describe{
#'     \item{\code{meta_network}}{Filtered PKN}
#'     \item{\code{tf_regulon}}{TF - target regulatory network}
#'     \item{\code{signaling_data_bin}}{Binarised signaling data} 
#'     \item{\code{metabolic_data}}{Metabolomics data}
#'     \item{\code{diff_expression_data_bin}}{Binarized gene expression data} 
#'     \item{\code{optimized_network}}{Initial optimized network if 
#'     \code{filter_tf_gene_interaction_by_optimization is TRUE}}
#'   }
#' @seealso \code{\link{meta_network}} for meta PKN,
#' \code{\link{load_tf_regulon_dorothea}} for tf regulon,
#' \code{\link{convert_genesymbols_to_entrezid}} for gene conversion,
#' \code{\link[CARNIVAL]{runCARNIVAL}}.
#' 
#' @examples
#' CARNIVAL_options <- cosmos::default_CARNIVAL_options()
#' CARNIVAL_options$solver <- "lpSolve"
#' test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = toy_network,
#' signaling_data = toy_signaling_input,
#' metabolic_data = toy_metabolic_input,
#' diff_expression_data = toy_RNA,
#' maximum_network_depth = 15,
#' remove_unexpressed_nodes = T,
#' CARNIVAL_options = CARNIVAL_options
#' )

preprocess_COSMOS_signaling_to_metabolism <- function(meta_network = meta_network,
                              tf_regulon = load_tf_regulon_dorothea(),
                              signaling_data,
                              metabolic_data,
                              diff_expression_data, 
                              diff_exp_threshold = 1,
                              maximum_network_depth = 8,
                              expressed_genes =  names(diff_expression_data)[!is.na(diff_expression_data)],
                              remove_unexpressed_nodes = TRUE,
                              filter_tf_gene_interaction_by_optimization = TRUE,
                              CARNIVAL_options = default_CARNIVAL_options()){
    
    out_data <- preprocess_COSMOS_core(meta_network = meta_network,
                                       tf_regulon = tf_regulon,
                                       signaling_data=signaling_data,
                                       metabolic_data=metabolic_data,
                                       input_layer = "signaling_data",
                                       output_layer = "metabolic_data",
                                       diff_expression_data=diff_expression_data, 
                                       diff_exp_threshold = diff_exp_threshold,
                                       maximum_network_depth = maximum_network_depth,
                                       expressed_genes =  expressed_genes,
                                       remove_unexpressed_nodes = remove_unexpressed_nodes,
                                       filter_tf_gene_interaction_by_optimization = filter_tf_gene_interaction_by_optimization,
                                       CARNIVAL_options = CARNIVAL_options)
    
    out_data$history  = c(out_data$history,"created by preprocess_COSMOS_signaling_to_metaboism\n")
    
    return(out_data)
}
