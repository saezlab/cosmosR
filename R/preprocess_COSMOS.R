#' preprocess_COSMOS
#' 
#' core function for preprocessing. 
#' runs checks on the input data and simplifies the prior knowledge network.
#' Simplification includes the removal of (1) nodes that are not reachable from 
#' signaling nodes and  (2) interactions between transcription factors (TFs) and
#' target genes if the target gene does not respond or the response is
#' contradictory with the change in the transcription factor activity. 
#' Optionally, further TF activities are estimated via network optimization via
#' CARNIVAL and the interactions between TF and genes are filtered again. 
#' 
#' @param meta_network prior knowledge network. By default COSMOS use a PKN 
#' derived from Omnipath, STITCHdb and Recon3D. See details on the function 
#' \cite{\code{load_meta_pkn}}.
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
#'  t-values. We use the `diff_exp_threshold` parameter to decide which genes 
#'  changed significantly. Genes with NA values are considered none expressed
#'  and they will be removed from the TF-gene expression interactions. 
#' @param diff_exp_threshold threshold parameter (default 1) used to binarize
#'  the values of \dQuote{\code{diff_expression_data}}. 
#' @param maximum_network_depth integer > 0 (default: 8). Nodes that are further 
#' than `maximum_network_depth` steps from the signaling nodes on the directed
#' graph of the PKN are considered non-reachable and are removed. 
#' @param remove_unexpressed_nodes if TRUE (default) removes nodes from the PKN 
#' that are not expressed, see input `expressed_genes`.
#' @param expressed_genes character vector. Names of nodes that are expressed. By 
#' default we consider all the nodes that appear in \code{diff_expression_data} with
#' a numeric value (i.e. nodes with NA are removed) 
#' @param filter_tf_gene_interaction_by_optimization (default:TRUE), if TRUE then runs 
#' a network optimization that estimates TF activity not included in the inputs
#' and checks the consistency between the estimated activity and change in gene 
#' expression. Removes interactions where TF and gene expression are inconsistent 
#' @param CARNIVAL_options list that controls the options of CARNIVAL. See details 
#'  in \code{\link{default_CARNIVAL_options}}. 
#' @param input_layer either signaling_data or metabolic_data. Influences the way
#' network pruning and CARNIVAL is ran. 
#' @param output_layer either signaling_data or metabolic_data. Influences the way
#' network pruning and CARNIVAL is ran.
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
#' @noRd
preprocess_COSMOS_core <- function(meta_network,
                                   tf_regulon,
                                   signaling_data,
                                   metabolic_data,
                                   input_layer = c("signaling_data","metabolic_data"),
                                   output_layer = c("metabolic_data","signaling_data"),
                                   diff_expression_data, 
                                   diff_exp_threshold,
                                   maximum_network_depth,
                                   expressed_genes,
                                   remove_unexpressed_nodes,
                                   filter_tf_gene_interaction_by_optimization,
                                   CARNIVAL_options){
    
    ## Checking COSMOS input format
    check_COSMOS_inputs(meta_network,
                        tf_regulon,
                        signaling_data,
                        metabolic_data,
                        diff_expression_data)
    
    check_gene_names(signaling_data,diff_expression_data,expressed_genes)
    
    # Check overlap among node names in the inputs
    check_network_data_coverage(meta_network,
                                tf_regulon,
                                signaling_data,
                                metabolic_data,
                                diff_expression_data)
    
    # Preprocess PKN ----
    
    # filter nodes that do not appear in Gene-expression data. 
    # if the gene has a t-value --> it is expressed. 
    if(remove_unexpressed_nodes){
        meta_network <- filter_pkn_expressed_genes(expressed_genes_entrez = expressed_genes,
                                                   meta_pkn = meta_network)
        # After modifying the PKN, inputs/measured nodes might get lost, we remove
        signaling_data = filter_input_nodes_not_in_pkn(pkn = meta_network,
                                                       data = signaling_data)
        metabolic_data = filter_input_nodes_not_in_pkn(pkn = meta_network,
                                                       data = metabolic_data)
    }
    
    if(input_layer == "signaling_data" & output_layer == "metabolic_data"){
        # Cut non-reachable nodes from inputs
        meta_network <- keep_controllable_neighbours(
            network = meta_network,
            n_steps =  maximum_network_depth, 
            input_nodes = names(signaling_data))
        
        # After modifying the PKN, inputs/measured nodes might get lost, we remove
        signaling_data = filter_input_nodes_not_in_pkn(pkn = meta_network,
                                                       data = signaling_data)
        metabolic_data = filter_input_nodes_not_in_pkn(pkn = meta_network,
                                                       data = metabolic_data)
        # Cut non-observable nodes from inputs
        meta_network <- keep_observable_neighbours(
            network = meta_network,
            n_steps =  maximum_network_depth, 
            observed_nodes = names(metabolic_data))
        
        # After modifying the PKN, inputs/measured nodes might get lost, we remove
        signaling_data = filter_input_nodes_not_in_pkn(pkn = meta_network,
                                                       data = signaling_data)
        metabolic_data = filter_input_nodes_not_in_pkn(pkn = meta_network,
                                                       data = metabolic_data)
        
    }else if(input_layer == "metabolic_data" & output_layer == "signaling_data"){
        # Cut non-reachable nodes from inputs
        meta_network <- keep_controllable_neighbours(
            network = meta_network,
            n_steps =  maximum_network_depth, 
            input_nodes = names(metabolic_data))
        
        # After modifying the PKN, inputs/measured nodes might get lost, we remove
        signaling_data = filter_input_nodes_not_in_pkn(pkn = meta_network,
                                                       data = signaling_data)
        metabolic_data = filter_input_nodes_not_in_pkn(pkn = meta_network,
                                                       data = metabolic_data)
        # Cut non-observable nodes from inputs
        meta_network <- keep_observable_neighbours(
            network = meta_network,
            n_steps =  maximum_network_depth, 
            observed_nodes = names(signaling_data))
        
        # After modifying the PKN, inputs/measured nodes might get lost, we remove
        signaling_data = filter_input_nodes_not_in_pkn(pkn = meta_network,
                                                       data = signaling_data)
        metabolic_data = filter_input_nodes_not_in_pkn(pkn = meta_network,
                                                       data = metabolic_data)
    }else {
        stop("invalid input_layer and output_layer combination")
    }
    
    
    # Filter TF -> target interaction from PKN if target expression not changing
    diff_expression_data_bin <- binarize_with_sign(diff_expression_data,
                                                   threshold = diff_exp_threshold)
    
    
    filtered_meta_network <- filter_transcriptional_regulations(
        network = meta_network, 
        gene_expression_binarized = diff_expression_data_bin,
        signaling_data  = signaling_data,
        tf_regulon = tf_regulon)
    
    # run CARNIVAL from signaling to metabolism,
    # this may estimate the activity of other TF-s.
    if(filter_tf_gene_interaction_by_optimization){
        
        check_CARNIVAL_options(CARNIVAL_options)
        
        if(input_layer == "signaling_data" & output_layer == "metabolic_data"){
            
            CARNIVAL_results = runCARNIVAL_wrapper(network = meta_network,
                                                   input_data = sign(signaling_data),
                                                   measured_data = metabolic_data,
                                                   options = CARNIVAL_options)
            
            # get the estimated activity of TFs from CARNIVAL results
            estimated_TF_activity <- get_TF_activity_from_CARNIVAL(CARNIVAL_results,
                                                                   tf_regulon$tf)
            
            filtered_meta_network <- filter_transcriptional_regulations(
                network = filtered_meta_network, 
                gene_expression_binarized = diff_expression_data_bin,
                signaling_data  = estimated_TF_activity,
                tf_regulon=tf_regulon[,c("tf","sign","target")])
            
        }else if(input_layer == "metabolic_data" & output_layer == "signaling_data"){
            
            CARNIVAL_results = runCARNIVAL_wrapper(network = meta_network,
                                                   input_data = sign(metabolic_data),
                                                   measured_data = signaling_data,
                                                   options = CARNIVAL_options)
            
            # get the estimated activity of TFs from CARNIVAL results
            estimated_TF_activity <- get_TF_activity_from_CARNIVAL(CARNIVAL_results,
                                                                   tf_regulon$tf)
            
            filtered_meta_network <- filter_transcriptional_regulations(
                network = filtered_meta_network, 
                gene_expression_binarized = diff_expression_data_bin,
                signaling_data  = estimated_TF_activity,
                tf_regulon=tf_regulon[,c("tf","sign","target")])
            
        }
        
    }else{
        CARNIVAL_results = list()
    }
    
    out_data <-  cosmos_data(meta_network = filtered_meta_network,
                             tf_regulon = tf_regulon,
                             signaling_data = signaling_data,
                             metabolic_data = metabolic_data,
                             expression_data = diff_expression_data)
    
    out_data$signaling_data_bin = sign(signaling_data)
    out_data$metabolic_data_bin = sign(metabolic_data)
    out_data$diff_expression_data_bin = diff_expression_data_bin
    out_data$optimized_network = CARNIVAL_results
    
    return(out_data)
}
