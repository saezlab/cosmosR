#' run_COSMOS_signaling2metabolism
#' 
#' runs COSMOS from signaling to metabolism
#' 
#' @param meta_network prior knowledge network. By default COSMOS use a PKN 
#' derived from Omnipath, STITCHdb and Recon3D. See details on the function 
#' \cite{\code{load_meta_pkn}}
#' @param signaling_data numerical vector, where names are signaling nodes 
#' in the PKN and values are from \{1, 0, -1\}. Continuous data will be 
#' discretized using the \link\code{sign} function.  
#' @param metabolic_data numerical vector, where names are metabolic nodes 
#' in the PKN and values are continuous values. These values are compared to 
#' with the the simulation [? range of values]
#' @param expression_data (optional) numerical vector, where names are gene
#' names and values are from \{1,-1\}
#' @param solver_path argument passed to \link\code{CARNIVAL::runCARNIVAL}
#' @param solver argument passed to \link\code{CARNIVAL::runCARNIVAL}
#' @param time_limit argument passed to \link\code{CARNIVAL::runCARNIVAL}
#' @param test_run (default:FALSE), if TRUE then runs the pipeline without the
#' last optimization step. Useful for diagnostics. 
#' @export
#' @import dorothea biomaRt igraph dplyr
run_COSMOS_signaling2metabolism <- function(meta_network = load_meta_pkn(),
                           tf_regulon = load_tf_regulon_dorothea(),
                           signaling_data,
                           metabolic_data,
                           expression_data,
                           solver_path = NULL, 
                           solver = "cplex",
                           time_limit = 3600,
                           test_run = FALSE){
    
    ## Checking COSMOS input format
    check_COSMOS_inputs(meta_network,
                        tf_regulon,
                        signaling_data,
                        metabolic_data,
                        expression_data)
    
    # Check overlap among node names in the inputs
    check_network_data_coverage(meta_network,
                                tf_regulon,
                                signaling_data,
                                metabolic_data,
                                expression_data)
    
    # Preprocess PKN
    # - cut unreachable nodes from inputs
    meta_network <- keep_downstream_neighbours(
        network = meta_network,
        n_steps =  8,
        starting_nodes = names(signaling_data))
    
    
    # Filter TF -> target interaction from PKN if target expression not changing
    filtered_meta_network <- filter_transcriptional_regulations(
        network = meta_network[,1:3], 
        gene_expression_binarized = binarize_with_sign(expression_data,
                                                       threshold = 1),
        signaling_data  = signaling_data,
        tf_regulon=tf_regulon[,c("tf","sign","target")])
    
    # run CARNIVAL from signaling to metabolism,
    # this may estimate the activity of other TF-s.
    if(!test_run){
    CARNIVAL_results = runCARNIVAL_wrapper(network = meta_network,
                        input_data = sign(signaling_data),
                        measured_data = measured_data,
                        solver_path = solver_path,
                        solver = solver,
                        timelimit = timelimit,
                        mipGAP = 0.2)
    
    }else(
        return(TRUE)
    )
    
    
    
}









# source("~/Dropbox/COSMOS/scripts/revision_COSMOS_functions.R")
# 
# 
# my_stat <- "t" 
# my_threshold <- 1
# 
# ##Dorothea/viper
# dorothea <- format_dorothea()
# 
# #Load inputs
# 
# meta_network <- read_csv("./data/meta_network_carnival_ready_exch_solved_fullomni_metfiltered.csv")
# 
# signaling_input_carnival <- as.data.frame(read_csv("data/for_attila/signaling_input_carnival_raw_new_doro.csv"))
# signaling_input <- unlist(signaling_input_carnival)
# 
# 
# 
# 
# metab_input_carnival <- as.data.frame(read_csv("data/for_attila/metab_input_carnival_raw.csv"))
# metabolic_data <- unlist(metab_input_carnival)
# 
# 
# ttop_RNA <- as.data.frame(read_csv("./data/for_attila/ttop_tumorvshealthy.csv"))
# expression_data = ttop_RNA$t
# names(expression_data) = paste0("X",ttop_RNA$ID)
# 
# input_names = names(signaling_input)
# 
# #Preprocess pkn
# 
# signaling_input_carnival <- filter_inputs(signaling_input_carnival, meta_network)
# 
# metab_input_carnival <- filter_inputs(metab_input_carnival, meta_network)
# 
# meta_network <- downstream_neighbours(meta_network, 8, names(signaling_input_carnival))
# ###############
# 
# 
# 
# old_pkn <- filter_TF_sign(meta_network = meta_network[,1:3],
#                           ttop_RNA = ttop_RNA,
#                           inputs = signaling_input_carnival, 
#                           TF_targets = dorothea,my_stat = "t",my_threshold = 1)
# 
# 
# new_pkn <- filter_transcriptional_regulations(network = meta_network[,1:3], 
#                                               gene_expression_binarized = binarize_with_sign(expression_data,threshold = 1),
#                                               signaling_data  =signaling_data,
#                                               tf_regulon=tf_regulon[,c("tf","sign","target")])
# 
# 
# 
# ###############
# 
# CARNIVAL_Result <- runCARNIVAL(inputObj = sign(signaling_input_carnival),
#                                measObj = metab_input_carnival,
#                                netObj = meta_network,
#                                solverPath = "~/Documents/cplex",
#                                solver = "cplex",
#                                timelimit = 3600,
#                                mipGAP = 0.2)
# 
# save.image("carni_forwardrun_res_fullomni_fullfilter.RData")
# 
# TF_signs <- get_TF_sign_from_CARNI(CARNIVAL_Result, dorothea)
# 
# meta_network <- filter_TF_sign(meta_network = meta_network, ttop_RNA = ttop_RNA, inputs = TF_signs, TF_targets = dorothea)
# 
# CARNIVAL_Result_rerun <- runCARNIVAL(inputObj = sign(signaling_input_carnival),
#                                      measObj = metab_input_carnival,
#                                      netObj = meta_network,
#                                      solverPath = "~/Documents/cplex",
#                                      solver = "cplex",
#                                      timelimit = 7200,
#                                      mipGAP = 0.2)
# 
# save.image("carni_forwardrun_res_fullomni_fullfilter.RData")
# 
# meta_network <- as.data.frame(
#     read_csv("~/Dropbox/COSMOS/support/meta_network_carnival_ready_exch_solved_fullomni_metfiltered_expfiltered.csv"))
# 
# ### FOR RELIEF COSMOS CHEM TO SIG
# 
# # input_network <- meta_network[meta_network$source %in% names(metab_input_carnival),]
# # 
# # meta_network$edgeID <- paste0(meta_network$source, meta_network$interaction, meta_network$target)
# # input_network$edgeID <- paste0(input_network$source, input_network$interaction, input_network$target)
# # 
# # meta_network <- meta_network[!(meta_network$edgeID %in% input_network$edgeID),]
# # 
# # meta_network <- meta_network[,-4]
# # input_network <- input_network[,-4]
# # 
# # input_network_chem_to_sig <- input_network
# # input_network_chem_to_sig$source <- gsub("___[a-z]____", "",input_network_chem_to_sig$source)
# # 
# # ##For later
# # input_network_conector <- input_network
# # input_network_conector$target <- input_network_chem_to_sig$source
# # input_network_conector$interaction <- 1
# # input_network_conector <- unique(input_network_conector)
# # ##
# # 
# # input_network_chem_to_comp <- input_network
# # input_network_chem_to_comp$target <- input_network_chem_to_comp$source
# # input_network_chem_to_comp$source <- input_network_chem_to_sig$source
# # input_network_chem_to_comp$interaction <- 1
# # input_network_chem_to_comp <- unique(input_network_chem_to_comp)
# # 
# # input_network_chem_to_sig <- unique(as.data.frame(rbind(input_network_chem_to_sig,input_network_chem_to_comp)))
# # 
# # meta_network <- as.data.frame(rbind(input_network_chem_to_sig, meta_network))
# # 
# # metab_input_carnival_reduced <- as.data.frame(t(metab_input_carnival))
# # metab_input_carnival_reduced$metab <- gsub("___[a-z]____", "",row.names(metab_input_carnival_reduced))
# # metab_input_carnival_reduced <- unique(metab_input_carnival_reduced)
# # row.names(metab_input_carnival_reduced) <- metab_input_carnival_reduced$metab
# # metab_input_carnival_reduced <- metab_input_carnival_reduced[,-2,drop = F]
# # 
# # metab_input_carnival_reduced <- as.data.frame(t(metab_input_carnival_reduced))
# # 
# # meta_network <- downstream_neighbours(meta_network, 6, names(metab_input_carnival_reduced))
# 
# ### FOR RELIEF COSMOS CHEM TO SIG (END)
# 
# meta_network <- downstream_neighbours(meta_network, 6, names(metab_input_carnival))
# 
# ######
# ######
# ######
# meta_network <- filter_TF_sign(meta_network = meta_network, ttop_RNA = ttop_RNA, inputs = signaling_input_carnival, TF_targets = dorothea)
# 
# meta_network <- filter_TF_sign(meta_network = meta_network, ttop_RNA = ttop_RNA, inputs = TF_signs, TF_targets = dorothea)
# 
# ####
# ####
# ####
# 
# # CARNIVAL_Result_2 <- runCARNIVAL(inputObj = sign(metab_input_carnival_reduced), 
# #                                  measObj = signaling_input_carnival,
# #                                  netObj = meta_network,
# #                                  solverPath = "~/Documents/cplex", 
# #                                  solver = "cplex", 
# #                                  timelimit = 7200,
# #                                  mipGAP = 0.2)
# 
# CARNIVAL_Result_2 <- runCARNIVAL(inputObj = sign(metab_input_carnival),
#                                  measObj = signaling_input_carnival,
#                                  netObj = meta_network,
#                                  solverPath = "~/Documents/cplex",
#                                  solver = "cplex",
#                                  timelimit = 14400,
#                                  mipGAP = 0.2)
# 
# save.image("carni_doublerun_res_fullomni_met_expfiltered.RData")
# 
# ttop_tumorvshealthy <- as.data.frame(read_csv("~/Dropbox/COSMOS/results/transcriptomic/ttop_tumorvshealthy.csv"))
# 
# RNAseq_entrez_to_symbol <- as.data.frame(read_delim("~/Dropbox/COSMOS/support/RNAseq_entrez_to_symbol",
#                                                     "\t", escape_double = FALSE, col_types = cols(`yourlist:M20191127216DA2B77BFBD2E6699CA9B6D1C41EB259129CL` = col_character()),
#                                                     trim_ws = TRUE)) #from uniprot 20191127
# 
# names(RNAseq_entrez_to_symbol)[1] <- "ID"
# 
# ttop_tumorvshealthy <- merge(ttop_tumorvshealthy, RNAseq_entrez_to_symbol[,c(1,6)])
# ttop_tumorvshealthy <- ttop_tumorvshealthy[,c(8,2:7)]
# names(ttop_tumorvshealthy)[1] <- "ID"
# ttop_tumorvshealthy$ID <- gsub(" .*","",ttop_tumorvshealthy$ID)