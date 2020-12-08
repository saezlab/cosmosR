library(cosmos)
# library(stringr)

my_options <- default_CARNIVAL_options()
my_options$solverPath <- "~/Documents/cplex"

#### FORWARD

test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = toy_sif,
                                                      signaling_data = toy_signaling_input_carnival_vec,
                                                      metabolic_data = toy_metab_input_carnival_vec,
                                                      diff_expression_data = toy_RNA,
                                                      maximum_network_depth = 15,
                                                      remove_unexpressed_nodes = T,
                                                      CARNIVAL_options = my_options
                                                      
)

test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                      CARNIVAL_options = my_options)

names(test_result_for)[1] <- 'weightedSIF'
names(test_result_for$weightedSIF) <- c("Node1","Sign","Node2","Weight")
names(test_result_for)[3] <- 'nodesAttributes'

metab_to_pubchem_vec <- metab_to_pubchem$name
names(metab_to_pubchem_vec) <- metab_to_pubchem$pubchem

test_result_for <- format_COSMOS_res(test_result_for,
                                     metab_mapping = metab_to_pubchem_vec,
                                     measured_nodes = unique(c(names(toy_metab_input_carnival_vec),names(toy_signaling_input_carnival_vec))),
                                     omnipath_ptm = omnipath_ptm)

View(test_result_for[[1]]) #SIF
View(test_result_for[[2]]) #ATTRIBUTES
#### BACKWARD

test_back <- preprocess_COSMOS_metabolism_to_signaling(meta_network = toy_sif,
                                                       signaling_data = toy_signaling_input_carnival_vec,
                                                       metabolic_data = toy_metab_input_carnival_vec,
                                                       diff_expression_data = toy_RNA,
                                                       maximum_network_depth = 15,
                                                       remove_unexpressed_nodes = F,
                                                       CARNIVAL_options = my_options
                                                       
)

test_result_back <- run_COSMOS_metabolism_to_signaling(data = test_back,
                                                       CARNIVAL_options = my_options)

names(test_result_back)[1] <- 'weightedSIF'
names(test_result_back$weightedSIF) <- c("Node1","Sign","Node2","Weight")
names(test_result_back)[3] <- 'nodesAttributes'

test_result_back <- format_COSMOS_res(test_result_back,
                                      metab_mapping = metab_to_pubchem_vec,
                                      measured_nodes = unique(c(names(toy_metab_input_carnival_vec),names(toy_signaling_input_carnival_vec))),
                                      omnipath_ptm = omnipath_ptm)

View(test_result_back[[1]]) #SIF
View(test_result_back[[2]]) #ATTRIBUTES
