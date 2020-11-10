


test_that("test run COSMOS signaling to metabolism", {
  
  meta_network <- cosmos:::meta_network_test
  signaling_data <- cosmos:::signaling_input_test
  expression_data <- cosmos:::expression_data_test
  metabolic_data <- cosmos:::metabolic_data_test
  
  
  res <- preprocess_COSMOS_signaling_to_metabolism(signaling_data = signaling_data,
                           meta_network = meta_network,
                           metabolic_data = metabolic_data,
                           diff_expression_data = expression_data,
                           remove_unexpressed_nodes = TRUE,
                           filter_tf_gene_interaction_by_optimization = FALSE)
  
  CARNIVAL_options = cosmos::default_CARNIVAL_options()
  CARNIVAL_options$solverPath = "/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex"
  
  res_network = run_COSMOS_signaling_to_metabolism(data = res,
                                     CARNIVAL_options = CARNIVAL_options,
                                     test_run=TRUE)
  
                                        
  expect_length(res_network, 5)
  expect_true(all(c("aggregated_network",
                    "N_networks",
                    "aggregated_network_node_attributes",
                    "individual_networks",
                    "individual_networks_node_attributes") %in% names(res_network)))
 
})
