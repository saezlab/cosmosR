


test_that("test run COSMOS signaling to metabolism", {
  
  meta_network <- cosmos:::meta_network_test
  signaling_data <- cosmos:::signaling_input_test
  expression_data <- cosmos:::expression_data_test
  metabolic_data <- cosmos:::metabolic_data_test
  
  
  res <- preprocess_COSMOS(signaling_data = signaling_data,
                           meta_network = meta_network,
                           metabolic_data = metabolic_data,
                           diff_expression_data = expression_data,
                           filter_tf_gene_interaction_by_optimization = FALSE)
  
  signaling_data <- res$signaling_data_bin
  meta_network <- res$meta_network
  diff_expression_data <- res$diff_expression_data_bin
  metabolic_data <- res$metabolic_data
  
  
  res_network = run_COSMOS_metabolism_to_signaling(meta_network = meta_network,
                                     metabolic_data = metabolic_data,
                                     signaling_data = signaling_data,
                                     test_run=TRUE)
  
                                        
  expect_length(res_network, 5)
  expect_true(all(c("aggregated_network",
                    "N_networks",
                    "aggregated_network_node_attributes",
                    "individual_networks",
                    "individual_networks_node_attributes") %in% names(res_network)))
 
})
