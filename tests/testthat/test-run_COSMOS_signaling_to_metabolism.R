


test_that("test run COSMOS signaling to metabolism", {
  
  
  meta_network <- cosmosR:::meta_network_test
  signaling_data <- cosmosR:::signaling_input_test
  expression_data <- cosmosR:::expression_data_test
  metabolic_data <- cosmosR:::metabolic_data_test
    
  
  res <- preprocess_COSMOS_signaling_to_metabolism(signaling_data = signaling_data,
                           meta_network = meta_network,
                           metabolic_data = metabolic_data,
                           diff_expression_data = expression_data,
                           remove_unexpressed_nodes = FALSE,
                           maximum_network_depth = 15,
                           filter_tf_gene_interaction_by_optimization = FALSE)
  
  CARNIVAL_options = CARNIVAL::defaultLpSolveCarnivalOptions()
  
  
  #cplex_file <- "/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex"
  #skip_if(!file.exists(cplex_file),"CPLEX optimization based test skipped.")
  
  #CARNIVAL_options$solverPath = cplex_file
  res_network = run_COSMOS_signaling_to_metabolism(data = res,
                                     CARNIVAL_options = CARNIVAL_options)
  
                                        
  expect_length(res_network, 5)
  expect_true(all(c("weightedSIF",
                    "N_networks",
                    "nodesAttributes",
                    "individual_networks",
                    "individual_networks_node_attributes") %in% names(res_network)))
 
})
