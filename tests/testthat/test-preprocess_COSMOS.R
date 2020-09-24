
test_that("test run COSMOS", {
    
    
    meta_network_test <- cosmos:::meta_network_test
    signaling_input_test <- cosmos:::signaling_input_test
    expression_data_test <- cosmos:::expression_data_test
    metabolic_data_test <- cosmos:::metabolic_data_test
    
    res <- preprocess_COSMOS(meta_network = meta_network_test,
                             signaling_data = signaling_input_test,
                             metabolic_data =  metabolic_data_test,
                             expression_data = expression_data_test,
                             filter_tf_gene_interaction_by_optimization = FALSE)
    
    expect_length(res, 6)
    expect_true(all(c("meta_network",
                      "tf_regulon",
                      "signaling_data_bin",
                      "metabolic_data",
                      "expression_data", 
                      "optimized_network") %in% names(res)))
})
