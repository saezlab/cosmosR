




test_that("test run COSMOS", {
  
    
    meta_network_test <- cosmos:::meta_network_test
    signaling_input_test <- cosmos:::signaling_input_test
    expression_data_test <- cosmos:::expression_data_test
    metabolic_data_test <- cosmos:::metabolic_data_test
    
    
                                        
    expect_true(run_COSMOS_signaling2metabolism(meta_network = meta_network_test,
                                                signaling_data = signaling_input_test,
                                                metabolic_data =  metabolic_data_test,
                                                expression_data = expression_data_test,
                                                test_run = TRUE))
    
})
