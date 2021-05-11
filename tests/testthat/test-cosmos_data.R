


test_that("construct cosmos_data", {
    
    meta_network_test <- cosmosR:::meta_network_test
    signaling_input_test <- cosmosR:::signaling_input_test
    expression_data_test <- cosmosR:::expression_data_test
    metabolic_data_test <- cosmosR:::metabolic_data_test
    
    cos <- new_cosmos_data(meta_network = meta_network_test,
                           tf_regulon = load_tf_regulon_dorothea(),
                           signaling_data = signaling_input_test,
                           expression_data = expression_data_test,
                           metabolic_data = metabolic_data_test)
    expect_true(is(cos,"cosmos_data"))
})



test_that("validate cosmos_data", {
    
    cos <- new_cosmos_data(meta_network = meta_network_test,
                           tf_regulon = load_tf_regulon_dorothea(),
                           signaling_data = signaling_input_test,
                           expression_data = expression_data_test,
                           metabolic_data = metabolic_data_test)
    expect_true({
        validate_cosmos_data(cos)
        TRUE
    })
})




