


test_that("construct cosmos_data", {
    
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




