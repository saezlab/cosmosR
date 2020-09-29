
test_that("test data", {
    
    # should not work with dataframe
    test_data = tibble(ID = c("A","B"),activity = c(0.1,5))
    expect_error(check_COSMOS_inputs(signaling_data = test_data))
    
    # should pass with named vector
    test_data = c(X0125=4,X123684=-4)
    expect_true(check_COSMOS_inputs(signaling_data = test_data))
    
})

test_that("test meta network", {
    
    test_data = tibble(source = c("A","B"),interaction = c(-1,1), target = c("C","D"))
    expect_true(check_COSMOS_inputs(meta_network = test_data))
    
})


test_that("test regulon network", {
    
    test_data = tibble(tf = c("A","B"),sign = c(-1,1), target = c("C","D"))
    expect_true(check_COSMOS_inputs(tf_regulon = test_data))
    
})
