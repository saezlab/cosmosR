test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})



test_that("test data", {
    
    test_data = tibble(ID = c("A","B"),activity = c(0.1,5))
    
    expect_error(check_COSMOS_inputs(signaling_data = test_data))
    
    
    test_data = c(A=4,B=4)
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
