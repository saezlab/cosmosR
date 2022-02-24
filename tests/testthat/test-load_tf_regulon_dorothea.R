


test_that("tf conversion", {
    
    reg = load_tf_regulon_dorothea(confidence = c("A","B","C"))
    # make sure A,B,C still contains a lot of interactions
    expect_true(nrow(reg) > 1)
    # should return 3 columns, correct names 
    expect_true(ncol(reg) == 3)
    expect_true(all(c("tf","sign","target") %in% colnames(reg)))
    # should not contain missing values
    expect_true(all(complete.cases(reg)))
    
})
