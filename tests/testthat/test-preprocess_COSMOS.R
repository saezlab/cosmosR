
test_that("test cosmos preprocessing (signaling to metabolism)", {
    
    meta_network_test <- cosmos:::meta_network_test
    signaling_input_test <- cosmos:::signaling_input_test
    expression_data_test <- cosmos:::expression_data_test
    metabolic_data_test <- cosmos:::metabolic_data_test
    
    res <- preprocess_COSMOS_signaling_to_metabolism(meta_network = meta_network_test,
                             signaling_data = signaling_input_test,
                             metabolic_data =  metabolic_data_test,
                             diff_expression_data = expression_data_test,
                             filter_tf_gene_interaction_by_optimization = FALSE)
    # check list:
    expect_length(res, 10)
    expect_true(all(c("meta_network",
                      "tf_regulon",
                      "signaling_data",
                      "signaling_data_bin",
                      "metabolic_data",
                      "metabolic_data_bin",
                      "expression_data",
                      "diff_expression_data_bin", 
                      "optimized_network") %in% names(res)))
    
    
    # checking network
    expect_equal(ncol(res$meta_network), 3)
    expect_equal(nrow(res$meta_network), 44005)
    expect_true(all(colnames(res$meta_network) %in% c("source","interaction","target")))
    
    # checking tf_regulon
    expect_true(all(c("tf","target","sign") %in% colnames(res$tf_regulon)))
    
    # checking signaling
    expect_true(is.vector(res$signaling_data_bin))
    expect_equal(length(res$signaling_data_bin),103)
    expect_equal(sum(grepl("^X",names(res$signaling_data_bin))),103)
    expect_true(all(res$signaling_data_bin %in% c(-1,0,1)))
    
    # checking metabolic
    expect_true(is.vector(res$metabolic_data))
    expect_equal(length(res$metabolic_data),35)
    expect_equal(sum(grepl("^XMetab",names(res$metabolic_data))),35)
    
    # check diff_expression_data_bin
    expect_true(is.vector(res$diff_expression_data_bin))
    expect_equal(length(res$diff_expression_data_bin),15919)
    expect_equal(sum(grepl("^X",names(res$diff_expression_data_bin))),15919)
    expect_true(all(res$diff_expression_data_bin %in% c(-1,0,1)))
})


test_that("test cosmos preprocessing (signaling to metabolism)", {
    
    
    meta_network_test <- cosmos:::meta_network_test
    signaling_input_test <- cosmos:::signaling_input_test
    expression_data_test <- cosmos:::expression_data_test
    metabolic_data_test <- cosmos:::metabolic_data_test
    
    res <- preprocess_COSMOS_metabolism_to_signaling(meta_network = meta_network_test,
                                                     signaling_data = signaling_input_test,
                                                     metabolic_data =  metabolic_data_test,
                                                     diff_expression_data = expression_data_test,
                                                     filter_tf_gene_interaction_by_optimization = FALSE)
    # check list:
    expect_length(res, 10)
    expect_true(all(c("meta_network",
                      "tf_regulon",
                      "signaling_data",
                      "signaling_data_bin",
                      "metabolic_data",
                      "metabolic_data_bin",
                      "expression_data",
                      "diff_expression_data_bin", 
                      "optimized_network") %in% names(res)))
    
    
    # checking network
    expect_equal(ncol(res$meta_network), 3)
    expect_equal(nrow(res$meta_network), 42103)
    expect_true(all(colnames(res$meta_network) %in% c("source","interaction","target")))
    
    # checking tf_regulon
    expect_true(all(c("tf","target","sign") %in% colnames(res$tf_regulon)))
    
    # checking signaling
    expect_true(is.vector(res$signaling_data_bin))
    expect_equal(length(res$signaling_data_bin),107)
    expect_equal(sum(grepl("^X",names(res$signaling_data_bin))),107)
    expect_true(all(res$signaling_data_bin %in% c(-1,0,1)))
    
    # checking metabolic
    expect_true(is.vector(res$metabolic_data))
    expect_equal(length(res$metabolic_data),29)
    expect_equal(sum(grepl("^XMetab",names(res$metabolic_data))),29)
    
    # check diff_expression_data_bin
    expect_true(is.vector(res$diff_expression_data_bin))
    expect_equal(length(res$diff_expression_data_bin),15919)
    expect_equal(sum(grepl("^X",names(res$diff_expression_data_bin))),15919)
    expect_true(all(res$diff_expression_data_bin %in% c(-1,0,1)))
})
