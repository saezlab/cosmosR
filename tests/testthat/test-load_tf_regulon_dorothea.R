


test_that("the retired DoRothEA loader gives an actionable migration error", {
    expect_error(
        load_tf_regulon_dorothea(confidence = c("A", "B", "C")),
        "retired"
    )
})
