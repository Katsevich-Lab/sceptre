test_that("example data use consistent cell IDs", {
    data("lowmoi_example_data")
    data("highmoi_example_data")

    for (example_data in list(lowmoi_example_data, highmoi_example_data)) {
        cell_ids <- colnames(example_data$response_matrix)
        expect_identical(colnames(example_data$grna_matrix), cell_ids)
        expect_identical(rownames(example_data$extra_covariates), cell_ids)
    }
})
