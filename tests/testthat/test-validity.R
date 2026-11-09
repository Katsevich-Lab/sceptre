test_that("sceptre_objects at each step of the pipeline are valid", {
    data(highmoi_example_data)
    data(grna_target_data_frame_highmoi)
    sceptre_object <- import_data(
        response_matrix = highmoi_example_data$response_matrix,
        grna_matrix = highmoi_example_data$grna_matrix,
        grna_target_data_frame = grna_target_data_frame_highmoi,
        moi = "high",
        extra_covariates = highmoi_example_data$extra_covariates,
        response_names = highmoi_example_data$gene_names
    )
    expect_true(methods::validObject(sceptre_object))
    sceptre_object <- sceptre_object |>
        set_analysis_parameters(
            positive_control_pairs = construct_positive_control_pairs(
                sceptre_object
            )
        ) |>
        assign_grnas(method = "thresholding") |>
        run_qc()
    expect_true(methods::validObject(sceptre_object))
    expect_true(methods::validObject(methods::new("sceptre_object")))
})

test_that("validity reports inconsistent slots", {
    data(lowmoi_example_data)
    sceptre_object <- import_data(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
        moi = "low"
    )
    invalid <- sceptre_object
    invalid@cells_in_use <- ncol(get_response_matrix(sceptre_object)) + 1L
    expect_error(methods::validObject(invalid), "`cells_in_use`")
    invalid <- sceptre_object
    invalid@low_moi <- c(TRUE, FALSE)
    expect_error(methods::validObject(invalid), "`low_moi`")
    invalid <- sceptre_object
    invalid@side_code <- 2L
    expect_error(methods::validObject(invalid), "`side_code`")
    invalid <- sceptre_object
    invalid@multiple_testing_alpha <- 1.5
    expect_error(methods::validObject(invalid), "`multiple_testing_alpha`")
    invalid <- sceptre_object
    invalid@grna_target_data_frame <- data.frame(grna = "g")
    expect_error(methods::validObject(invalid), "`grna_target_data_frame`")
    invalid <- sceptre_object
    SummarizedExperiment::assayNames(invalid@data) <- "umis"
    expect_error(methods::validObject(invalid), "assay \"counts\"")
})
