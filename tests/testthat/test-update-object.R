test_that("updateObject() converts an object of the 0.99.1 layout", {
    data(lowmoi_example_data)
    sceptre_object <- import_data(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
        extra_covariates = lowmoi_example_data$extra_covariates,
        moi = "low"
    ) |>
        set_analysis_parameters()
    old <- old_layout_object(sceptre_object)
    expect_false(methods::.hasSlot(old, "data"))

    updated <- updateObject(old)
    expect_true(methods::validObject(updated))
    expect_identical(
        get_response_matrix(updated),
        get_response_matrix(sceptre_object)
    )
    expect_identical(get_grna_matrix(updated), get_grna_matrix(sceptre_object))
    expect_identical(
        get_cell_covariates(updated),
        get_cell_covariates(sceptre_object)
    )
    expect_identical(updated@formula_object, sceptre_object@formula_object)
    expect_identical(updated@functs_called, sceptre_object@functs_called)
    expect_identical(updated@seed, 4L)
    expect_identical(updated@version, package_version("0.99.1"))
    expect_identical(updateObject(sceptre_object), sceptre_object)

    # an object with `data` but no `version` (a pre-release 0.99.2) is stamped
    # 0.99.1 too
    pre_release <- sceptre_object
    attr(pre_release, "version") <- NULL
    expect_identical(
        updateObject(pre_release)@version,
        package_version("0.99.1")
    )
})

test_that("a sceptre_object records the version of sceptre that created it", {
    data(lowmoi_example_data)
    sceptre_object <- import_data(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
        moi = "low"
    )
    expect_identical(sceptre_object@version, utils::packageVersion("sceptre"))
    expect_true(methods::validObject(sceptre_object))

    # an object created by a later version of sceptre is refused
    future <- sceptre_object
    future@version <- package_version("99.0.0")
    expect_error(get_response_matrix(future), "upgrade")
    expect_error(set_analysis_parameters(future), "upgrade")
})

test_that("an object of the 0.99.1 layout is refused until it is updated", {
    data(lowmoi_example_data)
    sceptre_object <- import_data(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
        moi = "low"
    )
    old <- old_layout_object(sceptre_object)
    expect_error(get_response_matrix(old), "updateObject()", fixed = TRUE)
    expect_error(get_cell_covariates(old), "updateObject()", fixed = TRUE)
    expect_error(set_analysis_parameters(old), "updateObject()", fixed = TRUE)
})

test_that("updateObject() asks to re-import an object without its matrices", {
    data(lowmoi_example_data)
    sceptre_object <- import_data(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
        moi = "low"
    )
    old <- old_layout_object(sceptre_object)
    attr(old, "response_matrix") <- list(list())
    expect_error(updateObject(old), "Re-import the data")
})

test_that("an updated object whose covariates lack cell names can be written", {
    # early versions could keep the cell names on the matrices only
    data(lowmoi_example_data)
    sceptre_object <- import_data(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
        moi = "low"
    )
    expect_false(is.null(colnames(get_response_matrix(sceptre_object))))
    old <- old_layout_object(sceptre_object)
    covariates <- attr(old, "covariate_data_frame")
    rownames(covariates) <- NULL
    attr(old, "covariate_data_frame") <- covariates
    updated <- updateObject(old)
    directory <- withr::local_tempdir()
    write_sceptre_object(updated, directory)
    reread <- read_sceptre_object(directory)
    expect_identical(response_values(reread), response_values(updated))
})
