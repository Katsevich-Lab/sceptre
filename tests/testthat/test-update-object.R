# An object of 0.99.1 and earlier: the matrices in one-element lists and the
# covariates in their own slot, in place of the `data` slot.
old_layout_object <- function(sceptre_object) {
    old <- sceptre_object
    attr(old, "response_matrix") <- list(get_response_matrix(sceptre_object))
    attr(old, "grna_matrix") <- list(get_grna_matrix(sceptre_object))
    attr(old, "covariate_data_frame") <- get_cell_covariates(sceptre_object)
    attr(old, "data") <- NULL
    attr(old, "seed") <- NULL
    attr(old, "rng_protocol") <- NULL
    old
}

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
    expect_identical(updated@rng_protocol, "mt19937-1")
    expect_identical(updateObject(sceptre_object), sceptre_object)
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
