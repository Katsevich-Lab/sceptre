test_that("conversion sums repeated triplets and reads stored FALSE as 0", {
    repeated <- Matrix::sparseMatrix(
        i = c(1, 1, 2),
        j = c(1, 1, 2),
        x = c(1, 2, 3),
        dims = c(2, 2),
        repr = "T"
    )
    stored_false <- Matrix::sparseMatrix(
        i = c(1, 2),
        j = c(1, 2),
        x = c(TRUE, FALSE),
        dims = c(2, 2)
    )
    for (row_access in c(TRUE, FALSE)) {
        actual <- set_matrix_accessibility(repeated, row_access)
        expect_true(methods::validObject(actual))
        expect_identical(as.matrix(actual), matrix(c(3, 0, 0, 3), 2, 2))
        actual <- set_matrix_accessibility(stored_false, row_access)
        expect_identical(as.matrix(actual), matrix(c(1, 0, 0, 0), 2, 2))
    }
})

# `convert_covariate_df_to_design_matrix` is used in `set_analysis_parameters`
# to take `sceptre_object@covariate_data_frame` and `sceptre_object@formula_object`
# and make the corresponding model matrix.
test_that("convert_covariate_df_to_design_matrix", {
    set.seed(187)
    n <- 12
    covariate_data_frame <- data.frame(
        x = rnorm(n),
        y = 1:n,
        z = factor(rep(0:1, each = n / 2), levels = 0:1)
    )
    fmla <- formula("~ x*z + log(y + 1)")

    ## testing errors with Inf or NA values
    FAIL_bad_values_covariate_data_frame <- covariate_data_frame
    FAIL_bad_values_covariate_data_frame[1, 1] <- Inf

    expect_error(
        convert_covariate_df_to_design_matrix(
            FAIL_bad_values_covariate_data_frame,
            formula_object = fmla
        ),
        regex = "contains entries that are -Inf, Inf, or NA"
    )

    FAIL_bad_values_covariate_data_frame[1, 1] <- NA
    expect_error(
        convert_covariate_df_to_design_matrix(
            FAIL_bad_values_covariate_data_frame,
            formula_object = fmla
        ),
        regex = "Some rows of `covariate_data_frame` were lost"
    )

    ## now testing if it's the formula that causes these
    FAIL_fmla <- formula("~ log(y - 1)") # should lead to an  -Inf from log(0)
    expect_error(
        convert_covariate_df_to_design_matrix(
            covariate_data_frame,
            formula_object = FAIL_fmla
        ),
        regex = "has been applied contains entries that are -Inf, Inf, or NA"
    )

    ## testing low rank
    FAIL_low_Rank_covariate_data_frame <- covariate_data_frame |>
        dplyr::mutate(low_rank_col = 2 * x - log(y + 1))
    # the formula doesn't use the low rank column so there should be no error currently
    expect_no_error(
        convert_covariate_df_to_design_matrix(
            FAIL_low_Rank_covariate_data_frame,
            formula_object = fmla
        )
    )
    expect_error(
        convert_covariate_df_to_design_matrix(
            FAIL_low_Rank_covariate_data_frame,
            formula_object = formula("~ x*z + log(y + 1) + low_rank_col")
        ),
        regex = "contains redundant information"
    )

    ## testing correctness
    results <- convert_covariate_df_to_design_matrix(
        covariate_data_frame,
        formula_object = fmla
    )
    expect_equal(
        ncol(results),
        5 # with intercept
    )
    expect_true(
        all(results[, 1] == 1)
    )
    expect_equal(
        results[, "x"] |> as.numeric(),
        covariate_data_frame$x
    )
    expect_equal(
        results[, "z1"] |> as.numeric(),
        covariate_data_frame$z |> as.character() |> as.numeric()
    )
    expect_equal(
        results[, "log(y + 1)"] |> as.numeric(),
        covariate_data_frame$y |> log1p()
    )
    expect_equal(
        results[, "x:z1"] |> as.numeric(),
        (covariate_data_frame$z |> as.character() |> as.numeric()) *
            covariate_data_frame$x
    )
})
