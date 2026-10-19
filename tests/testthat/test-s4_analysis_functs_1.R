# most of the content of `import_data` is subfunctions
# which are tested elsewhere
test_that("import_data", {
    set.seed(12321)
    num_cells <- 29
    num_responses <- 17
    grna_target_data_frame <- make_mock_grna_target_data(c(1, 3, 1), 1, 1, 6)
    grna_matrix <- matrix(
        rpois(nrow(grna_target_data_frame) * num_cells, 1),
        ncol = num_cells
    ) |>
        `rownames<-`(grna_target_data_frame$grna_id) |>
        set_matrix_accessibility()
    response_matrix <- matrix(
        rpois(num_responses * num_cells, 1),
        ncol = num_cells
    ) |>
        as("TsparseMatrix")
    extra_covariates <- data.frame(x = rep("aaa", num_cells))

    sceptre_object_low_no_ec_with_response_names <- import_data(
        response_matrix = response_matrix,
        grna_matrix = grna_matrix,
        grna_target_data_frame = grna_target_data_frame,
        moi = "low",
        response_names = paste0("rrr", 1:num_responses)
    )
    sceptre_object_high_with_ec <- import_data(
        response_matrix = response_matrix,
        grna_matrix = grna_matrix,
        grna_target_data_frame = grna_target_data_frame,
        moi = "high",
        extra_covariates = extra_covariates
    )

    expect_identical(slotNames(sceptre_object_high_with_ec),
        c("data", "spec", "preparation", "fits", "results", "last_function_called",
          "functs_called", "nuclear", "nf_pipeline", "elements_to_analyze", "schema_version"))
    for (object in list(sceptre_object_low_no_ec_with_response_names,
                        sceptre_object_high_with_ec)) {
        expect_identical(object@schema_version, 2L)
        expect_s4_class(object@data, "ScreenData")
        expect_s4_class(object@spec, "AnalysisSpec")
        expect_s4_class(object@preparation, "PreparationRecord")
        expect_s4_class(object@fits, "FitRecord")
        expect_s4_class(object@results, "ResultRecord")
    }

    ##### slots set in section 4 of `import_data`
    expect_equal(
        get_response_matrix(sceptre_object_high_with_ec),
        set_matrix_accessibility(response_matrix)
    )
    expect_equal(
        get_grna_matrix(sceptre_object_high_with_ec),
        set_matrix_accessibility(grna_matrix)
    )

    expect_true(sceptre_object_low_no_ec_with_response_names@data@low_moi)
    expect_false(sceptre_object_high_with_ec@data@low_moi)

    expect_equal(
        sceptre_object_low_no_ec_with_response_names@data@covariate_data_frame |>
            rownames(),
        as.character(seq_len(num_cells))
    )

    expect_equal(
        names(
            sceptre_object_low_no_ec_with_response_names@data@covariate_data_frame
        ),
        c(
            "response_n_nonzero",
            "response_n_umis",
            "grna_n_nonzero",
            "grna_n_umis"
        )
    )
    expect_equal(
        names(sceptre_object_high_with_ec@data@covariate_data_frame),
        c(
            "response_n_nonzero",
            "response_n_umis",
            "grna_n_nonzero",
            "grna_n_umis",
            "x"
        )
    )
    expect_equal(
        sceptre_object_high_with_ec@data@covariate_data_frame$x,
        extra_covariates$x
    )

    ##### slots set in section 5 of `import_data`
    expect_equal(
        sceptre_object_high_with_ec@last_function_called,
        "import_data"
    )
    expect_equal(
        names(sceptre_object_high_with_ec@functs_called),
        c(
            "import_data",
            "set_analysis_parameters",
            "assign_grnas",
            "run_qc",
            "run_calibration_check",
            "run_power_check",
            "run_discovery_analysis"
        )
    )
    expect_equal(
        sceptre_object_high_with_ec@functs_called |> as.logical(),
        rep(c(TRUE, FALSE), c(1, 6))
    )
})


# Check the public parameter-setting behavior here. Input-validation, saved-fit
# reuse and guide-grouping cases also have dedicated tests in this suite.

test_that("set_analysis_parameters", {
    set.seed(12333)
    grna_target_data_frame <- make_mock_grna_target_data(
        num_guides_per_target = c(1, 4, 7),
        chr_distances = 1,
        chr_starts = 1,
        num_nt_guides = 5
    )
    num_cells <- 23
    num_responses <- 18
    num_targets <- nrow(grna_target_data_frame)
    response_matrix <- make_mock_response_matrices(
        num_responses,
        num_cells,
        patterns = "column"
    )
    grna_matrix <- make_mock_grna_matrices(
        grna_target_data_frame,
        num_cells,
        non_nt_patterns = "column",
        nt_patterns = "row"
    )

    disc_pairs_empty <- matrix("", 0, 2) |>
        as.data.frame() |>
        setNames(c("grna_target", "response_id"))

    sceptre_object_post_all_defaults_low_moi <- import_data(
        response_matrix = response_matrix,
        grna_matrix = grna_matrix,
        grna_target_data_frame = grna_target_data_frame,
        extra_covariates = data.frame(),
        moi = "low"
    ) |>
        set_analysis_parameters(
            discovery_pairs = disc_pairs_empty
        )

    sceptre_object_post_all_defaults_high_moi <- import_data(
        response_matrix = response_matrix,
        grna_matrix = grna_matrix,
        grna_target_data_frame = grna_target_data_frame,
        moi = "high"
    ) |>
        set_analysis_parameters(
            discovery_pairs = disc_pairs_empty,
        )

    ## 0. if we call `set_analysis_parameters` again on the results nothing should change
    expect_equal(
        sceptre_object_post_all_defaults_low_moi,
        sceptre_object_post_all_defaults_low_moi |>
            set_analysis_parameters(discovery_pairs = disc_pairs_empty)
    )

    ## 1. do default values get set correctly based on MOI?

    # complement should be TRUE for high MOI and FALSE for low
    expect_true(
        sceptre_object_post_all_defaults_high_moi@spec@control_group_complement
    )
    expect_false(
        sceptre_object_post_all_defaults_low_moi@spec@control_group_complement
    )

    # permutations should be used for low MOI and CRT for high
    expect_false(sceptre_object_post_all_defaults_high_moi@spec@run_permutations)
    expect_true(sceptre_object_post_all_defaults_low_moi@spec@run_permutations)

    ## 2. does the formula behave correctly?
    expect_equal(
        as.character(sceptre_object_post_all_defaults_high_moi@spec@formula_object)[
            2
        ],
        "log(response_n_nonzero + 1) + log(response_n_umis + 1) + log(grna_n_nonzero) + log(grna_n_umis)"
    )
    sceptre_object_nonsense_formula <- import_data(
        response_matrix = response_matrix,
        grna_matrix = grna_matrix,
        grna_target_data_frame = grna_target_data_frame,
        extra_covariates = data.frame(aaa = rnorm(num_cells)), # needed for the formula
        moi = "low"
    ) |>
        set_analysis_parameters(
            discovery_pairs = disc_pairs_empty,
            formula_object = formula("~ aaa")
        )
    expect_equal(
        as.character(sceptre_object_nonsense_formula@spec@formula_object)[2],
        "aaa"
    )
})
