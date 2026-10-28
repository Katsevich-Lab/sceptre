test_that("empirical p-values cap only the doubled two-sided tail at one", {
    draws <- fisher_yates_sampler(n_tot = 5L, M = 1L, B = 4L)
    scores <- list(
        rep(0, 5),
        c(0, -2, -1, 1, 2),
        c(0, 0, 0, 0, 1),
        c(1, -1, 0, 0, 1),
        c(-2, -2, -2, -1, 0),
        c(2, 0, 0, 1, 2)
    )
    capped <- 0L
    for (a in scores) {
        for (side in -1:1) {
            actual <- run_low_level_test_full_v4(
                y = rep(2, 5),
                mu = rep(2, 5),
                a = a,
                w = rep(1, 5),
                D = matrix(0, 1L, 5L),
                trt_idxs = 1L,
                n_trt = 1L,
                use_all_cells = FALSE,
                synthetic_idxs = draws,
                B1 = 4L,
                B2 = 0L,
                B3 = 0L,
                fit_parametric_curve = FALSE,
                return_resampling_dist = TRUE,
                side_code = side
            )
            null <- actual$resampling_dist
            observed <- actual$z_orig
            expect_identical(observed, a[1L])
            # inclusive add-one tail probabilities, computed independently
            left <- (1 + sum(null <= observed)) / (1 + length(null))
            right <- (1 + sum(null >= observed)) / (1 + length(null))
            old <- if (side == -1L) {
                left
            } else if (side == 1L) {
                right
            } else {
                2 * min(left, right)
            }
            expect_identical(actual$p, min(1, old))
            if (old <= 1) {
                expect_identical(actual$p, old)
            }
            capped <- capped + (old > 1)
        }
    }
    expect_gt(capped, 0L)
})

test_that("a discovery analysis with all-tied scores returns a p-value of one", {
    inputs <- refactor_make_inputs()
    inputs$response_matrix[1L, ] <- 3
    object <- refactor_import(inputs, refactor_case_definitions()$low_nt_perm)
    object <- set_analysis_parameters(
        object,
        discovery_pairs = data.frame(response_id = "r1", grna_target = "A"),
        formula_object = ~1,
        side = "both",
        resampling_mechanism = "permutations",
        resampling_approximation = "no_approximation",
        control_group = "complement"
    )
    object <- assign_grnas(
        object,
        method = "thresholding",
        threshold = 40,
        print_progress = FALSE
    )
    object <- suppressMessages(run_qc(
        object,
        response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1),
        n_nonzero_trt_thresh = 7L,
        n_nonzero_cntrl_thresh = 7L
    ))
    object <- suppressMessages(run_discovery_analysis(
        object,
        print_progress = FALSE
    ))
    result <- get_result(object, "run_discovery_analysis")
    expect_identical(nrow(result), 1L)
    expect_true(result$pass_qc)
    expect_identical(result$p_value, 1)
    expect_false(result$significant)
})
