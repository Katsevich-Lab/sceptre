test_that("native empirical p values cap only the doubled two-sided tail", {
    bridge <- refactor_native_bridge(test_path("fixtures", "refactor", "native-bridge.cpp"))
    banks <- list(rep(0, 4), c(-2, -1, 1, 2), c(0, 0, 0, 1),
        c(-1, 0, 0, 1), c(-2, -2, -1, 0))
    set.seed(719)
    state <- .Random.seed
    for (draws in banks) for (observed in -2:2) for (side in -1:1) {
        case <- list(y = rep(2, 5), mu = rep(2, 5), a = c(observed, draws),
            w = rep(1, 5), D = matrix(0, 1, 5), trt_idxs = 1L, n_trt = 1L,
            use_all_cells = FALSE, draws = as.list(1:4), B1 = 4L, B2 = 0L, B3 = 0L,
            fit_parametric_curve = FALSE, return_resampling_dist = TRUE, side_code = side)
        # Independent inclusive counts; preserve the original arithmetic below 1.
        left <- (1 + sum(draws <= observed)) / (1 + length(draws))
        right <- (1 + sum(draws >= observed)) / (1 + length(draws))
        old <- if (side == -1L) left else if (side == 1L) right else 2 * min(left, right)
        actual <- refactor_native_run_case(case, bridge)
        expect_identical(actual$p, min(1, old))
        if (side != 0L || old <= 1) expect_identical(actual$p, old)
        expect_identical(actual$z_orig, as.double(observed))
        expect_identical(actual$resampling_dist, as.double(draws))
        expect_identical(actual$stage, 1L)
        expect_identical(.Random.seed, state)
    }
})

constant_score_screen <- function() {
    inputs <- refactor_make_inputs()
    inputs$response_matrix[1L, ] <- 3
    object <- refactor_import(inputs, refactor_case_definitions()$low_nt_perm)
    object <- refactor_capture(set_analysis_parameters(object,
        discovery_pairs = data.frame(response_id = "r1", grna_target = "A"),
        formula_object = ~ 1, side = "both", resampling_mechanism = "permutations",
        resampling_approximation = "no_approximation", control_group = "complement"))$value
    object <- refactor_capture(assign_grnas(object, method = "thresholding",
        threshold = 40, print_progress = FALSE))$value
    refactor_capture(run_qc(object, response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1), n_nonzero_trt_thresh = 7L,
        n_nonzero_cntrl_thresh = 7L))$value
}

test_that("public all-tied finite comparison returns p one and a complete result", {
    object <- constant_score_screen()
    set.seed(913)
    state <- .Random.seed
    direct <- refactor_capture(run_discovery_analysis(object,
        print_progress = FALSE, output_amount = 3L))$value
    result <- get_result(direct, "run_discovery_analysis")
    expect_identical(nrow(result), 1L)
    expect_true(result$pass_qc)
    expect_identical(result$p_value, 1)
    expect_identical(result$stage, 1L)
    expect_true(is.finite(result$z_orig))
    null <- unlist(result[, grep("^z_null_", names(result)), with = FALSE], use.names = FALSE)
    expect_true(all(null == result$z_orig))
    expect_false(result$significant)
    manifest <- prepare_execution(object, "run_discovery_analysis", 2L, 3L)
    proposals <- lapply(manifest$units, function(unit)
        execute_work_unit(object, manifest, unit$id))
    collected <- collect_execution(object, manifest, proposals)
    expect_identical(refactor_plain(get_result(collected, "run_discovery_analysis")),
        refactor_plain(result))
    expect_identical(.Random.seed, state)
})
