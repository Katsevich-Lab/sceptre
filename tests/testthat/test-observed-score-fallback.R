invalid_score_case <- function(a = c(0, 0), w = c(0, 1), D = matrix(0, 1, 2)) {
    list(y = c(1, 3), mu = c(2, 2), a = a, w = w, D = D,
        trt_idxs = 1L, n_trt = 1L, use_all_cells = FALSE,
        draws = rep(list(1L), 3), B1 = 1L, B2 = 1L, B3 = 1L,
        fit_parametric_curve = FALSE, return_resampling_dist = TRUE, side_code = 0L)
}

test_that("invalid observed scores return p one before null ranking, without consuming RNG", {
    bridge <- refactor_native_bridge(test_path("fixtures", "refactor", "native-bridge.cpp"))
    cases <- list(invalid_score_case(), invalid_score_case(a = c(1, 0)),
        invalid_score_case(a = c(-1, 0)), invalid_score_case(w = c(-1, 1)),
        invalid_score_case(a = c(NaN, 0)), invalid_score_case(a = c(Inf, 0), w = c(1, 1)))
    set.seed(319)
    state <- .Random.seed
    for (case in cases) for (side in -1:1) for (crt in c(FALSE, TRUE)) {
        case$return_resampling_dist <- TRUE
        case$side_code <- side
        case$use_all_cells <- crt
        before <- serialize(case, NULL)
        actual <- refactor_native_run_case(case, bridge)
        expect_identical(actual$p, 1)
        expect_false(is.finite(actual$z_orig))
        expect_identical(actual$stage, 0L)
        expect_identical(actual$sn_params, rep(NA_real_, 3))
        expect_identical(actual$resampling_dist, numeric())
        expect_identical(serialize(case, NULL), before)
        expect_identical(.Random.seed, state)
        case$return_resampling_dist <- FALSE
        compact <- refactor_native_run_case(case, bridge)
        expect_identical(compact, actual[names(compact)])
    }
    bad <- invalid_score_case()
    bad$w <- 1
    expect_error(refactor_native_run_case(bad, bridge), "aligned cell dimensions")
})

test_that("stage zero rows survive all output levels and all-invalid worker tables", {
    bridge <- refactor_native_bridge(test_path("fixtures", "refactor", "native-bridge.cpp"))
    fallback <- refactor_native_run_case(invalid_score_case(), bridge)
    case <- invalid_score_case(a = c(1, 0), w = c(1, 1))
    finite <- refactor_native_run_case(case, bridge)
    rows <- data.table::data.table(response_id = c("r1", "r2"), grna_group = "A")
    for (output in 1:3) for (native in list(list(fallback, fallback), list(fallback, finite))) {
        result <- refactor_private("pair_results_from_native")(data.table::copy(rows), native, output)
        expect_identical(nrow(result), 2L)
        expect_identical(result$response_id, rows$response_id)
        expect_identical(result$p_value[1], 1)
        if (output >= 2L) expect_identical(result$stage[1], 0L)
        if (output == 3L) expect_true(all(is.na(result[1, grep("^z_null_", names(result)), with = FALSE])))
    }
})

invalid_score_screen <- function() {
    inputs <- refactor_make_inputs()
    case <- refactor_case_definitions()$low_complement_perm
    # All A-treated cells, and only those cells, are in batch A. r3 remains
    # variable, passes pair QC, and has a nonfinite score after adjustment.
    inputs$extra_covariates$batch <- factor(ifelse(
        colSums(inputs$low_grna_matrix[1:2, , drop = FALSE]) > 0, "A", "B"))
    inputs$remove_cells <- integer()
    refactor_capture(refactor_prepare(inputs, case))$value
}

test_that("public confounded pair is retained as p one in direct and distributed families", {
    object <- invalid_score_screen()
    before <- serialize(object, NULL)
    set.seed(401)
    state <- .Random.seed
    for (output in 1:3) {
        direct <- refactor_capture(run_discovery_analysis(object,
            print_progress = FALSE, output_amount = output))$value
        result <- get_result(direct, "run_discovery_analysis")
        invalid <- result[result$response_id == "r3" & result$grna_target == "A", ]
        expect_identical(nrow(result), 5L)
        expect_identical(nrow(invalid), 1L)
        expect_true(invalid$pass_qc)
        expect_identical(invalid$p_value, 1)
        expect_false(invalid$significant)
        if (output >= 2L) {
            expect_identical(invalid$stage, 0L)
            expect_false(is.finite(invalid$z_orig))
        }
        manifest <- prepare_execution(object, "run_discovery_analysis", 2L, output)
        proposals <- lapply(manifest$units, function(unit)
            execute_work_unit(object, manifest, unit$id))
        collected <- collect_execution(object, manifest, rev(proposals))
        expect_identical(refactor_plain(get_result(collected, "run_discovery_analysis")),
            refactor_plain(result))
        expect_identical(.Random.seed, state)
    }
    expect_identical(serialize(object, NULL), before)
})
