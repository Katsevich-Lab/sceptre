test_that("prepared native inputs satisfy declared stages and finite-parameter contracts", {
    inputs <- readRDS(test_path("fixtures", "refactor", "native-inputs.rds"))
    bridge <- refactor_native_bridge(test_path("fixtures", "refactor", "native-bridge.cpp"))
    for (name in names(inputs$cases)) {
        case <- inputs$cases[[name]]
        before <- serialize(case, NULL)
        actual <- refactor_native_run_case(case, bridge)
        expect_identical(actual$stage, case$expected_stage, info = name)
        expect_identical(all(is.finite(actual$sn_params)), case$expected_sn, info = name)
        fc <- mean(case$y[case$trt_idxs]) / mean(case$mu[case$trt_idxs])
        se <- sqrt(sum((case$y[case$trt_idxs] - fc * case$mu[case$trt_idxs])^2) /
                   sum(case$mu[case$trt_idxs])^2)
        expect_equal(actual$fc, fc, tolerance = 1e-8, info = name)
        expect_equal(actual$se, se, tolerance = 1e-8, info = name)
        expect_identical(serialize(case, NULL), before)
    }
})

test_that("native fallback uses its selected segment including zero B3", {
    inputs <- readRDS(test_path("fixtures", "refactor", "native-inputs.rds"))
    bridge <- refactor_native_bridge(test_path("fixtures", "refactor", "native-bridge.cpp"))
    for (name in c("threshold_right", "failed_sn_stage3", "failed_sn_B3_zero", "empirical_B3_zero")) {
        case <- inputs$cases[[name]]
        actual <- refactor_native_run_case(case, bridge)
        selected <- if (case$B3 > 0L) {
            case$B1 + case$B2 + seq_len(case$B3)
        } else if (case$fit_parametric_curve) {
            case$B1 + seq_len(case$B2)
        } else seq_len(case$B1)
        # These prepared cases have one selected cell, unit weight and D = 0.
        selected_scores <- vapply(case$draws[selected], function(draw) case$a[draw + 1L], numeric(1))
        expect_equal(actual$resampling_dist, selected_scores, tolerance = 1e-8, info = name)
    }
})

test_that("native distribution can be omitted without changing the inference", {
    fixture <- test_path("fixtures", "refactor")
    case <- readRDS(file.path(fixture, "native-inputs.rds"))$cases$accepted_sn_right
    bridge <- refactor_native_bridge(file.path(fixture, "native-bridge.cpp"))
    full <- refactor_native_run_case(case, bridge)
    case$return_resampling_dist <- FALSE
    compact <- refactor_native_run_case(case, bridge)
    expect_identical(names(compact), setdiff(names(full), "resampling_dist"))
    expect_equal(compact, full[names(compact)], tolerance = 1e-8)
})

test_that("production samplers retain draw counts, ranges and membership laws", {
    fixture <- test_path("fixtures", "refactor")
    inputs <- readRDS(file.path(fixture, "native-inputs.rds"))
    bridge <- refactor_native_bridge(file.path(fixture, "native-bridge.cpp"))
    set.seed(719L)
    state <- .Random.seed
    active <- refactor_native_candidate_samplers()
    for (name in names(active)) {
        config <- inputs$samplers[[name]]
        draws <- refactor_native_sample(config, active[[name]], bridge)
        expect_identical(refactor_native_sample(config, active[[name]], bridge), draws,
                         info = paste(name, "repeated-call determinism"))
        expect_identical(length(draws), config$args$B, info = name)
        expect_true(all(vapply(draws, function(draw) is.integer(draw) && !anyNA(draw) &&
                                   !anyDuplicated(draw), logical(1))), info = name)
        if (grepl("^crt", name)) {
            expect_true(all(vapply(draws, function(draw) all(draw >= 0L & draw < 5L) &&
                                       !(0L %in% draw) && 4L %in% draw, logical(1))), info = name)
        } else if (name == "hybrid") {
            expect_true(all(lengths(draws) == config$args$M))
            for (size in config$args$m:config$args$M) {
                expect_true(all(vapply(draws, function(draw) all(draw[seq_len(size)] >= 0L &
                                       draw[seq_len(size)] < config$args$N + size), logical(1))))
            }
        } else {
            expect_true(all(lengths(draws) == config$args$M))
            expect_true(all(unlist(draws) >= 0L & unlist(draws) < config$args$n_tot))
        }
    }
    expect_identical(.Random.seed, state)
})

test_that("test-only bank construction rejects invalid native indices", {
    bridge <- refactor_native_bridge(test_path("fixtures", "refactor", "native-bridge.cpp"))
    expect_error(bridge$refactor_native_pointer(list(-1L), 5L), "zero-based cell range")
    expect_error(bridge$refactor_native_pointer(list(5L), 5L), "zero-based cell range")
    expect_error(bridge$refactor_native_pointer(list(NA_integer_), 5L), "zero-based cell range")
    expect_error(bridge$refactor_native_pointer(list(c(1L, 1L)), 5L), "must be unique")
    expect_error(bridge$refactor_native_pointer(list(1.5), 5L), "integer vectors")
})
