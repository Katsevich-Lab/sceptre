experimental_logistic_fixture <- function(n = 500L, p = 6L, seed = 1006L) {
    withr::with_seed(seed, {
        D <- cbind(1, matrix(stats::rnorm(n * (p - 1L)), n, p - 1L))
        colnames(D) <- c(
            "(Intercept)",
            if (p > 1L) paste0("z", seq_len(p - 1L))
        )
        beta <- c(-1.2, if (p > 1L) seq(-.35, .35, length.out = p - 1L))
        y <- stats::rbinom(n, 1, stats::plogis(D %*% beta))
        list(D = D, y = y, idx = which(y == 1L))
    })
}

test_that("candidate targets the logistic likelihood without replacing reference", {
    set.seed(6107)
    for (p in c(1L, 6L, 19L)) {
        f <- experimental_logistic_fixture(700L, p, 1000L + p)
        before <- .Random.seed
        ref <- stats::glm.fit(f$D, f$y, family = stats::binomial())
        tight <- stats::glm.fit(
            f$D,
            f$y,
            family = stats::binomial(),
            control = stats::glm.control(epsilon = 1e-12, maxit = 100L)
        )
        ordinary <- fit_perturbation_model(f$idx, f$D)
        candidate <- fit_perturbation_model(f$idx, f$D, "fast_logistic")
        expect_identical(ordinary$probabilities, ref$fitted.values)
        expect_identical(candidate$fit$diagnostics$method_used, "fast_logistic")
        expect_true(candidate$fit$diagnostics$converged)
        expect_lt(max(abs(candidate$probabilities - tight$fitted.values)), 1e-7)
        expect_lt(
            max(abs(candidate$fit$fitted_coefs - tight$coefficients)),
            1e-6
        )
        expect_identical(
            candidate$probabilities,
            reconstruct_perturbation_probabilities(candidate$fit, f$D)
        )
        expect_identical(.Random.seed, before)
    }
})

test_that("unsupported and difficult designs fall back with reference warnings", {
    f <- experimental_logistic_fixture()
    cases <- list(
        list(D = f$D[, -1L, drop = FALSE], idx = f$idx),
        list(D = cbind(f$D, duplicate = f$D[, 2]), idx = f$idx),
        list(D = f$D, idx = integer()),
        list(D = f$D, idx = seq_len(nrow(f$D))),
        list(
            D = cbind('(Intercept)' = 1, z = seq(-1, 1, length.out = 500)),
            idx = 251:500
        )
    )
    for (case in cases) {
        warnings <- character()
        ref <- withCallingHandlers(
            fit_perturbation_glm(case$idx, case$D),
            warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        )
        actual_warnings <- character()
        result <- withCallingHandlers(
            fit_perturbation_model(case$idx, case$D, "fast_logistic"),
            warning = function(w) {
                actual_warnings <<- c(actual_warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        )
        expect_identical(result$probabilities, ref$fitted.values)
        expect_identical(actual_warnings, warnings)
        if (!is.null(result$fit)) {
            expect_identical(result$fit$diagnostics$method_used, "glm.fit")
            expect_true(nzchar(result$fit$diagnostics$fallback_reason))
        }
    }
    expect_error(fit_perturbation_model(f$idx, f$D, "typo"), "method")
    # Repeated treatment positions retain the original indicator semantics.
    expect_identical(
        fit_perturbation_model(c(f$idx, f$idx), f$D, "fast_logistic"),
        fit_perturbation_model(f$idx, f$D, "fast_logistic")
    )
})

test_that("selector survives serialization and preserves historical default identities", {
    o <- suppressMessages(refactor_prepare(
        refactor_make_inputs(),
        refactor_case_definitions()$high_complement_crt
    ))
    old <- o
    attr(old@spec, "perturbation_fit_method") <- NULL
    restored <- migrate_sceptre_object(unserialize(serialize(old, NULL)))
    expect_identical(
        analysis_spec_snapshot(restored),
        analysis_spec_snapshot(o)
    )
    expect_identical(
        perturbation_cache_identity(restored),
        perturbation_cache_identity(o)
    )
    candidate <- set_perturbation_fit_method(o, "fast_logistic")
    expect_true(stage_complete(candidate, "run_qc"))
    expect_identical(
        get_perturbation_fit_method(unserialize(serialize(candidate, NULL))),
        "fast_logistic"
    )
    expect_false(identical(
        perturbation_cache_identity(candidate),
        perturbation_cache_identity(o)
    ))
    expect_identical(
        analysis_spec_snapshot(set_perturbation_fit_method(candidate)),
        analysis_spec_snapshot(o)
    )
    expect_error(set_perturbation_fit_method(o, character()), "scalar")
    expect_error(set_perturbation_fit_method(o, "other"), "method")
})

test_that("candidate cold and warm CRT preserve results draws and method certificate", {
    for (control in c("complement", "nt_cells")) {
        case <- refactor_case_definitions()$high_complement_crt
        if (control == "nt_cells") {
            case$moi <- "low"
            case$control <- control
        }
        o <- suppressMessages(refactor_prepare(refactor_make_inputs(), case))
        o <- set_perturbation_fit_method(o, "fast_logistic")
        run <- function(x) {
            withr::with_seed(
                817L,
                suppressMessages(run_discovery_analysis(
                    x,
                    print_progress = FALSE,
                    parallel = FALSE,
                    output_amount = 3L
                ))
            )
        }
        cold <- run(o)
        warm <- run(cold)
        expect_identical(get_result(cold), get_result(warm))
        expect_identical(
            cold@fits@perturbation_fits,
            warm@fits@perturbation_fits
        )
        expect_true(length(cold@fits@perturbation_fits) > 0L)
        expect_identical(
            cold@results@records$run_discovery_analysis$provenance$perturbation_fit_method,
            perturbation_fit_method_id("fast_logistic")
        )
        ref <- set_perturbation_fit_method(cold)
        expect_length(ref@fits@perturbation_fits, 0L)
        expect_false(result_record_is_current(ref, "run_discovery_analysis"))
    }
})


test_that("historical schema-three physical objects retain completed results", {
    old <- readRDS(test_path(
        "fixtures",
        "saved-objects",
        "sceptre-schema3-list-plans-2a0938e.rds"
    ))
    expect_null(attr(old@spec, "perturbation_fit_method", exact = TRUE))
    before <- serialize(old, NULL)
    current <- migrate_sceptre_object(old)
    expect_identical(get_perturbation_fit_method(current), "glm.fit")
    expect_true(result_record_is_current(current, "run_discovery_analysis"))
    expect_identical(
        current@results@records$run_discovery_analysis$identity,
        old@results@records$run_discovery_analysis$identity
    )
    expect_identical(serialize(old, NULL), before)
})

test_that("distributed manifests reject a preselected experimental object", {
    o <- suppressMessages(refactor_prepare(
        refactor_make_inputs(),
        refactor_case_definitions()$high_complement_crt
    ))
    o <- set_perturbation_fit_method(o, "fast_logistic")
    expect_error(
        prepare_execution(o, "run_discovery_analysis"),
        "direct R execution"
    )
})

test_that("malformed native prepared payloads fail before intercept access", {
    bad <- list(
        X = matrix(numeric(), 2, 0),
        centers = numeric(),
        scales = numeric(),
        column_sums = numeric(),
        crossprod = matrix(numeric(), 0, 0),
        conditioned = TRUE,
        valid = FALSE,
        rank = 0L,
        reason = "invalid"
    )
    expect_error(
        fit_grna_logistic_prepared_cpp(bad, 1L),
        "two rows and one column"
    )
})


test_that("the selector is valid before analysis parameters are set", {
    input <- refactor_make_inputs()
    o <- suppressMessages(refactor_import(
        input,
        refactor_case_definitions()$high_complement_crt
    ))
    candidate <- set_perturbation_fit_method(o, "fast_logistic")
    expect_true(validObject(candidate))
    expect_identical(get_perturbation_fit_method(candidate), "fast_logistic")
})


test_that("CRT calibration uses the selected fitter transiently", {
    o <- suppressMessages(refactor_prepare(
        refactor_make_inputs(),
        refactor_case_definitions()$high_complement_crt
    ))
    o <- set_perturbation_fit_method(o, "fast_logistic")
    fit <- fit_fast_perturbation_glm
    calls <- 0L
    testthat::local_mocked_bindings(fit_fast_perturbation_glm = function(...) {
        calls <<- calls + 1L
        fit(...)
    })
    before <- o@fits@perturbation_fits
    testthat::local_mocked_bindings(lookup_perturbation_fit = function(...) {
        stop("Calibration must not look up logistic caches")
    })
    completed <- withr::with_seed(
        817L,
        suppressMessages(run_calibration_check(
            o,
            n_calibration_pairs = 3L,
            calibration_group_size = 2L,
            print_progress = FALSE,
            parallel = FALSE
        ))
    )
    expect_gt(calls, 0L)
    expect_identical(completed@fits@perturbation_fits, before)
    expect_true(result_record_is_current(completed, "run_calibration_check"))
})
