fit_cache_crt_fixture <- function(control = "complement") {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    cases <- readRDS(test_path("fixtures", "refactor", "cases.rds"))
    case <- if (control == "nt_cells") cases$low_nt_perm else cases$high_complement_crt
    case$mechanism <- "crt"
    suppressMessages(refactor_prepare(inputs, case))
}

fit_cache_discover <- function(object) {
    withr::with_seed(817L, refactor_capture({
        object <- run_discovery_analysis(object, print_progress = FALSE,
            parallel = FALSE, output_amount = 3L)
        list(object = object, rng = .Random.seed)
    }))
}

fit_cache_group_inputs <- function(object, group) {
    pairs <- dplyr::filter(object@preparation@discovery_pairs_with_info, pass_qc)
    context <- prepare_analysis_context(object, "discovery_analysis", pairs)
    if (context$shared_response_fit) {
        list(source_cells = context$common_view@source_cells,
            design = context$design, trt = make_treatment_index_lookup(context$index)(group)$trt_idxs)
    } else {
        view <- make_fit_view(context$index, group)
        list(source_cells = view@source_cells,
            design = context$design[view@fit_to_analysis, , drop = FALSE], trt = view@trt_fit_pos)
    }
}

test_that("both CRT control modes reuse real fits while drawing fresh full banks", {
    original_fit <- fit_perturbation_glm
    calls <- 0L
    testthat::local_mocked_bindings(fit_perturbation_glm = function(...) {
        calls <<- calls + 1L
        original_fit(...)
    })
    for (control in c("complement", "nt_cells")) {
        object <- fit_cache_crt_fixture(control)
        expect_false(object@spec@run_permutations)
        expect_identical(object@spec@control_group_complement, control == "complement")
        expect_identical(object@data@low_moi, control == "nt_cells")
        groups <- unique(dplyr::filter(object@preparation@discovery_pairs_with_info, pass_qc)$grna_group)
        before <- calls
        cold <- fit_cache_discover(object)
        fitted <- cold$value$object
        expect_identical(calls - before, length(groups))
        expect_setequal(names(fitted@fits@perturbation_fits), groups)
        before <- calls
        warm <- fit_cache_discover(fitted)
        expect_identical(calls, before)
        expect_identical(warm$value$rng, cold$value$rng)
        expect_identical(refactor_plain(get_result(warm$value$object, "run_discovery_analysis")),
            refactor_plain(get_result(fitted, "run_discovery_analysis")))
        expect_identical(warm$warnings, character())
        expect_identical(cold$warnings, character())
        expect_identical(warm$value$object@fits@perturbation_fits, fitted@fits@perturbation_fits)
    }
})

test_that("NT perturbation fits bind actual group cells without repeating common vectors", {
    object <- fit_cache_discover(fit_cache_crt_fixture("nt_cells"))$value$object
    expect_true(object@data@low_moi)
    expect_false(object@spec@control_group_complement)
    expect_false(object@spec@run_permutations)
    expect_identical(object@fits@perturbation_identity$cells, object@preparation@analysis_cells)
    expect_true(is.list(object@fits@perturbation_identity$runtime))
    pairs <- dplyr::filter(object@preparation@discovery_pairs_with_info, pass_qc)
    context <- prepare_analysis_context(object, "discovery_analysis", pairs)
    nt_source <- object@preparation@analysis_cells[context$index$nt_to_analysis]
    for (group in names(object@fits@perturbation_fits)) {
        fit <- object@fits@perturbation_fits[[group]]
        inputs <- fit_cache_group_inputs(object, group)
        expect_identical(tail(inputs$source_cells, length(nt_source)), nt_source)
        indicator <- integer(nrow(inputs$design))
        indicator[inputs$trt] <- 1L
        reference <- stats::glm.fit(x = inputs$design, y = indicator, family = stats::binomial())
        expect_identical(reconstruct_perturbation_probabilities(fit, inputs$design), reference$fitted.values)
        expect_setequal(names(fit), c("fitted_coefs", "diagnostics", "identity"))
        expect_length(fit$fitted_coefs, ncol(inputs$design))
        expect_setequal(names(fit$identity), c("cells", "design", "treatment"))
        expect_true(all(lengths(fit$identity) == 1L))
        for (signature in fit$identity) expect_match(signature, "^[[:xdigit:]]{64}$")
        expect_setequal(names(fit$diagnostics), c("rank", "converged", "boundary", "iterations", "warnings"))
    }
})

test_that("group labels cannot substitute for the actual perturbation fitting inputs", {
    object <- fit_cache_discover(fit_cache_crt_fixture("nt_cells"))$value$object
    group <- names(object@fits@perturbation_fits)[1L]
    inputs <- fit_cache_group_inputs(object, group)
    identity <- perturbation_fit_identity(inputs$source_cells, inputs$design, inputs$trt)
    expect_identical(lookup_perturbation_fit(object@fits@perturbation_fits, group, identity),
        object@fits@perturbation_fits[[group]])
    changed <- list(
        perturbation_fit_identity(rev(inputs$source_cells), inputs$design, inputs$trt),
        perturbation_fit_identity(inputs$source_cells,
            replace(inputs$design, 1L, inputs$design[1L] + 1), inputs$trt),
        perturbation_fit_identity(inputs$source_cells, inputs$design,
            c(inputs$trt, setdiff(seq_len(nrow(inputs$design)), inputs$trt)[1L])))
    for (scope in changed) expect_null(lookup_perturbation_fit(object@fits@perturbation_fits, group, scope))
    expect_identical(perturbation_fit_identity(inputs$source_cells, inputs$design,
        c(inputs$trt, inputs$trt)), identity)
})

test_that("assignment and eligibility changes invalidate the two model families independently", {
    object <- fit_cache_discover(fit_cache_crt_fixture())$value$object
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    response <- object@fits@response_fits
    perturbation <- object@fits@perturbation_fits
    expect_gt(length(response), 0L)
    expect_gt(length(perturbation), 0L)
    reporting <- set_analysis_parameters(object,
        discovery_pairs = inputs$discovery_pairs, positive_control_pairs = inputs$positive_control_pairs,
        formula_object = object@spec@formula_object, side = "left", control_group = "complement",
        resampling_mechanism = "crt", resampling_approximation = "skew_normal", multiple_testing_alpha = 0.2)
    expect_identical(reporting@fits@response_fits, response)
    expect_identical(reporting@fits@perturbation_fits, perturbation)
    reassigned <- suppressMessages(assign_grnas(object, method = "thresholding", threshold = 45,
        print_progress = FALSE))
    expect_identical(reassigned@fits@response_fits, response)
    expect_length(reassigned@fits@perturbation_fits, 0L)
    expect_length(reassigned@fits@perturbation_identity, 0L)
    eligibility <- suppressMessages(run_qc(object, response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1), n_nonzero_trt_thresh = 30L,
        n_nonzero_cntrl_thresh = 7L, additional_cells_to_remove = inputs$remove_cells))
    expect_identical(eligibility@preparation@analysis_cells, object@preparation@analysis_cells)
    expect_false(identical(eligibility@preparation@discovery_pairs_with_info$pass_qc,
        object@preparation@discovery_pairs_with_info$pass_qc))
    expect_identical(eligibility@fits@response_fits, response)
    expect_identical(eligibility@fits@perturbation_fits, perturbation)
})

test_that("calibration fits afresh and preserves completed discovery perturbation fits", {
    object <- fit_cache_discover(fit_cache_crt_fixture())$value$object
    saved <- object@fits@perturbation_fits
    identity <- object@fits@perturbation_identity
    original_fit <- fit_perturbation_glm
    calls <- 0L
    testthat::local_mocked_bindings(fit_perturbation_glm = function(...) {
        calls <<- calls + 1L
        original_fit(...)
    })
    for (run in 1:2) {
        before <- calls
        object <- suppressMessages(run_calibration_check(object, n_calibration_pairs = 4L,
            calibration_group_size = 1L, print_progress = FALSE, parallel = FALSE))
        expect_gt(calls, before)
        expect_identical(object@fits@perturbation_fits, saved)
        expect_identical(object@fits@perturbation_identity, identity)
    }
})

test_that("a failed CRT worker commits neither newly fitted model family", {
    object <- fit_cache_crt_fixture()
    before <- serialize(object, NULL)
    calls <- 0L
    original_fit <- fit_perturbation_glm
    testthat::local_mocked_bindings(fit_perturbation_glm = function(...) {
        calls <<- calls + 1L
        original_fit(...)
    }, sample_crt_indices = function(...) stop("Controlled CRT worker failure"))
    expect_error(suppressMessages(run_discovery_analysis(object, print_progress = FALSE,
        parallel = FALSE)), "Controlled CRT worker failure")
    expect_identical(calls, 1L)
    expect_identical(serialize(object, NULL), before)
    expect_length(object@fits@response_fits, 0L)
    expect_length(object@fits@perturbation_fits, 0L)
    expect_false(workflow_status(object)[["run_discovery_analysis"]])
})

test_that("public NT CRT refits incompatible complement entries despite matching labels", {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    complement <- fit_cache_crt_fixture()
    complement <- suppressMessages(run_qc(complement, response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1), n_nonzero_trt_thresh = 7L,
        n_nonzero_cntrl_thresh = 7L, additional_cells_to_remove = c(inputs$remove_cells, 25L)))
    complement <- fit_cache_discover(complement)$value$object
    target <- fit_cache_crt_fixture("nt_cells")
    expect_true(target@data@low_moi)
    expect_false(target@spec@control_group_complement)
    expect_false(target@spec@run_permutations)
    expect_gt(length(complement@fits@perturbation_fits), 0L)
    expect_false(identical(complement@fits@perturbation_identity, perturbation_cache_identity(target)))
    reference <- fit_cache_discover(target)
    groups <- names(reference$value$object@fits@perturbation_fits)
    expect_setequal(names(complement@fits@perturbation_fits), groups)
    original_fit <- fit_perturbation_glm
    calls <- 0L
    testthat::local_mocked_bindings(fit_perturbation_glm = function(...) {
        calls <<- calls + 1L
        original_fit(...)
    })
    for (force_common_identity in c(FALSE, TRUE)) {
        transplanted <- target
        transplanted@fits@perturbation_fits <- complement@fits@perturbation_fits
        transplanted@fits@perturbation_identity <- if (force_common_identity) {
            perturbation_cache_identity(target)
        } else complement@fits@perturbation_identity
        before <- calls
        actual <- fit_cache_discover(transplanted)
        expect_identical(calls - before, length(groups))
        expect_identical(actual$value$rng, reference$value$rng)
        expect_identical(refactor_plain(get_result(actual$value$object, "run_discovery_analysis")),
            refactor_plain(get_result(reference$value$object, "run_discovery_analysis")))
        expect_identical(actual$value$object@fits@perturbation_fits,
            reference$value$object@fits@perturbation_fits)
    }
})
