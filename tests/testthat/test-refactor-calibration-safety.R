test_that("native producers retain safe zero and endpoint domains", {
    bridge <- refactor_native_bridge(test_path("fixtures", "refactor", "native-bridge.cpp"))
    materialize <- bridge$refactor_native_materialize
    fy <- refactor_private("fisher_yates_sampler")
    hybrid <- refactor_private("hybrid_fisher_iwor_sampler")
    expect_identical(materialize(fy(0L, 0L, 3L)), rep(list(integer()), 3L))
    expect_identical(materialize(fy(5L, 0L, 2L)), rep(list(integer()), 2L))
    expect_identical(materialize(fy(5L, 3L, 0L)), list())
    expect_identical(materialize(hybrid(2L, 1L, 3L, 0L)), list())
    sampler <- refactor_private("sample_crt_indices")
    expect_identical(materialize(sampler(numeric(), 3L)), rep(list(integer()), 3L))
    expect_identical(materialize(sampler(c(0, 1, 0, 1), 3L)), rep(list(c(1L, 3L)), 3L))
    expect_identical(materialize(sampler(c(0, 1), 0L)), list())
    for (probabilities in list(NA_real_, NaN, Inf, -0.1, 1.1))
        expect_error(sampler(probabilities, 0L), "finite and in")
    expect_error(sampler(0.5, -1L), "nonnegative")
    for (args in list(list(-1L, 0L, 0L), list(1L, -1L, 0L),
                     list(1L, 2L, 0L), list(1L, 1L, NA_integer_)))
        expect_error(do.call(fy, args), "nonnegative|must not exceed")
    for (args in list(list(0L, 1L, 1L, 0L), list(1L, 0L, 1L, 0L),
                     list(1L, 2L, 1L, 0L), list(1L, 1L, 1L, -1L),
                     list(.Machine$integer.max, 1L, 1L, 0L),
                     list(1L, 1L, .Machine$integer.max, 0L)))
        expect_error(do.call(hybrid, args), "requires|nonnegative|capacity")
})

test_that("calibration combinations reject unsafe sizes and rates before allocation", {
    enumerate <- refactor_private("iterate_over_combinations")
    sample_groups <- refactor_private("sample_calibration_groups")
    expect_identical(enumerate(3L, 2L, 3L), rbind(c(1L, 2L), c(0L, 2L), c(0L, 1L)))
    for (args in list(list(2L, 3L, 0L), list(2L, 0L, 1L),
                     list(3L, 2L, 2L), list(3L, 2L, NA_integer_)))
        expect_error(do.call(enumerate, args), "group size|combination count|nonnegative")
    for (rate in c(0, NA_real_, NaN, Inf, -1, 1.1))
        expect_error(sample_groups(2L, 1, 6, 4L, 1, 100L, rate), "p_hat")
    # Clamp the double estimate to the finite represented count before casting.
    expect_equal(nrow(sample_groups(2L, 1, 6, 4L, 1, 100L, .Machine$double.xmin)), 6L)
    expect_error(sample_groups(1000L, 1, Inf, 2000L, 1, 100L,
                               .Machine$double.xmin), "Selected.*capacity")
    for (request in c(0, 1.5, Inf, .Machine$integer.max + 1))
        expect_error(sample_groups(2L, request, 6, 4L, 1, 100L, 1), "positive integer")
})

test_that("candidate calibration sampling checks dimensions counts and wide products", {
    sample_pairs <- refactor_private("sample_calibration_pairs")
    call <- function(counts = matrix(1L, 1L, 1L), totals = 2L,
                     groups = matrix(0L, 1L, 1L), genes = 1L, request = 1L) {
        sample_pairs(counts, totals, groups, genes, request, 0L, 0L)
    }
    expect_identical(call()$response_idxs, 1L)
    expect_error(call(totals = integer()), "matching gene dimensions")
    expect_error(call(groups = matrix(integer(), 0L, 1L)), "nonempty")
    expect_error(call(groups = matrix(1L, 1L, 1L)), "zero-based")
    expect_error(call(counts = matrix(NA_integer_, 1L, 1L)), "not NA")
    expect_error(call(counts = matrix(3L, 1L, 1L)), "exceeds its total")
    expect_error(call(counts = matrix(.Machine$integer.max, 1L, 1L),
                      totals = .Machine$integer.max, groups = matrix(c(0L, 0L), 1L)),
                 "treatment count exceeds native")
    # These small inputs represent > INT_MAX candidate pairs, without allocating
    # the former multi-billion-element permutation or executing the old kernel.
    expect_error(call(counts = matrix(1L, 1L, 46341L), totals = rep(2L, 46341L),
                      groups = matrix(0L, 46341L, 1L), genes = 46341L),
                 "Candidate pair count exceeds")
})

test_that("R budget conversion is finite integral and checked before narrowing", {
    total <- refactor_private("checked_resampling_total")
    budget <- refactor_private("compute_empirical_budget")
    expect_identical(total(499L, 0L, 0L), 499L)
    expect_identical(budget(2L, 0.05, 0L, 499L, 0L), 400L)
    for (value in list(NA_integer_, Inf, -1L, 0.5, numeric(), c(1, 2), .Machine$integer.max + 1))
        expect_error(total(value, 0L, 0L), "finite integer")
    expect_error(total(.Machine$integer.max, 1L, 0L), "B1 \\+ B2 \\+ B3")
    expect_error(budget(1L, 0, 1L, 0L, 0L), "Derived B3")
    expect_error(budget(.Machine$integer.max, 0.05, 0L, 0L, 0L), "Derived B3")
    expect_error(budget(1L, 1, 1L, .Machine$integer.max, 0L), "B1 \\+ B2 \\+ B3")
})

test_that("public calibration uses actual NT groups surviving cellwise QC", {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    case <- readRDS(test_path("fixtures", "refactor", "cases.rds"))$low_complement_perm
    object <- suppressMessages(refactor_prepare(inputs, case))
    groups <- object@grna_assignments_raw$indiv_nt_grna_idxs
    expect_true(length(groups) >= 2L)
    removed_group <- names(groups)[[1L]]
    # Remove source cells through the public QC path; metadata intentionally
    # still describes the guide that no longer has an assigned analysis cell.
    object <- suppressMessages(sceptre::run_qc(object,
        response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
        n_nonzero_trt_thresh = 0L, n_nonzero_cntrl_thresh = 0L,
        additional_cells_to_remove = unique(c(inputs$remove_cells, groups[[1L]]))))
    expect_false(removed_group %in% names(object@grna_assignments$indiv_nt_grna_idxs))
    expect_true(removed_group %in% object@grna_target_data_frame$grna_id)
    active <- length(object@grna_assignments$indiv_nt_grna_idxs)
    expect_error(suppressMessages(sceptre::run_calibration_check(object,
        calibration_group_size = active + 1L, n_calibration_pairs = 1L)),
        "active non-targeting gRNA count after QC")
    one_active <- suppressMessages(sceptre::run_qc(object,
        response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
        n_nonzero_trt_thresh = 0L, n_nonzero_cntrl_thresh = 0L,
        additional_cells_to_remove = unique(c(inputs$remove_cells,
            unlist(groups[-length(groups)], use.names = FALSE)))))
    expect_length(one_active@grna_assignments$indiv_nt_grna_idxs, 1L)
    expect_gt(refactor_private("compute_calibration_group_size")(
        one_active@grna_target_data_frame), 1L)
    # Omit group size deliberately: the metadata-derived default is checked
    # against active assignments, rather than silently reduced after QC.
    expect_error(suppressMessages(sceptre::run_calibration_check(one_active,
        n_calibration_pairs = 1L)), "active non-targeting gRNA count after QC")
    for (request in list(0L, -1L, 1.5, NA_real_, Inf, .Machine$integer.max + 1))
        expect_error(suppressMessages(sceptre::run_calibration_check(object,
            calibration_group_size = 1L, n_calibration_pairs = request)),
            "zero negative control pairs|finite integer")
})

test_that("overlapping complement calibration bounds its union bank by cells", {
    bridge <- refactor_native_bridge(test_path("fixtures", "refactor", "native-bridge.cpp"))
    generate <- refactor_private("get_synthetic_permutation_idxs")
    assignments <- list(indiv_nt_grna_idxs = list(nt1 = 1:4, nt2 = 2:5))
    draws <- bridge$refactor_native_materialize(generate(assignments, 7L, TRUE, TRUE, 2L, 5L))
    expect_identical(lengths(draws), rep(5L, 7L))
    expect_true(all(vapply(draws, function(x) identical(sort(x), 0:4), logical(1))))
    expect_error(generate(assignments, 1L, TRUE, TRUE, 3L, 5L), "active non-targeting")
})

test_that("public high-MOI calibration accepts overlapping NT unions", {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    case <- readRDS(test_path("fixtures", "refactor", "cases.rds"))$high_perm
    nt <- inputs$grna_target_data_frame$grna_id[
        inputs$grna_target_data_frame$grna_target == "non-targeting"]
    n <- ncol(inputs$high_grna_matrix)
    common <- seq_len(n - 6L)
    inputs$high_grna_matrix[nt, ] <- 0
    inputs$high_grna_matrix[nt, common] <- 50
    target <- setdiff(rownames(inputs$high_grna_matrix), nt)[[1L]]
    inputs$high_grna_matrix[target, setdiff(seq_len(n), common)] <- 50
    object <- suppressMessages(refactor_prepare(inputs, case))
    object <- suppressMessages(sceptre::run_qc(object,
        response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
        n_nonzero_trt_thresh = 0L, n_nonzero_cntrl_thresh = 0L,
        additional_cells_to_remove = inputs$remove_cells))
    assignments <- object@grna_assignments$indiv_nt_grna_idxs
    expect_true(sum(sort(lengths(assignments), decreasing = TRUE)[1:2]) >
                    length(object@cells_in_use))
    result <- suppressMessages(sceptre::run_calibration_check(object,
        calibration_group_size = 2L, n_calibration_pairs = 1L, print_progress = FALSE))
    expect_equal(nrow(result@calibration_result), 1L)
    pairs <- result@negative_control_pairs
    guides <- strsplit(pairs$grna_group[[1L]], "&", fixed = TRUE)[[1L]]
    selected <- unique(unlist(assignments[guides], use.names = FALSE))
    source <- object@cells_in_use[selected]
    complement <- setdiff(object@cells_in_use, source)
    response <- inputs$response_matrix[pairs$response_id[[1L]], ]
    expect_equal(pairs$n_nonzero_trt, sum(response[source] > 0))
    expect_equal(pairs$n_nonzero_cntrl, sum(response[complement] > 0))
})

test_that("native adapter checks selected prefixes without restricting CRT cardinality", {
    bridge <- refactor_native_bridge(test_path("fixtures", "refactor", "native-bridge.cpp"))
    run <- refactor_private("run_association_test")
    args <- list(y = c(2, 1, 3), mu = c(1, 1, 1), a = c(2, -1, 0),
                 w = rep(1, 3), D = matrix(0, 0L, 3L), trt_idxs = 1L,
                 n_trt = 1L, use_all_cells = FALSE,
                 synthetic_idxs = bridge$refactor_native_pointer(list(1L), 3L),
                 B1 = 1L, B2 = 0L, B3 = 0L, fit_parametric_curve = FALSE,
                 return_resampling_dist = TRUE, side_code = 1L)
    expect_equal(do.call(run, args)$stage, 1)
    changed <- function(...) do.call(run, utils::modifyList(args, list(...)))
    expect_error(changed(y = 1), "aligned cell dimensions")
    expect_error(changed(D = matrix(0, 0L, 2L)), "aligned cell dimensions")
    expect_error(changed(trt_idxs = 0L), "outside.*cell range")
    expect_error(changed(trt_idxs = 4L), "outside.*cell range")
    expect_error(changed(n_trt = 2L), "must equal")
    expect_error(changed(B1 = -1L), "nonnegative")
    expect_error(changed(B1 = .Machine$integer.max, B2 = 1L), "capacity")
    expect_error(changed(B1 = 2L), "bank length")
    empty <- bridge$refactor_native_pointer(list(integer()), 3L)
    expect_error(changed(synthetic_idxs = empty), "shorter.*prefix")
    expect_equal(changed(synthetic_idxs = empty, use_all_cells = TRUE)$stage, 1)
    zero_bank <- refactor_private("fisher_yates_sampler")(3L, 1L, 0L)
    expect_length(changed(synthetic_idxs = zero_bank, B1 = 0L)$resampling_dist, 0L)
    # A hybrid suffix can exceed the smaller current fit view; its selected
    # prefix is valid for N + n_trt cells and must not be rejected by a bank scan.
    hybrid <- refactor_private("hybrid_fisher_iwor_sampler")(2L, 1L, 5L, 1L)
    expect_equal(changed(synthetic_idxs = hybrid)$stage, 1)
    # Short B2 is irrelevant to a stage-1 result, including an unused bank tail.
    expect_equal(changed(B2 = 1L, fit_parametric_curve = TRUE)$stage, 1)
    reached <- bridge$refactor_native_pointer(rep(list(1L), 149L), 3L)
    expect_error(changed(synthetic_idxs = reached, B1 = 50L, B2 = 99L,
                         fit_parametric_curve = TRUE), "at least 100")
})
