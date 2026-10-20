run_plan_fixture <- function(case = "low_complement_perm") {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    cases <- readRDS(test_path("fixtures", "refactor", "cases.rds"))
    suppressMessages(refactor_prepare(inputs, cases[[case]]))
}

test_that("run budgets and selected families are transient and immutable", {
    object <- run_plan_fixture()
    before <- serialize(object, NULL)
    plan <- make_run_plan(object, "run_discovery_analysis")
    expect_identical(serialize(object, NULL), before)
    expect_s4_class(plan, "RunPlan")
    expect_true(validObject(plan))
    expect_identical(plan@B1, 499L)
    expect_identical(plan@B2, 0L)
    expect_identical(plan@B3, compute_empirical_budget(
        max(pair_count(object, "discovery", passing = TRUE),
            pair_count(object, "positive_control", passing = TRUE)),
        object@spec@multiple_testing_alpha, object@spec@side_code, 499L, 0L))
    expect_false(any(c("B1", "B2", "B3", "negative_control_pairs",
        "calibration_group_size", "n_calibration_pairs", "requested") %in%
        methods::slotNames(object@spec)))
    expect_identical(plan@pairs, object@preparation@discovery_pairs_with_info)
    expect_silent(validate_run_plan(object, plan))
    changed <- plan
    changed@B3 <- plan@B3 + 1L
    expect_error(validate_run_plan(object, changed), "no longer matches")
    changed <- plan
    changed@pairs <- plan@pairs[nrow(plan@pairs):1, , drop = FALSE]
    expect_error(validate_run_plan(object, changed), "no longer matches")
})

test_that("run plans bind actual data, design, ordered cells and assignments", {
    original <- run_plan_fixture()
    plan <- make_run_plan(original, "run_discovery_analysis")
    for (mutation in c("data", "design", "cells", "assignments", "alpha")) {
        object <- original
        if (mutation == "data") object@data@identity <- "different-content"
        if (mutation == "design") object@preparation@design_matrix[1, 1] <- 17
        if (mutation == "cells") object@preparation@analysis_cells <- rev(object@preparation@analysis_cells)
        if (mutation == "assignments") object@preparation@analysis_groups$grna_group_idxs[[1]] <- integer()
        if (mutation == "alpha") object@spec@multiple_testing_alpha <- 0.2
        expect_error(validate_run_plan(object, plan), "no longer matches", info = mutation)
    }
})

test_that("analysis contexts use explicit run budgets with the full requested family", {
    object <- run_plan_fixture()
    plan <- make_run_plan(object, "discovery_analysis", budgets = c(B1 = 9L, B2 = 0L, B3 = 0L))
    expect_true(any(!plan@pairs$pass_qc))
    eligible <- dplyr::filter(plan@pairs, pass_qc)
    context <- prepare_analysis_context(object, "discovery_analysis", eligible, plan = plan)
    expect_identical(context$test_options$B1, 9L)
    expect_identical(context$test_options$total_resamples, 9L)
    expect_identical(context$pairs, eligible)
    expect_error(prepare_analysis_context(object, "discovery_analysis",
        eligible[0, , drop = FALSE], plan = plan), "eligible family")
    expect_error(prepare_analysis_context(object, "power_check", eligible, plan = plan), "Analysis type")
})

test_that("calibration preparation returns a plan without persisting execution settings", {
    object <- run_plan_fixture()
    before <- serialize(object, NULL)
    prepared <- suppressMessages(run_calibration_check_pt_1(object,
        n_calibration_pairs = 2L, calibration_group_size = 1L))
    expect_identical(serialize(prepared$object, NULL), before)
    expect_named(prepared, c("object", "plan"))
    expect_identical(prepared$plan@n_calibration_pairs, 2L)
    expect_identical(prepared$plan@calibration_group_size, 1L)
    expect_equal(nrow(prepared$plan@pairs), 2L)
    expect_null(prepared$object@results@records$run_calibration_check)
    expect_silent(validate_run_plan(prepared$object, prepared$plan))
})

test_that("completed records own the table, frozen choices and realized run plan", {
    object <- run_plan_fixture()
    plan <- make_run_plan(object, "run_discovery_analysis")
    table <- dplyr::rename(plan@pairs, grna_target = grna_group)
    table$p_value <- 0.4
    table$significant <- FALSE
    completed <- commit_result_record(object, "run_discovery_analysis", table,
        plan@pairs, output_amount = 1, plan = plan)
    expect_identical(methods::slotNames(completed@results), c("records", "import_metadata"))
    record <- completed@results@records$run_discovery_analysis
    expect_false(any(c("complete", "status", "calibration_result", "discovery_result") %in% names(record)))
    expect_identical(record$table, table)
    expect_identical(record$provenance$plan, plan)
    expect_identical(record$provenance$spec, analysis_spec_snapshot(object))
    expect_true(result_record_is_current(completed, "run_discovery_analysis"))
    expect_identical(get_result(completed, "run_discovery_analysis"), table)
    updated <- invalidate_result_records(completed, "actual design changed")
    expect_identical(updated@results@records$run_discovery_analysis$invalidated_by, "actual design changed")
    expect_error(get_result(updated, "run_discovery_analysis"), "No completed current")
    expect_warning(old <- get_result(updated, "run_discovery_analysis", historical = TRUE), "historical")
    expect_identical(old$p_value, rep(0.4, nrow(table)))
    rerun <- commit_result_record(updated, "run_discovery_analysis",
        transform(table, p_value = 0.7), plan@pairs, 1, plan = plan)
    expect_identical(rerun@results@records$run_discovery_analysis$history[[1]]$table, table)
    expect_true(result_record_is_current(rerun, "run_discovery_analysis"))
})

test_that("failed reruns preserve committed records and prepared plans", {
    object <- run_plan_fixture()
    prepared <- suppressMessages(run_calibration_check_pt_1(object,
        n_calibration_pairs = 1L, calibration_group_size = 1L))
    table <- dplyr::rename(prepared$plan@pairs, grna_target = grna_group)
    table$p_value <- 0.3
    table$significant <- FALSE
    completed <- commit_result_record(object, "run_calibration_check", table,
        prepared$plan@pairs, 1, plan = prepared$plan)
    before <- serialize(completed, NULL)
    testthat::local_mocked_bindings(analyze_association_pairs = function(...) stop("worker failed"))
    expect_error(run_calibration_check_pt_2(completed, prepared$plan,
        print_progress = FALSE), "worker failed")
    expect_identical(serialize(completed, NULL), before)
    expect_identical(get_result(completed, "run_calibration_check"), table)
})


test_that("committing a result requires valid preparation and the complete family", {
    object <- run_plan_fixture()
    plan <- make_run_plan(object, "run_discovery_analysis")
    table <- dplyr::rename(plan@pairs, grna_target = grna_group)
    table$p_value <- 0.5
    table$significant <- FALSE
    expect_error(commit_result_record(object, "run_discovery_analysis", table[-1, ],
        plan@pairs, 1, plan), "complete run-plan hypothesis family")
    expect_error(commit_result_record(object, "run_discovery_analysis", rbind(table, table[1, ]),
        plan@pairs, 1, plan), "complete run-plan hypothesis family")
    subset_plan <- make_run_plan(object, "run_discovery_analysis", plan@pairs[-1, ])
    expect_error(commit_result_record(object, "run_discovery_analysis", table[-1, ],
        subset_plan@pairs, 1, subset_plan), "complete prepared hypothesis family")
    object@preparation@stages$pair_qc <- NULL
    expect_error(commit_result_record(object, "run_discovery_analysis", table,
        plan@pairs, 1, plan), "run_qc")
})


test_that("a complete result requires full-family adjustment rather than shard decisions", {
    object <- run_plan_fixture()
    plan <- make_run_plan(object, "run_discovery_analysis")
    table <- dplyr::rename(plan@pairs, grna_target = grna_group)
    table$p_value <- rep(1, nrow(table))
    table$p_value[1] <- 0.04
    # A singleton shard rejects 0.04 at alpha 0.1; the complete five-pair family
    # does not. Such shard decisions must not be committed as adjusted output.
    table$significant <- table$p_value < object@spec@multiple_testing_alpha
    expect_error(commit_result_record(object, "run_discovery_analysis", table,
        plan@pairs, 1, plan), "adjusted across the complete hypothesis family")
    table$significant <- stats::p.adjust(table$p_value,
        method = object@spec@multiple_testing_method) < object@spec@multiple_testing_alpha
    expect_silent(commit_result_record(object, "run_discovery_analysis", table,
        plan@pairs, 1, plan))
    table$p_value[1] <- Inf
    expect_error(commit_result_record(object, "run_discovery_analysis", table,
        plan@pairs, 1, plan), "numeric probabilities")
})


test_that("run-plan constructors reject hypotheses absent from prepared data", {
    object <- run_plan_fixture()
    pairs <- object@preparation@discovery_pairs_with_info
    pairs$response_id[1] <- "not-an-imported-response"
    expect_error(make_run_plan(object, "run_discovery_analysis", pairs), "imported screen")
    pairs <- object@preparation@discovery_pairs_with_info
    pairs$grna_group[which(pairs$pass_qc)[1]] <- "unprepared-guide"
    expect_error(make_run_plan(object, "run_discovery_analysis", pairs), "prepared assignments")
    prepared <- suppressMessages(run_calibration_check_pt_1(object,
        n_calibration_pairs = 1L, calibration_group_size = 1L))
    pairs <- prepared$plan@pairs
    pairs$grna_group[1] <- "unprepared-NT"
    expect_error(make_run_plan(object, "run_calibration_check", pairs), "prepared non-targeting")
})


test_that("run plans validate canonical count IDs independently of display names", {
    object <- run_plan_fixture()
    store <- get_response_matrix(object)
    canonical <- rownames(store)
    object@data@response_names <- paste0("display-symbol-", seq_along(canonical))
    expect_false(any(object@data@response_names %in% canonical))
    plan <- make_run_plan(object, "run_discovery_analysis")
    expect_true(all(plan@pairs$response_id %in% canonical))
    pairs <- plan@pairs
    pairs$response_id[1] <- object@data@response_names[1]
    expect_error(make_run_plan(object, "run_discovery_analysis", pairs), "imported screen")

    descriptor <- list(path = "/private/tmp/closed-plan-test.h5",
        content_md5 = paste(rep("a", 32), collapse = ""),
        dimensions = as.integer(dim(store)), dimnames = dimnames(store), integer_id = 1L)
    object@data@response_matrix <- methods::new("BackedCountStore",
        descriptor = descriptor, handle = NULL)
    expect_identical(rownames(object@data@response_matrix), canonical)
    expect_identical(make_run_plan(object, "run_discovery_analysis")@pairs, plan@pairs)
})

test_that("RunPlan structural validity rejects malformed execution values", {
    object <- run_plan_fixture()
    plan <- make_run_plan(object, "run_discovery_analysis")
    for (value in list(integer(), c(1L, 2L), NA_integer_, -1L)) {
        malformed <- plan
        malformed@B1 <- value
        expect_error(validObject(malformed))
    }
    malformed <- plan
    malformed@B1 <- .Machine$integer.max
    malformed@B2 <- 1L
    malformed@B3 <- 0L
    expect_error(validObject(malformed))
    malformed <- plan
    malformed@n_calibration_pairs <- 1L
    expect_error(validObject(malformed))
    malformed <- plan
    malformed@pairs$pass_qc[1L] <- NA
    expect_error(validObject(malformed))
    malformed <- plan
    malformed@analysis <- "discovery_analysis"
    expect_error(validObject(malformed))
    # Empty noncalibration scopes are legitimate numerical work blocks; the
    # complete-family result guard decides whether they can be committed.
    empty <- make_run_plan(object, "run_discovery_analysis", plan@pairs[0, , drop = FALSE])
    expect_true(validObject(empty))
    expect_identical(empty@n_calibration_pairs, integer())
    expect_identical(empty@calibration_group_size, integer())
})

test_that("RunPlan payload preserves the original ordered identity contract", {
    object <- run_plan_fixture()
    plan <- make_run_plan(object, "discovery_analysis")
    expected <- list(analysis = "run_discovery_analysis", pairs = plan@pairs,
        B1 = plan@B1, B2 = plan@B2, B3 = plan@B3,
        n_calibration_pairs = NULL, calibration_group_size = NULL)
    expect_identical(run_plan_payload(plan), expected)
    expect_identical(plan@identity$execution, digest::digest(expected, algo = "sha256"))
    expect_identical(plan@identity, run_plan_identity(object, plan))
})

test_that("RunPlan owns its selected family independently of mutable caller tables", {
    object <- run_plan_fixture()
    pairs <- data.table::as.data.table(object@preparation@discovery_pairs_with_info)
    original <- data.table::copy(pairs$response_id)
    plan <- make_run_plan(object, "run_discovery_analysis", pairs)
    data.table::set(pairs, i = 1L, j = "response_id", value = "changed-after-planning")
    expect_identical(plan@pairs$response_id, original)
    expect_silent(validate_run_plan(object, plan))
})
