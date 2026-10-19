# Shared deterministic inputs and observations for package regression tests.
# Historical old-versus-new helpers live only in tools/refactor/consistency.
refactor_fixture_dir <- function() {
    getOption("sceptre.refactor.fixture_dir", testthat::test_path("fixtures", "refactor"))
}
refactor_private <- function(name) getFromNamespace(name, "sceptre")

refactor_make_inputs <- function() {
    withr::with_seed(20261016, {
        n <- 160L
        guide_ids <- c("a1", "a2", "b1", "b2", "nt1", "nt2", "nt3", "nt4")
        targets <- c("A", "A", "B", "B", rep("non-targeting", 4L))
        guide <- matrix(0, length(guide_ids), n,
                        dimnames = list(guide_ids, paste0("cell_", seq_len(n))))
        assigned <- rep(seq_along(guide_ids), length.out = n)
        guide[cbind(assigned, seq_len(n))] <- 50
        guide[1L, which(assigned == 1L)[seq_len(4L)]] <- 43
        response <- matrix(stats::rpois(8L * n, 3), 8L, n,
                           dimnames = list(paste0("r", seq_len(8L)), colnames(guide)))
        response[1L, assigned %in% 1:2] <- response[1L, assigned %in% 1:2] + 20
        response[2L, assigned %in% 3:4] <- 0
        response[3L, ] <- 0
        response[3L, assigned == 1L] <- 5
        response[3L, assigned >= 5L] <- 3
        response[4L, ] <- 0
        high <- guide
        high[1L, seq(10L, n, 16L)] <- 50
        high[3L, seq(5L, n, 16L)] <- 50
        high[6L, seq(3L, n, 16L)] <- 50
        list(
            response_matrix = response,
            low_grna_matrix = guide,
            high_grna_matrix = high,
            grna_target_data_frame = data.frame(grna_id = guide_ids, grna_target = targets),
            extra_covariates = data.frame(
                batch = factor(rep(c("b", "a", "b", "b", "a"), length.out = n)),
                z = sin(seq_len(n)), row.names = colnames(guide)
            ),
            discovery_pairs = data.frame(
                grna_target = c("B", "A", "A", "B", "A"),
                response_id = c("r5", "r3", "r1", "r2", "r4")
            ),
            positive_control_pairs = data.frame(grna_target = c("A", "B"), response_id = c("r1", "r2")),
            remove_cells = c(2L, 17L, 49L, 103L, 158L)
        )
    })
}

refactor_case_definitions <- function() {
    list(
        low_nt_perm = list(moi = "low", control = "nt_cells", mechanism = "permutations", integration = "union", side = "both", approximation = "skew_normal", output = 3L),
        low_complement_perm = list(moi = "low", control = "complement", mechanism = "permutations", integration = "union", side = "left", approximation = "no_approximation", output = 2L),
        high_complement_crt = list(moi = "high", control = "complement", mechanism = "crt", integration = "union", side = "right", approximation = "skew_normal", output = 3L),
        low_complement_crt = list(moi = "low", control = "complement", mechanism = "crt", integration = "union", side = "both", approximation = "no_approximation", output = 1L),
        high_perm = list(moi = "high", control = "complement", mechanism = "permutations", integration = "union", side = "both", approximation = "no_approximation", output = 3L),
        low_bonferroni = list(moi = "low", control = "nt_cells", mechanism = "permutations", integration = "bonferroni", side = "left", approximation = "skew_normal", output = 2L),
        low_singleton = list(moi = "low", control = "complement", mechanism = "permutations", integration = "singleton", side = "right", approximation = "no_approximation", output = 1L)
    )
}

refactor_expand_inputs <- function(inputs, copies = 8L) {
    n <- ncol(inputs$response_matrix)
    ids <- paste0("cell_", seq_len(n * copies))
    for (name in c("response_matrix", "low_grna_matrix", "high_grna_matrix")) {
        inputs[[name]] <- inputs[[name]][, rep(seq_len(n), copies), drop = FALSE]
        colnames(inputs[[name]]) <- ids
    }
    inputs$extra_covariates <- inputs$extra_covariates[rep(seq_len(n), copies), , drop = FALSE]
    rownames(inputs$extra_covariates) <- ids
    inputs$remove_cells <- as.integer(unlist(lapply(seq_len(copies) - 1L, function(i) inputs$remove_cells + i * n)))
    inputs
}

refactor_import <- function(inputs, case, backed = FALSE, directory = NULL) {
    sceptre::import_data(
        response_matrix = inputs$response_matrix,
        grna_matrix = inputs[[paste0(case$moi, "_grna_matrix")]],
        grna_target_data_frame = inputs$grna_target_data_frame,
        extra_covariates = inputs$extra_covariates, moi = case$moi,
        use_ondisc = backed, directory_to_write = directory
    )
}

refactor_prepare <- function(inputs, case, backed = FALSE, directory = NULL, assignment_threshold = 40) {
    object <- refactor_import(inputs, case, backed, directory)
    object <- sceptre::set_analysis_parameters(
        object, discovery_pairs = inputs$discovery_pairs,
        positive_control_pairs = inputs$positive_control_pairs,
        formula_object = stats::as.formula("~ batch + z", env = baseenv()),
        side = case$side, grna_integration_strategy = case$integration,
        control_group = case$control, resampling_mechanism = case$mechanism,
        resampling_approximation = case$approximation
    )
    object <- sceptre::assign_grnas(object, method = "thresholding", threshold = assignment_threshold, print_progress = FALSE)
    sceptre::run_qc(object,
        response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
        n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L,
        additional_cells_to_remove = inputs$remove_cells
    )
}

# Keep table types and row order, but omit data.table's internal self-reference
# pointer, which is not useful when comparing saved values.
refactor_plain <- function(x) {
    if (is.data.frame(x)) {
        return(list(class = class(x), names = names(x), row_names = rownames(x),
                    columns = lapply(x, refactor_plain)))
    }
    if (is.list(x)) return(lapply(x, refactor_plain))
    x
}

refactor_observe <- function(object) {
    prep <- object@preparation
    observed <- list(analysis_cells = prep@analysis_cells,
        analysis_groups = prep@analysis_groups,
        cell_removal_metrics = prep@cell_removal_metrics,
        nt_nonzero_matrix = prep@nt_nonzero_matrix, control_nonzero_totals = prep@control_nonzero_totals,
        discovery_pairs_with_info = prep@discovery_pairs_with_info,
        positive_control_pairs_with_info = prep@positive_control_pairs_with_info,
        side_code = object@spec@side_code,
        workflow_status = refactor_private("workflow_status")(object),
        response_precomputations = object@fits@response_fits)
    for (analysis in c("calibration", "power", "discovery")) {
        method <- switch(analysis, calibration = "run_calibration_check",
            power = "run_power_check", discovery = "run_discovery_analysis")
        table <- object@results@records[[method]]$table
        if (is.null(table)) table <- data.frame()
        observed[[paste0(analysis, "_result")]] <- table
    }
    lapply(observed, refactor_plain)
}

refactor_capture <- function(code) {
    warnings <- character()
    value <- withCallingHandlers(suppressMessages(code), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    list(value = value, warnings = warnings)
}

refactor_run_case <- function(inputs, case, parallel = FALSE, backed = FALSE, directory = NULL) {
    withr::with_seed(271828, refactor_capture({
        object <- refactor_prepare(inputs, case, backed, directory)
        prepared <- refactor_observe(object)
        object <- sceptre::run_calibration_check(object, n_calibration_pairs = 4L,
            calibration_group_size = 2L, print_progress = FALSE,
            parallel = parallel, n_processors = 2L, output_amount = case$output)
        calibrated <- refactor_observe(object)
        object <- sceptre::run_power_check(object, print_progress = FALSE,
            parallel = parallel, n_processors = 2L, output_amount = case$output)
        object <- sceptre::run_discovery_analysis(object, print_progress = FALSE,
            parallel = parallel, n_processors = 2L, output_amount = case$output)
        list(prepared = prepared, calibrated = calibrated, final = refactor_observe(object),
             rng_after = get(".Random.seed", envir = .GlobalEnv))
    }))
}

# Calculate eligibility directly from the toy count matrices, without calling
# SCEPTRE's assignment, indexing or QC helpers. These inputs use threshold 40,
# remove five named cells and require seven nonzero counts in each group.
workflow_expected_pairs <- function(inputs, case, pairs) {
    guides <- inputs[[paste0(case$moi, "_grna_matrix")]]
    metadata <- inputs$grna_target_data_frame
    cells <- setdiff(seq_len(ncol(guides)), inputs$remove_cells)
    nt <- which(colSums(guides[metadata$grna_target == "non-targeting", , drop = FALSE] >= 40) > 0)
    records <- lapply(seq_len(nrow(pairs)), function(i) {
        ids <- which(metadata$grna_target == pairs$grna_target[i])
        groups <- if (case$integration == "union") list(ids) else as.list(ids)
        passes <- vapply(groups, function(rows) {
            treatment <- intersect(cells, which(colSums(guides[rows, , drop = FALSE] >= 40) > 0))
            control <- if (case$control == "nt_cells") intersect(cells, nt) else setdiff(cells, treatment)
            y <- inputs$response_matrix[pairs$response_id[i], ]
            sum(y[treatment] > 0) >= 7L && sum(y[control] > 0) >= 7L
        }, logical(1))
        if (case$integration == "bonferroni") passes <- any(passes)
        paste(pairs$response_id[i], pairs$grna_target[i], passes, sep = "/")
    })
    sort(unlist(records, use.names = FALSE))
}

expect_workflow_results <- function(actual, inputs, case) {
    expect_identical(actual$warnings, character(), info = "unexpected workflow warnings")
    for (table in c("calibration_result", "power_result", "discovery_result")) {
        result <- actual$value$final[[table]]$columns
        expect_true(all(c("response_id", "grna_target", "pass_qc", "p_value") %in% names(result)))
        expect_false(anyNA(result$pass_qc))
        expect_true(all(result$response_id %in% rownames(inputs$response_matrix)))
        expect_true(all(is.finite(result$p_value[result$pass_qc])))
        expect_true(all(result$p_value[result$pass_qc] >= 0 & result$p_value[result$pass_qc] <= 1),
                    info = "eligible p-values on these toy inputs must lie in [0, 1]")
        expect_true(all(is.na(result$p_value[!result$pass_qc])))
        if (table == "calibration_result") {
            expect_length(result$response_id, 4L)
            expect_true(all(result$pass_qc))
        } else {
            pairs <- inputs[[if (table == "power_result") "positive_control_pairs" else "discovery_pairs"]]
            observed <- sort(paste(result$response_id, result$grna_target, result$pass_qc, sep = "/"))
            expect_identical(observed, workflow_expected_pairs(inputs, case, pairs),
                             info = paste(table, "pair identities and independent QC"))
            # The input generator adds 20 counts to A-targeted cells for r1.
            # This is a strong positive effect, regardless of inference method.
            positive <- result$response_id == "r1" & result$grna_target == "A"
            expect_true(any(positive & result$pass_qc))
            if (case$side == "left") {
                expect_true(all(result$p_value[positive] > .9))
            } else {
                expect_true(all(result$p_value[positive] < .05))
            }
            if ("fold_change" %in% names(result)) expect_true(all(result$fold_change[positive] > 2))
        }
    }
}

refactor_cache_observations <- function(inputs, case) {
    refactor_capture({
        object <- refactor_prepare(inputs, case)
        cold <- sceptre::run_discovery_analysis(object, print_progress = FALSE, output_amount = 2L)
        warm <- sceptre::run_discovery_analysis(cold, print_progress = FALSE, output_amount = 2L)
        cleared <- cold
        cleared@fits@response_fits <- list()
        cleared <- sceptre::run_discovery_analysis(cleared, print_progress = FALSE, output_amount = 2L)
        changed <- sceptre::set_analysis_parameters(cold,
            discovery_pairs = inputs$discovery_pairs[rev(seq_len(nrow(inputs$discovery_pairs))), ],
            positive_control_pairs = inputs$positive_control_pairs,
            formula_object = stats::as.formula("~ batch + z", env = baseenv()),
            control_group = case$control, resampling_mechanism = case$mechanism,
            resampling_approximation = case$approximation, side = "right")
        side_cache <- changed@fits@response_fits
        changed <- sceptre::run_discovery_analysis(changed, print_progress = FALSE, output_amount = 2L)
        formula_changed <- sceptre::set_analysis_parameters(cold,
            discovery_pairs = inputs$discovery_pairs,
            formula_object = stats::as.formula("~ batch", env = baseenv()),
            control_group = case$control, resampling_mechanism = case$mechanism)
        reassigned <- sceptre::assign_grnas(cold, method = "thresholding", threshold = 45, print_progress = FALSE)
        reassigned_cache <- reassigned@fits@response_fits
        reassigned <- sceptre::run_qc(reassigned,
            response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
            n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L,
            additional_cells_to_remove = inputs$remove_cells)
        reassigned <- sceptre::run_discovery_analysis(reassigned, print_progress = FALSE, output_amount = 2L)
        reassigned_cold <- refactor_prepare(inputs, case, assignment_threshold = 45)
        reassigned_cold <- sceptre::run_discovery_analysis(reassigned_cold, print_progress = FALSE, output_amount = 2L)
        qc_changed <- sceptre::run_qc(cold,
            response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
            n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L,
            additional_cells_to_remove = c(inputs$remove_cells, 25L))
        list(cold = refactor_observe(cold), warm = refactor_observe(warm),
             cleared = refactor_observe(cleared), changed_side = refactor_observe(changed),
             retained_side_cache = side_cache,
             formula_cache = formula_changed@fits@response_fits,
             reassigned_cache = reassigned_cache,
             reassigned = refactor_observe(reassigned), reassigned_cold = refactor_observe(reassigned_cold),
             changed_qc_cache = qc_changed@fits@response_fits)
    })
}

refactor_index_observations <- function() {
    object <- methods::new("sceptre_object")
    object <- refactor_private("set_response_matrix")(object,
        set_matrix_accessibility(matrix(101:112, nrow = 1L)))
    object@data@low_moi <- TRUE
    object@spec@control_group_complement <- FALSE
    object@preparation@analysis_cells <- c(2L, 4L, 5L, 8L, 11L)
    raw <- list(grna_group_idxs = list(A = c(5L, 4L)),
                indiv_nt_grna_idxs = list(ntA = 8L, ntB = c(2L, 11L)))
    # Use the saved assignment format; build the combined NT-cell list below.
    object@preparation@source_groups <- raw
    before <- serialize(raw, NULL)
    out <- refactor_private("update_grna_assignments_given_qc")(object)
    again <- refactor_private("update_grna_assignments_given_qc")(object)
    groups <- out@preparation@analysis_groups
    combined <- c(groups$grna_group_idxs$A, groups$all_nt_idxs)
    list(assignments = groups, repeated = again@preparation@analysis_groups,
         caller_unchanged = identical(serialize(raw, NULL), before),
         treatment_source = object@preparation@analysis_cells[groups$grna_group_idxs$A],
         nt_source = object@preparation@analysis_cells[groups$all_nt_idxs],
         combined_source = object@preparation@analysis_cells[combined],
         extracted_response = (101:112)[object@preparation@analysis_cells][combined],
         extracted_covariates = cbind(intercept = 1, cell = 101:112)[object@preparation@analysis_cells, , drop = FALSE][combined, , drop = FALSE])
}

refactor_edge_observations <- function(inputs, case) {
    outcome <- function(code) {
        tryCatch(refactor_capture(code), error = function(e) list(error_class = class(e), error = conditionMessage(e)))
    }
    duplicate <- inputs
    duplicate$discovery_pairs <- rbind(inputs$discovery_pairs, inputs$discovery_pairs[1L, ])
    no_pairs <- inputs
    no_pairs$discovery_pairs <- data.frame(grna_target = "A", response_id = "r4")
    no_pairs$positive_control_pairs <- inputs$positive_control_pairs[FALSE, ]
    list(
        duplicate_requests = outcome(refactor_prepare(duplicate, case)),
        no_passing_pairs = outcome({
            object <- refactor_prepare(no_pairs, case)
            refactor_observe(sceptre::run_discovery_analysis(object, print_progress = FALSE))
        }),
        intercept_only_nt = outcome({
            object <- refactor_prepare(inputs, readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_nt_perm)
            object <- sceptre::set_analysis_parameters(object, discovery_pairs = inputs$discovery_pairs,
                formula_object = stats::as.formula("~ 1", env = baseenv()), control_group = "nt_cells",
                resampling_mechanism = "permutations")
            refactor_observe(sceptre::run_discovery_analysis(object, print_progress = FALSE))
        })
    )
}
