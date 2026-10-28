# Frozen behavior suite: shared harness.
#
# Defines the inputs, the cases and how each case is observed. The cases call
# only exported functions. Observations have two parts:
#
# - public: what a user can see (exported getters, get_result(), printed
#   summaries, written outputs, errors, warnings and the RNG state);
# - internal: `sceptre_object` slots, stored under elements named `internal`.
#
# tests/helper-frozen.R snapshots the two parts separately, so a restructuring
# that moves slots can waive the internal snapshots explicitly. Nothing here
# writes outside `tempdir()`.

consistency_input_seed <- 20261016L
consistency_case_seed <- 271828L

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

# Synthetic 160-cell, 8-response, 8-guide data. The guide design, covariates,
# pairs and removed cells follow the refactor characterization inputs (#234).
# Responses are negative binomial (size 2) rather than Poisson, so that every
# response has a well-determined dispersion: with Poisson counts the estimate
# is near the 1000 cap and moves with floating-point round-off, which made
# results differ between platforms.
consistency_synthetic_inputs <- function() {
    withr::with_seed(consistency_input_seed, {
        n <- 160L
        counts <- function(k, mu) stats::rnbinom(k, size = 2, mu = mu)
        guide_ids <- c("a1", "a2", "b1", "b2", "nt1", "nt2", "nt3", "nt4")
        targets <- c("A", "A", "B", "B", rep("non-targeting", 4L))
        guide <- matrix(
            0,
            length(guide_ids),
            n,
            dimnames = list(guide_ids, paste0("cell_", seq_len(n)))
        )
        assigned <- rep(seq_along(guide_ids), length.out = n)
        guide[cbind(assigned, seq_len(n))] <- 50
        guide[1L, which(assigned == 1L)[seq_len(4L)]] <- 43
        response <- matrix(
            counts(8L * n, 3),
            8L,
            n,
            dimnames = list(paste0("r", seq_len(8L)), colnames(guide))
        )
        # r1: strong effect of target A; r2: silenced by target B;
        # r3: expressed only in a1 and NT cells; r4: all zero (fails QC).
        response[1L, assigned %in% 1:2] <- response[1L, assigned %in% 1:2] +
            20
        response[2L, assigned %in% 3:4] <- 0
        response[3L, ] <- 0
        response[3L, assigned == 1L] <- counts(sum(assigned == 1L), 5)
        response[3L, assigned >= 5L] <- counts(sum(assigned >= 5L), 3)
        response[4L, ] <- 0
        high <- guide
        high[1L, seq(10L, n, 16L)] <- 50
        high[3L, seq(5L, n, 16L)] <- 50
        high[6L, seq(3L, n, 16L)] <- 50
        list(
            response_matrix = response,
            low_grna_matrix = guide,
            high_grna_matrix = high,
            grna_target_data_frame = data.frame(
                grna_id = guide_ids,
                grna_target = targets
            ),
            extra_covariates = data.frame(
                batch = factor(rep(c("b", "a", "b", "b", "a"), length.out = n)),
                # Rounded so that the input is bit-identical on every
                # platform (libm results can differ in the last bit).
                z = round(sin(seq_len(n)), 6L),
                row.names = colnames(guide)
            ),
            discovery_pairs = data.frame(
                grna_target = c("B", "A", "A", "B", "A"),
                response_id = c("r5", "r3", "r1", "r2", "r4")
            ),
            positive_control_pairs = data.frame(
                grna_target = c("A", "B"),
                response_id = c("r1", "r2")
            ),
            remove_cells = c(2L, 17L, 49L, 103L, 158L)
        )
    })
}

# The package's own example data sets, as loaded from the package under test.
# Their digests are recorded so that a changed input is reported as such.
consistency_example_data <- function() {
    env <- new.env(parent = emptyenv())
    utils::data(
        list = c(
            "lowmoi_example_data",
            "highmoi_example_data",
            "grna_target_data_frame_highmoi"
        ),
        package = "sceptre",
        envir = env
    )
    low <- env$lowmoi_example_data
    low$grna_target_data_frame <- as.data.frame(low$grna_target_data_frame)
    list(
        low = low,
        high = env$highmoi_example_data,
        high_targets = env$grna_target_data_frame_highmoi
    )
}

# MD5 of an object, from serialization format 2 without its 14-byte header
# (which records the writing R version), so it depends on the object only.
consistency_md5 <- function(x) {
    file <- tempfile("sceptre-consistency-")
    on.exit(unlink(file), add = TRUE)
    writeBin(serialize(x, NULL, xdr = TRUE, version = 2L)[-seq_len(14L)], file)
    unname(tools::md5sum(file))
}

consistency_extdata <- function(...) {
    path <- system.file("extdata", ..., package = "sceptre")
    if (!nzchar(path)) {
        stop(
            "Example file not found in the installed package: ",
            file.path(...)
        )
    }
    path
}

# Every n-th row, in order; deterministic subsets that consume no RNG draws.
consistency_every_nth <- function(df, n) {
    df <- as.data.frame(df)
    rownames(df) <- NULL
    df[seq(1L, nrow(df), by = n), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Observation
# ---------------------------------------------------------------------------

# A complete, serializable form of a value. A data frame (or data.table)
# becomes list(class, names, row_names, columns), keeping its classes, row
# names and column types (factors with their levels) but not data.table's
# self-reference pointer. Formulas become text.
consistency_plain <- function(x) {
    if (is.data.frame(x)) {
        return(list(
            class = class(x),
            names = names(x),
            row_names = attr(x, "row.names"),
            columns = lapply(unclass(x), consistency_plain)
        ))
    }
    if (inherits(x, "formula")) {
        return(paste(deparse(x, width.cutoff = 500L), collapse = " "))
    }
    if (is.list(x)) {
        return(lapply(x, consistency_plain))
    }
    x
}

consistency_matrix_table <- function(m) {
    m <- methods::as(m, "TsparseMatrix")
    out <- data.frame(
        row = as.integer(m@i + 1L),
        col = as.integer(m@j + 1L),
        value = as.numeric(m@x)
    )
    out <- out[order(out$col, out$row), , drop = FALSE]
    rownames(out) <- NULL
    out
}

consistency_matrix_info <- function(m) {
    list(class = class(m), dim = dim(m), dimnames = dimnames(m))
}

# Imported data: public getters, plus the slots that record the import.
consistency_observe_import <- function(object) {
    response <- sceptre::get_response_matrix(object)
    grna <- sceptre::get_grna_matrix(object)
    list(
        response_matrix = consistency_matrix_info(response),
        grna_matrix = consistency_matrix_info(grna),
        tables = list(
            response_matrix = consistency_plain(consistency_matrix_table(
                response
            )),
            grna_matrix = consistency_plain(consistency_matrix_table(grna)),
            cell_covariates = consistency_plain(sceptre::get_cell_covariates(
                object
            ))
        ),
        internal = consistency_internal_part(list(
            covariate_names = object@covariate_names,
            low_moi = object@low_moi,
            tables = list(
                grna_target_data_frame = consistency_plain(
                    object@grna_target_data_frame
                )
            )
        ))
    )
}

# Assignments from the exported getter, one row per (guide, cell), in guide
# order and then cell order.
consistency_assignment_table <- function(object) {
    m <- methods::as(sceptre::get_grna_assignments(object), "TsparseMatrix")
    order <- order(m@i, m@j)
    consistency_plain(data.frame(
        grna_id = rownames(m)[m@i[order] + 1L],
        cell = as.integer(m@j[order] + 1L)
    ))
}

consistency_internal_slots <- c(
    "cells_in_use",
    "cell_removal_metrics",
    "grna_assignments",
    "grnas_per_cell",
    "cells_w_zero_or_twoplus_grnas",
    "M_matrix",
    "n_nonzero_tot_vector",
    "negative_control_pairs",
    "n_ok_discovery_pairs",
    "n_ok_positive_control_pairs",
    "calibration_group_size",
    "n_calibration_pairs",
    "B1",
    "B2",
    "B3",
    "side_code",
    "formula_object",
    "control_group_complement",
    "run_permutations",
    "grna_integration_strategy",
    "resampling_approximation",
    "multiple_testing_method",
    "multiple_testing_alpha",
    "grna_assignment_method",
    "grna_assignment_hyperparameters",
    "cellwise_qc_thresholds",
    "mean_cells_per_grna",
    "functs_called",
    "last_function_called",
    "response_precomputations"
)

consistency_analyses <- c(
    calibration_result = "run_calibration_check",
    power_result = "run_power_check",
    discovery_result = "run_discovery_analysis"
)

# The value of `x`, or NULL when it calls a getter for a step that has not run
# yet. Asking the exported getters, rather than reading which functions were
# called from a slot, keeps the public outputs free of internal state.
consistency_if_run <- function(x) {
    tryCatch(x, error = function(e) {
        if (grepl("has not yet been (called|run)", conditionMessage(e))) {
            return(NULL)
        }
        stop(e)
    })
}

# The internal part of an observation, or NULL when the internal snapshots are
# waived. `x` is evaluated only when it is used, so a restructuring that removes
# slots can still check the public outputs.
consistency_internal_part <- function(x) {
    if (identical(Sys.getenv("SCEPTRE_FROZEN_WAIVE_INTERNAL"), "true")) {
        return(NULL)
    }
    x
}

# An analysis object: public tables (getters and get_result()) and internal
# slots, including the requested pairs with their pairwise QC.
consistency_observe <- function(object) {
    tables <- list(
        cell_covariates = consistency_plain(sceptre::get_cell_covariates(
            object
        ))
    )
    # assigning NULL drops the element, so steps that have not run are absent
    tables$grna_assignments <- consistency_if_run(consistency_assignment_table(
        object
    ))
    for (name in names(consistency_analyses)) {
        tables[[name]] <- consistency_if_run(consistency_plain(
            sceptre::get_result(object, consistency_analyses[[name]])
        ))
    }
    internal <- consistency_internal_part({
        internal <- lapply(consistency_internal_slots, function(name) {
            consistency_plain(methods::slot(object, name))
        })
        names(internal) <- consistency_internal_slots
        internal$tables <- list(
            discovery_pairs_with_info = consistency_plain(
                object@discovery_pairs_with_info
            ),
            positive_control_pairs_with_info = consistency_plain(
                object@positive_control_pairs_with_info
            )
        )
        internal
    })
    list(tables = tables, internal = internal)
}

consistency_summary_text <- function(object) {
    withr::with_options(list(crayon.enabled = FALSE, width = 80L), {
        utils::capture.output(print(object))
    })
}

# Evaluate `code`, muffling messages and recording warnings and errors.
consistency_capture <- function(code) {
    warnings <- character()
    value <- tryCatch(
        withCallingHandlers(
            suppressMessages(code),
            warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) {
            structure(
                list(message = conditionMessage(e)),
                class = "consistency_error"
            )
        }
    )
    error <- if (inherits(value, "consistency_error")) value$message else NULL
    if (!is.null(error)) {
        value <- NULL
    }
    list(value = value, warnings = warnings, error = error)
}

# ---------------------------------------------------------------------------
# Synthetic workflows (seven settings from #234)
# ---------------------------------------------------------------------------

consistency_workflow_settings <- function() {
    list(
        low_nt_perm = list(
            moi = "low",
            control = "nt_cells",
            mechanism = "permutations",
            integration = "union",
            side = "both",
            approximation = "skew_normal",
            output = 3L
        ),
        low_complement_perm = list(
            moi = "low",
            control = "complement",
            mechanism = "permutations",
            integration = "union",
            side = "left",
            approximation = "no_approximation",
            output = 2L
        ),
        high_complement_crt = list(
            moi = "high",
            control = "complement",
            mechanism = "crt",
            integration = "union",
            side = "right",
            approximation = "skew_normal",
            output = 3L
        ),
        low_complement_crt = list(
            moi = "low",
            control = "complement",
            mechanism = "crt",
            integration = "union",
            side = "both",
            approximation = "no_approximation",
            output = 1L
        ),
        high_perm = list(
            moi = "high",
            control = "complement",
            mechanism = "permutations",
            integration = "union",
            side = "both",
            approximation = "no_approximation",
            output = 3L
        ),
        low_bonferroni = list(
            moi = "low",
            control = "nt_cells",
            mechanism = "permutations",
            integration = "bonferroni",
            side = "left",
            approximation = "skew_normal",
            output = 2L
        ),
        low_singleton = list(
            moi = "low",
            control = "complement",
            mechanism = "permutations",
            integration = "singleton",
            side = "right",
            approximation = "no_approximation",
            output = 1L
        )
    )
}

consistency_import_synthetic <- function(inputs, setting) {
    sceptre::import_data(
        response_matrix = inputs$response_matrix,
        grna_matrix = inputs[[paste0(setting$moi, "_grna_matrix")]],
        grna_target_data_frame = inputs$grna_target_data_frame,
        extra_covariates = inputs$extra_covariates,
        moi = setting$moi
    )
}

consistency_qc_synthetic <- function(object, inputs, extra_removed = NULL) {
    sceptre::run_qc(
        object,
        response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1),
        n_nonzero_trt_thresh = 7L,
        n_nonzero_cntrl_thresh = 7L,
        additional_cells_to_remove = c(inputs$remove_cells, extra_removed)
    )
}

consistency_prepare_synthetic <- function(
    inputs,
    setting,
    threshold = 40,
    formula = stats::as.formula("~ batch + z", env = baseenv())
) {
    object <- consistency_import_synthetic(inputs, setting)
    object <- sceptre::set_analysis_parameters(
        object,
        discovery_pairs = inputs$discovery_pairs,
        positive_control_pairs = inputs$positive_control_pairs,
        formula_object = formula,
        side = setting$side,
        grna_integration_strategy = setting$integration,
        control_group = setting$control,
        resampling_mechanism = setting$mechanism,
        resampling_approximation = setting$approximation
    )
    object <- sceptre::assign_grnas(
        object,
        method = "thresholding",
        threshold = threshold,
        print_progress = FALSE
    )
    consistency_qc_synthetic(object, inputs)
}

consistency_run_workflow <- function(setting, parallel = FALSE) {
    inputs <- consistency_synthetic_inputs()
    consistency_capture({
        object <- consistency_prepare_synthetic(inputs, setting)
        prepared <- consistency_internal_part(object@response_precomputations)
        object <- sceptre::run_calibration_check(
            object,
            n_calibration_pairs = 4L,
            calibration_group_size = 2L,
            print_progress = FALSE,
            parallel = parallel,
            n_processors = 2L,
            output_amount = setting$output
        )
        object <- sceptre::run_power_check(
            object,
            print_progress = FALSE,
            parallel = parallel,
            n_processors = 2L,
            output_amount = setting$output
        )
        object <- sceptre::run_discovery_analysis(
            object,
            print_progress = FALSE,
            parallel = parallel,
            n_processors = 2L,
            output_amount = setting$output
        )
        out <- consistency_observe(object)
        if (!is.null(prepared)) {
            out$internal$prepared_precomputations <- prepared
        }
        out$summary <- consistency_summary_text(object)
        out
    })
}

# Order cached response fits by name, so that two-worker runs (whose caches
# can be filled in a different order) compare equal to serial runs.
consistency_sort_precomputations <- function(observation) {
    internal <- observation$value$internal
    for (name in c("response_precomputations", "prepared_precomputations")) {
        fits <- internal[[name]]
        if (length(fits)) internal[[name]] <- fits[sort(names(fits))]
    }
    observation$value$internal <- internal
    observation
}

# ---------------------------------------------------------------------------
# Import
# ---------------------------------------------------------------------------

consistency_case_import_matrix_low <- function() {
    inputs <- consistency_synthetic_inputs()
    consistency_capture({
        consistency_observe_import(consistency_import_synthetic(
            inputs,
            list(moi = "low")
        ))
    })
}

# Guide counts given per guide and summed to vectors (#234 aggregation input).
consistency_case_import_matrix_vector <- function() {
    counts <- methods::new(
        "dgTMatrix",
        i = c(0L, 0L, 1L, 2L, 3L, 1L, 3L),
        j = c(0L, 0L, 1L, 2L, 0L, 0L, 2L),
        x = c(.25, .75, 0, 1.5, 2, 3, .5),
        Dim = c(4L, 3L),
        Dimnames = list(paste0("g", 1:4), c("c-z", "c-a", "c-m"))
    )
    metadata <- data.frame(
        grna_id = c("g3", "g1", "g1", "g2", "g4"),
        grna_target = c("targetZ", "targetA", "targetA", "targetA", "targetZ"),
        vector_id = c("vZ", "vA", "vA", "vA", "vZ")
    )
    response <- matrix(
        c(1, 2, 3, 4, 5, 6),
        2L,
        dimnames = list(c("r-z", "r-a"), colnames(counts))
    )
    consistency_capture({
        consistency_observe_import(sceptre::import_data(
            response,
            counts,
            metadata,
            moi = "high"
        ))
    })
}

# Known bug (A1 #4): with matching dimnames, Matrix stores a symmetric dense
# input as one triangle and the lower-triangle value is lost (#234 example).
consistency_case_import_dense_symmetric <- function() {
    guides <- matrix(
        c(1, 0, 0, 1),
        2L,
        2L,
        dimnames = list(c("g1", "g2"), c("c1", "c2"))
    )
    targets <- data.frame(
        grna_id = c("g1", "g2"),
        grna_target = c("non-targeting", "target")
    )
    import <- function(lower) {
        response <- matrix(
            c(1, lower, 2, 3),
            2L,
            2L,
            dimnames = list(c("c1", "c2"), c("c1", "c2"))
        )
        consistency_observe_import(sceptre::import_data(
            response,
            guides,
            targets,
            "low",
            response_names = c("MT-first", "second")
        ))
    }
    consistency_capture(list(
        symmetric = import(2),
        near_symmetric = import(2 + 4 * .Machine$double.eps)
    ))
}

consistency_import_cellranger <- function(directories, ...) {
    sceptre::import_data_from_cellranger(
        directories = directories,
        moi = "high",
        grna_target_data_frame = consistency_example_data()$high_targets,
        ...
    )
}

consistency_case_import_cellranger_two_directories <- function() {
    consistency_capture({
        consistency_observe_import(consistency_import_cellranger(c(
            consistency_extdata("highmoi_example", "gem_group_1"),
            consistency_extdata("highmoi_example", "gem_group_2")
        )))
    })
}

consistency_case_import_cellranger_one_directory <- function() {
    consistency_capture({
        consistency_observe_import(consistency_import_cellranger(
            consistency_extdata("highmoi_example", "gem_group_1")
        ))
    })
}

# Known bug (a7b08441, A2): one directory drops the user's `batch` column.
# Only the covariates are observed; the counts are covered by the case above.
consistency_case_import_cellranger_one_directory_covariates <- function() {
    directory <- consistency_extdata("highmoi_example", "gem_group_1")
    n_cells <- length(readLines(file.path(directory, "barcodes.tsv.gz")))
    consistency_capture({
        object <- consistency_import_cellranger(
            directory,
            extra_covariates = data.frame(
                batch = factor(rep(c("x", "y"), length.out = n_cells)),
                lane = rep(c(1, 2, 3), length.out = n_cells)
            )
        )
        list(
            tables = list(
                cell_covariates = consistency_plain(
                    sceptre::get_cell_covariates(object)
                )
            ),
            internal = consistency_internal_part(list(
                covariate_names = object@covariate_names
            ))
        )
    })
}

consistency_case_import_parse <- function() {
    directory <- consistency_extdata("parse_example")
    targets <- data.frame(
        grna_id = c("guide_A", "guide_B", "guide_C"),
        grna_target = c("target-A", "target-B", "non-targeting")
    )
    consistency_capture({
        consistency_observe_import(sceptre::import_data_from_parse(
            gene_mat_fp = file.path(directory, "gene_mat.mtx"),
            grna_mat_fp = file.path(directory, "grna_mat.mtx"),
            all_genes_fp = file.path(directory, "all_genes.csv"),
            all_grnas_fp = file.path(directory, "all_grnas.csv"),
            moi = "low",
            grna_target_data_frame = targets
        ))
    })
}

# ---------------------------------------------------------------------------
# Guide assignment on the package example data
# ---------------------------------------------------------------------------

consistency_import_example <- function(moi) {
    data <- consistency_example_data()
    if (moi == "low") {
        sceptre::import_data(
            response_matrix = data$low$response_matrix,
            grna_matrix = data$low$grna_matrix,
            grna_target_data_frame = data$low$grna_target_data_frame,
            extra_covariates = data$low$extra_covariates,
            moi = "low"
        )
    } else {
        sceptre::import_data(
            response_matrix = data$high$response_matrix,
            grna_matrix = data$high$grna_matrix,
            grna_target_data_frame = data$high_targets,
            extra_covariates = data$high$extra_covariates,
            moi = "high",
            response_names = data$high$gene_names
        )
    }
}

consistency_observe_assignment <- function(object) {
    list(
        tables = list(grna_assignments = consistency_assignment_table(object)),
        internal = consistency_internal_part(list(
            method = object@grna_assignment_method,
            hyperparameters = consistency_plain(
                object@grna_assignment_hyperparameters
            ),
            grnas_per_cell = object@grnas_per_cell,
            cells_w_zero_or_twoplus_grnas = object@cells_w_zero_or_twoplus_grnas,
            import_grna_assignment_info = consistency_plain(
                object@import_grna_assignment_info
            ),
            mean_cells_per_grna = object@mean_cells_per_grna
        ))
    )
}

consistency_assignment_case <- function(moi, ...) {
    function() {
        consistency_capture({
            object <- consistency_import_example(moi) |>
                sceptre::set_analysis_parameters()
            object <- sceptre::assign_grnas(
                object,
                print_progress = FALSE,
                ...
            )
            consistency_observe_assignment(object)
        })
    }
}

# ---------------------------------------------------------------------------
# End-to-end workflows on the package example data, with package defaults
# ---------------------------------------------------------------------------

consistency_finish_example <- function(object, write_outputs = FALSE) {
    object <- sceptre::run_qc(object)
    object <- sceptre::run_calibration_check(object, print_progress = FALSE)
    object <- sceptre::run_power_check(object, print_progress = FALSE)
    object <- sceptre::run_discovery_analysis(object, print_progress = FALSE)
    out <- consistency_observe(object)
    out$summary <- consistency_summary_text(object)
    if (write_outputs) {
        directory <- tempfile("sceptre-consistency-outputs-")
        on.exit(unlink(directory, recursive = TRUE), add = TRUE)
        sceptre::write_outputs_to_directory(object, directory)
        out$written_files <- sort(list.files(directory))
        out$written_summary <- readLines(
            file.path(directory, "analysis_summary.txt"),
            warn = FALSE
        )
        out$written_results <- lapply(
            c(
                "results_run_calibration_check.rds",
                "results_run_power_check.rds",
                "results_run_discovery_analysis.rds"
            ),
            function(file) {
                consistency_plain(readRDS(file.path(directory, file)))
            }
        )
    }
    out
}

# Low MOI: trans pairs (every 50th), default assignment, QC and analyses.
consistency_case_example_low_defaults <- function() {
    consistency_capture({
        object <- consistency_import_example("low")
        positive <- sceptre::construct_positive_control_pairs(object)
        discovery <- sceptre::construct_trans_pairs(
            object,
            positive_control_pairs = data.frame(positive),
            pairs_to_exclude = "pc_pairs"
        ) |>
            consistency_every_nth(50L)
        object <- sceptre::set_analysis_parameters(
            object,
            discovery_pairs = discovery,
            positive_control_pairs = positive
        )
        object <- sceptre::assign_grnas(object, print_progress = FALSE)
        consistency_finish_example(object, write_outputs = TRUE)
    })
}

# High MOI: cis pairs (every 8th), mixture assignment, CRT, default QC.
consistency_case_example_high_defaults <- function() {
    consistency_capture({
        object <- consistency_import_example("high")
        positive <- sceptre::construct_positive_control_pairs(object)
        discovery <- sceptre::construct_cis_pairs(
            object,
            positive_control_pairs = positive,
            distance_threshold = 5e6
        ) |>
            consistency_every_nth(8L)
        object <- sceptre::set_analysis_parameters(
            object,
            discovery_pairs = discovery,
            positive_control_pairs = positive,
            side = "left"
        )
        object <- sceptre::assign_grnas(object, print_progress = FALSE)
        consistency_finish_example(object)
    })
}

# ---------------------------------------------------------------------------
# Pair constructors
# ---------------------------------------------------------------------------

consistency_case_pairs_constructors <- function() {
    consistency_capture({
        object <- consistency_import_example("high")
        positive <- sceptre::construct_positive_control_pairs(object)
        # Fresh copies: construct_trans_pairs() may modify its argument (#8,
        # recorded separately in `pairs_trans_caller_mutation`).
        copy <- function() data.frame(positive)
        list(
            positive_classes = class(positive),
            tables = list(
                positive_control = consistency_plain(positive),
                cis_default = consistency_plain(sceptre::construct_cis_pairs(
                    object,
                    positive_control_pairs = copy()
                )),
                cis_5mb = consistency_plain(sceptre::construct_cis_pairs(
                    object,
                    positive_control_pairs = copy(),
                    distance_threshold = 5e6
                )),
                trans_none = consistency_plain(sceptre::construct_trans_pairs(
                    object
                )),
                trans_pc_pairs = consistency_plain(
                    sceptre::construct_trans_pairs(
                        object,
                        positive_control_pairs = copy(),
                        pairs_to_exclude = "pc_pairs"
                    )
                ),
                trans_pc_targets = consistency_plain(
                    sceptre::construct_trans_pairs(
                        object,
                        positive_control_pairs = copy(),
                        pairs_to_exclude = "pairs_containing_pc_targets"
                    )
                )
            )
        )
    })
}

# Known bug (A2 #8): the caller's data.frame becomes a data.table.
consistency_case_pairs_trans_caller_mutation <- function() {
    consistency_capture({
        object <- consistency_import_example("high")
        caller <- data.frame(sceptre::construct_positive_control_pairs(object))
        before <- class(caller)
        invisible(sceptre::construct_trans_pairs(
            object,
            positive_control_pairs = caller,
            pairs_to_exclude = "pc_pairs"
        ))
        list(
            caller_classes_before = before,
            caller_classes_after = class(caller)
        )
    })
}

# Known bug (A2 #9): at midpoint 100 and distance 20, positions 90, 1000, 110
# should select the first and third responses.
consistency_case_pairs_cis_unsorted_positions <- function() {
    response <- matrix(
        1,
        3L,
        4L,
        dimnames = list(c("near", "far", "also_near"), paste0("c", 1:4))
    )
    grna <- matrix(
        c(9, 0, 0, 9, 9, 0, 0, 9),
        2L,
        dimnames = list(c("g1", "nt1"), paste0("c", 1:4))
    )
    targets <- data.frame(
        grna_id = c("g1", "nt1"),
        grna_target = c("t1", "non-targeting"),
        chr = c("chr1", NA),
        start = c(95L, NA),
        end = c(105L, NA)
    )
    positions <- data.frame(
        response_id = c("near", "far", "also_near"),
        chr = "chr1",
        position = c(90L, 1000L, 110L)
    )
    consistency_capture({
        object <- sceptre::import_data(response, grna, targets, moi = "high")
        list(
            tables = list(
                cis_pairs = consistency_plain(sceptre::construct_cis_pairs(
                    object,
                    distance_threshold = 20L,
                    response_position_data_frame = positions
                ))
            )
        )
    })
}

# ---------------------------------------------------------------------------
# Cached response fits and edge cases (#234)
# ---------------------------------------------------------------------------

consistency_discover <- function(object) {
    sceptre::run_discovery_analysis(
        object,
        print_progress = FALSE,
        output_amount = 2L
    )
}

consistency_case_cache_reuse <- function() {
    inputs <- consistency_synthetic_inputs()
    setting <- consistency_workflow_settings()$low_complement_perm
    formula <- stats::as.formula("~ batch + z", env = baseenv())
    discover <- consistency_discover
    consistency_capture({
        cold <- discover(consistency_prepare_synthetic(inputs, setting))
        warm <- discover(cold)
        cleared <- cold
        cleared@response_precomputations <- list()
        cleared <- discover(cleared)
        changed <- sceptre::set_analysis_parameters(
            cold,
            discovery_pairs = inputs$discovery_pairs[
                rev(seq_len(nrow(inputs$discovery_pairs))),
            ],
            positive_control_pairs = inputs$positive_control_pairs,
            formula_object = formula,
            control_group = setting$control,
            resampling_mechanism = setting$mechanism,
            resampling_approximation = setting$approximation,
            side = "right"
        )
        side_cache <- consistency_internal_part(
            changed@response_precomputations
        )
        changed <- discover(changed)
        formula_changed <- sceptre::set_analysis_parameters(
            cold,
            discovery_pairs = inputs$discovery_pairs,
            formula_object = stats::as.formula("~ batch", env = baseenv()),
            control_group = setting$control,
            resampling_mechanism = setting$mechanism
        )
        reassigned <- sceptre::assign_grnas(
            cold,
            method = "thresholding",
            threshold = 45,
            print_progress = FALSE
        )
        reassigned_cache <- consistency_internal_part(
            reassigned@response_precomputations
        )
        reassigned <- discover(consistency_qc_synthetic(reassigned, inputs))
        reassigned_cold <- discover(consistency_prepare_synthetic(
            inputs,
            setting,
            threshold = 45
        ))
        qc_changed <- consistency_qc_synthetic(cold, inputs, 25L)
        list(
            cold = consistency_observe(cold),
            warm = consistency_observe(warm),
            cleared = consistency_observe(cleared),
            changed_side = consistency_observe(changed),
            reassigned = consistency_observe(reassigned),
            reassigned_cold = consistency_observe(reassigned_cold),
            internal = consistency_internal_part(list(
                retained_side_cache = side_cache,
                formula_cache = formula_changed@response_precomputations,
                reassigned_cache = reassigned_cache,
                changed_qc_cache = qc_changed@response_precomputations
            ))
        )
    })
}

# Calibrate with NT controls, then switch to complement controls.
consistency_control_switch <- function(clear_fits) {
    inputs <- consistency_synthetic_inputs()
    setting <- consistency_workflow_settings()$low_nt_perm
    consistency_capture({
        object <- consistency_prepare_synthetic(inputs, setting)
        object <- sceptre::run_calibration_check(
            object,
            n_calibration_pairs = 4L,
            calibration_group_size = 2L,
            print_progress = FALSE,
            output_amount = 2L
        )
        changed <- sceptre::set_analysis_parameters(
            object,
            discovery_pairs = inputs$discovery_pairs,
            positive_control_pairs = inputs$positive_control_pairs,
            formula_object = stats::as.formula("~ batch + z", env = baseenv()),
            control_group = "complement",
            resampling_mechanism = "permutations"
        )
        n_retained_fits <- consistency_internal_part(
            length(changed@response_precomputations)
        )
        changed <- sceptre::assign_grnas(
            changed,
            method = "thresholding",
            threshold = 40,
            print_progress = FALSE
        )
        changed <- consistency_qc_synthetic(changed, inputs)
        if (clear_fits) {
            changed@response_precomputations <- list()
        }
        out <- consistency_observe(consistency_discover(changed))
        # The number of fits kept by the switch is what #12 was about (the
        # switch must discard them), so the uncleared run records it.
        if (!clear_fits && !is.null(n_retained_fits)) {
            out$internal$n_retained_fits <- n_retained_fits
        }
        out
    })
}

# Known bug (A2 #12): the switch keeps the NT-cell fits.
consistency_case_cache_control_switch <- function() {
    consistency_control_switch(clear_fits = FALSE)
}

# The same analysis with the fits cleared: what the fix should produce.
consistency_case_cache_control_switch_cleared <- function() {
    consistency_control_switch(clear_fits = TRUE)
}

consistency_case_edge_duplicate_pairs <- function() {
    inputs <- consistency_synthetic_inputs()
    inputs$discovery_pairs <- rbind(
        inputs$discovery_pairs,
        inputs$discovery_pairs[1L, ]
    )
    setting <- consistency_workflow_settings()$low_complement_perm
    consistency_capture({
        consistency_observe(consistency_prepare_synthetic(inputs, setting))
    })
}

consistency_case_edge_no_passing_pairs <- function() {
    inputs <- consistency_synthetic_inputs()
    inputs$discovery_pairs <- data.frame(grna_target = "A", response_id = "r4")
    inputs$positive_control_pairs <- inputs$positive_control_pairs[FALSE, ]
    setting <- consistency_workflow_settings()$low_complement_perm
    consistency_capture({
        object <- consistency_prepare_synthetic(inputs, setting)
        consistency_observe(sceptre::run_discovery_analysis(
            object,
            print_progress = FALSE
        ))
    })
}

# Known bug (A2 #11): `~ 1` with NT controls fails on the reference commit.
consistency_case_edge_intercept_only_nt <- function() {
    inputs <- consistency_synthetic_inputs()
    setting <- consistency_workflow_settings()$low_nt_perm
    consistency_capture({
        object <- consistency_prepare_synthetic(
            inputs,
            setting,
            formula = stats::as.formula("~ 1", env = baseenv())
        )
        consistency_observe(sceptre::run_discovery_analysis(
            object,
            print_progress = FALSE
        ))
    })
}

# ---------------------------------------------------------------------------
# Saved objects and public API
# ---------------------------------------------------------------------------

# A prepared low-MOI/complement/permutation object written by sceptre 0.99.0
# (21f9ea0), preserved from #234. It is a file-format input: never regenerate
# it. Its counts are #234's Poisson inputs, so its fits are the least stable
# across platforms (see README.md, "Cross-platform policy").
consistency_saved_object_file <- "sceptre-0.99.0-21f9ea0.rds"

consistency_case_saved_object_0_99_0 <- function(input_dir) {
    path <- file.path(
        input_dir,
        "saved-objects",
        consistency_saved_object_file
    )
    consistency_capture({
        object <- readRDS(path)
        loaded <- consistency_observe(object)
        object <- sceptre::run_discovery_analysis(
            object,
            print_progress = FALSE,
            output_amount = 2L
        )
        out <- consistency_observe(object)
        out$loaded <- loaded
        out$summary <- consistency_summary_text(object)
        out
    })
}

consistency_case_api <- function() {
    ns <- asNamespace("sceptre")
    exports <- getNamespaceExports(ns)
    exports <- sort(exports[!startsWith(exports, ".__")], method = "radix")
    describe_formals <- function(f) {
        if (!is.function(f)) {
            return(NULL)
        }
        lapply(formals(f), function(x) {
            paste(deparse(x, width.cutoff = 500L), collapse = "\n")
        })
    }
    methods_found <- lapply(c("plot", "print"), function(name) {
        lapply(methods::findMethods(name, where = ns), describe_formals)
    })
    names(methods_found) <- c("plot", "print")
    consistency_capture(list(
        exports = exports,
        formals = stats::setNames(
            lapply(exports, function(x) describe_formals(get(x, ns))),
            exports
        ),
        methods = methods_found,
        internal = list(
            slots = methods::getSlots(
                methods::getClass("sceptre_object", where = ns)
            )
        )
    ))
}

# ---------------------------------------------------------------------------
# Case registry
# ---------------------------------------------------------------------------

# Each case: group, description and a runner taking the fixture directory.
consistency_cases <- function() {
    case <- function(group, description, run) {
        list(group = group, description = description, run = run)
    }
    ignore <- function(f) function(input_dir) f()
    settings <- consistency_workflow_settings()
    workflow <- lapply(names(settings), function(name) {
        s <- settings[[name]]
        case(
            "workflow",
            paste0(
                "Synthetic 160-cell data, import to discovery: ",
                s$moi,
                " MOI, ",
                s$control,
                " controls, ",
                s$mechanism,
                ", ",
                s$integration,
                ", side ",
                s$side,
                ", ",
                s$approximation,
                ", output_amount ",
                s$output
            ),
            function(input_dir) consistency_run_workflow(s)
        )
    })
    names(workflow) <- paste0("workflow_", names(settings))
    c(
        list(
            api = case(
                "api",
                "Exports, their formals and S4 methods; sceptre_object slots",
                ignore(consistency_case_api)
            ),
            import_matrix_low = case(
                "import",
                "import_data() from dense synthetic matrices, low MOI",
                ignore(consistency_case_import_matrix_low)
            ),
            import_matrix_vector = case(
                "import",
                "import_data() with a vector_id column (guides summed)",
                ignore(consistency_case_import_matrix_vector)
            ),
            import_dense_symmetric = case(
                "import",
                "import_data() with a square symmetric dense response matrix",
                ignore(consistency_case_import_dense_symmetric)
            ),
            import_cellranger_two_directories = case(
                "import",
                "import_data_from_cellranger(), packaged two-directory example",
                ignore(consistency_case_import_cellranger_two_directories)
            ),
            import_cellranger_one_directory = case(
                "import",
                "import_data_from_cellranger(), one directory",
                ignore(consistency_case_import_cellranger_one_directory)
            ),
            import_cellranger_one_directory_covariates = case(
                "import",
                paste(
                    "import_data_from_cellranger(), one directory with",
                    "extra covariates `batch` and `lane` (covariates only)"
                ),
                ignore(
                    consistency_case_import_cellranger_one_directory_covariates
                )
            ),
            import_parse = case(
                "import",
                "import_data_from_parse(), packaged Parse example",
                ignore(consistency_case_import_parse)
            ),
            assign_low_maximum = case(
                "assignment",
                "assign_grnas() default (maximum), low-MOI example data",
                ignore(consistency_assignment_case("low"))
            ),
            assign_low_thresholding = case(
                "assignment",
                "assign_grnas(method = 'thresholding'), low-MOI example data",
                ignore(consistency_assignment_case(
                    "low",
                    method = "thresholding"
                ))
            ),
            assign_low_mixture = case(
                "assignment",
                "assign_grnas(method = 'mixture'), low-MOI example data",
                ignore(consistency_assignment_case("low", method = "mixture"))
            ),
            assign_high_mixture = case(
                "assignment",
                "assign_grnas() default (mixture), high-MOI example data",
                ignore(consistency_assignment_case("high"))
            ),
            assign_high_thresholding = case(
                "assignment",
                "assign_grnas(method = 'thresholding'), high-MOI example data",
                ignore(consistency_assignment_case(
                    "high",
                    method = "thresholding"
                ))
            )
        ),
        workflow,
        list(
            example_low_defaults = case(
                "workflow",
                paste(
                    "Low-MOI example data with package defaults: trans pairs,",
                    "maximum assignment, NT controls, permutations;",
                    "write_outputs_to_directory()"
                ),
                ignore(consistency_case_example_low_defaults)
            ),
            example_high_defaults = case(
                "workflow",
                paste(
                    "High-MOI example data with package defaults: cis pairs,",
                    "mixture assignment, complement controls, CRT"
                ),
                ignore(consistency_case_example_high_defaults)
            ),
            pairs_constructors = case(
                "pairs",
                "construct_positive_control/cis/trans_pairs(), high-MOI example",
                ignore(consistency_case_pairs_constructors)
            ),
            pairs_trans_caller_mutation = case(
                "pairs",
                "construct_trans_pairs(pairs_to_exclude = 'pc_pairs'): caller's table",
                ignore(consistency_case_pairs_trans_caller_mutation)
            ),
            pairs_cis_unsorted_positions = case(
                "pairs",
                "construct_cis_pairs() with an unsorted position table",
                ignore(consistency_case_pairs_cis_unsorted_positions)
            ),
            cache_reuse = case(
                "cache",
                "Reuse and invalidation of cached response fits",
                ignore(consistency_case_cache_reuse)
            ),
            cache_control_switch = case(
                "cache",
                "Switch from NT to complement controls after calibration",
                ignore(consistency_case_cache_control_switch)
            ),
            cache_control_switch_cleared = case(
                "cache",
                "The same switch with the cached fits cleared",
                ignore(consistency_case_cache_control_switch_cleared)
            ),
            edge_duplicate_pairs = case(
                "edge",
                "Duplicate discovery pairs are rejected",
                ignore(consistency_case_edge_duplicate_pairs)
            ),
            edge_no_passing_pairs = case(
                "edge",
                "Discovery when no pair passes pairwise QC",
                ignore(consistency_case_edge_no_passing_pairs)
            ),
            edge_intercept_only_nt = case(
                "edge",
                "Intercept-only formula with NT controls",
                ignore(consistency_case_edge_intercept_only_nt)
            ),
            saved_object_0_99_0 = case(
                "saved_object",
                "Released 0.99.0 saved object: load and run discovery",
                consistency_case_saved_object_0_99_0
            )
        )
    )
}

# Run `code` from a fixed seed, recording the RNG state it leaves behind.
consistency_seeded <- function(code_function) {
    withr::with_seed(consistency_case_seed, {
        observation <- code_function()
        observation$rng_after <- get(".Random.seed", envir = globalenv())
        observation
    })
}

# Run one case. Observations carry value, warnings, error and rng_after.
consistency_run_case <- function(name, input_dir) {
    cases <- consistency_cases()
    if (!name %in% names(cases)) {
        stop("Unknown consistency case: ", name)
    }
    data.table::setDTthreads(1L)
    consistency_seeded(function() cases[[name]]$run(input_dir))
}

# ---------------------------------------------------------------------------
# Public and internal parts of an observation
# ---------------------------------------------------------------------------

# The observation without its `internal` elements.
consistency_public <- function(x) {
    if (!is.list(x)) {
        return(x)
    }
    keep <- if (is.null(names(x))) {
        rep(TRUE, length(x))
    } else {
        names(x) != "internal"
    }
    out <- lapply(x[keep], consistency_public)
    attributes(out) <- attributes(x[keep])
    out
}

# Differences between two observations, as readable paths. Doubles agree when
# |actual - expected| <= tolerance * |expected| elementwise (an expected exact
# zero allows |actual| <= 1e-12, unless the tolerance is 0, which asks for
# equality), with the same NA/NaN/Inf positions. The
# comparison is relative at every magnitude, so a p-value of 1e-25 that moves
# to 1e-14 is a difference. Every other value must be identical.
consistency_differences <- function(
    actual,
    expected,
    tolerance,
    path = "observation"
) {
    if (!identical(class(actual), class(expected))) {
        return(paste0(path, ": class differs"))
    }
    if (is.list(expected)) {
        if (!identical(names(actual), names(expected))) {
            return(paste0(path, ": names differ"))
        }
        if (length(actual) != length(expected)) {
            return(paste0(path, ": length differs"))
        }
        other_attributes <- function(x) {
            a <- attributes(x)
            a[setdiff(names(a), "names")]
        }
        if (!identical(other_attributes(actual), other_attributes(expected))) {
            return(paste0(path, ": attributes differ"))
        }
        keys <- if (is.null(names(expected))) {
            paste0("[[", seq_along(expected), "]]")
        } else {
            paste0("$", names(expected))
        }
        return(unlist(lapply(seq_along(expected), function(k) {
            consistency_differences(
                actual[[k]],
                expected[[k]],
                tolerance,
                paste0(path, keys[[k]])
            )
        })))
    }
    if (!is.double(expected)) {
        if (identical(actual, expected)) {
            return(character())
        }
        return(paste0(path, ": value differs"))
    }
    if (!identical(attributes(actual), attributes(expected))) {
        return(paste0(path, ": attributes differ"))
    }
    special <- function(x) {
        ifelse(
            is.nan(x),
            "NaN",
            ifelse(
                is.na(x),
                "NA",
                ifelse(is.infinite(x), ifelse(x > 0, "Inf", "-Inf"), "finite")
            )
        )
    }
    if (!identical(special(actual), special(expected))) {
        return(paste0(path, ": NA/NaN/Inf pattern differs"))
    }
    finite <- is.finite(expected)
    a <- actual[finite]
    e <- expected[finite]
    gap <- abs(a - e)
    allowed <- ifelse(e == 0 & tolerance > 0, 1e-12, tolerance * abs(e))
    bad <- gap > allowed
    if (any(bad)) {
        worst <- which.max(gap / allowed)
        return(sprintf(
            "%s: %d value(s) differ (worst: actual %.17g, expected %.17g)",
            path,
            sum(bad),
            a[worst],
            e[worst]
        ))
    }
    character()
}

# The `internal` elements of an observation, keyed by their path.
consistency_internal <- function(x, path = "") {
    if (!is.list(x) || is.null(names(x))) {
        return(list())
    }
    found <- list()
    for (k in seq_along(x)) {
        name <- names(x)[[k]]
        child <- if (nzchar(path)) paste0(path, "$", name) else name
        if (identical(name, "internal")) {
            found[[child]] <- x[[k]]
        } else {
            found <- c(found, consistency_internal(x[[k]], child))
        }
    }
    found
}
