#' Check count matrices, cell labels and metadata before importing
#' @param response_matrix,grna_matrix Feature-by-cell matrices in the supported
#'   numeric response or numeric/logical guide formats, or open ondisc matrices.
#' @param grna_target_data_frame Guide-to-target metadata. Supplied guide names
#'   must occur in the guide matrix. Keep metadata order and repeated rows.
#' @param moi `"low"` or `"high"`.
#' @param extra_covariates Optional covariates with rows in the same cell order.
#' @return NULL on success, otherwise an error at the first failing check below.
#'   Check classes, dimensions, labels and metadata fields without scanning every
#'   count value, converting storage or reordering cells.
#' @noRd
check_import_data_inputs <- function(
    response_matrix,
    grna_matrix,
    grna_target_data_frame,
    moi,
    extra_covariates
) {
    # check column names of grna_target_data_frame
    colnames_present <- all(
        c("grna_id", "grna_target") %in% colnames(grna_target_data_frame)
    )
    if (!colnames_present) {
        stop(
            "The data frame `grna_target_data_frame` should have columns `grna_id` and `grna_target`. The `grna_target` column should specify the target of each individual `grna_id`."
        )
    }

    # verify that the row names are unique for both response and grna
    # modalities
    response_ids <- rownames(response_matrix)
    grna_ids <- rownames(grna_matrix)
    if (length(response_ids) != length(unique(response_ids))) {
        stop("The rownames of the `response_matrix` must be unique.")
    }
    if (length(grna_ids) != length(unique(grna_ids))) {
        stop("The rownames of the `grna_matrix` must be unique.")
    }

    # ensure that the ampersand symbol (&) is absent from the grna ids;
    # ensure that no gRNA is named "non-targeting"
    problematic_grna_ids <- grep(pattern = "&", x = grna_ids)
    if (length(problematic_grna_ids) >= 1) {
        stop(paste0(
            "The ampersand character (&) cannot be present in the gRNA IDs. The following gRNA IDs contain an ampersand: ",
            paste0(grna_ids[problematic_grna_ids], collapse = ", ")
        ))
    }
    if (any(grna_ids == "non-targeting")) {
        stop(
            "No individual gRNA can have the ID `non-targeting`. The string `non-targeting` is reserved for the `grna_target` column of the `grna_target_data_frame`."
        )
    }

    # check that the ids in the grna group data frame are a subset of the ids
    # in the grna matrix
    if (!all(grna_target_data_frame$grna_id %in% rownames(grna_matrix))) {
        stop(
            "The column `grna_id` of the `grna_target_data_frame` must be a subset of the row names of the grna expression matrix. The row names of the grna expression matrix are as follows: ",
            paste0(rownames(grna_matrix), collapse = ", ")
        )
    }

    # check type of input matrices
    # Check that one input matrix has a supported class
    # input_matrix: Matrix to check; `input_matrix_name` identifies it in errors.
    # allowed_matrix_classes: Supported class names. Do not convert the input.
    # Returns: Invisible NULL on success, otherwise an error naming the input matrix.
    check_matrix_class <- function(
        input_matrix,
        input_matrix_name,
        allowed_matrix_classes
    ) {
        ok_class <- vapply(
            X = allowed_matrix_classes,
            function(mat_class) methods::is(input_matrix, mat_class),
            FUN.VALUE = logical(1)
        ) |>
            any()
        if (!ok_class) {
            stop(paste0(
                "`",
                input_matrix_name,
                "` must be an object of class ",
                paste0(allowed_matrix_classes, collapse = ", "),
                "."
            ))
        }
    }
    check_matrix_class(
        response_matrix,
        "response_matrix",
        c("matrix", "dgTMatrix", "dgCMatrix", "dgRMatrix", "odm")
    )
    check_matrix_class(
        grna_matrix,
        "grna_matrix",
        c(
            "matrix",
            "dgTMatrix",
            "dgCMatrix",
            "dgRMatrix",
            "lgTMatrix",
            "lgCMatrix",
            "lgRMatrix",
            "odm"
        )
    )

    # check for agreement in number of cells
    check_ncells <- ncol(response_matrix) == ncol(grna_matrix)
    if (nrow(extra_covariates) >= 1L) {
        check_ncells <- check_ncells &&
            (ncol(response_matrix) == nrow(extra_covariates))
    }
    if (!check_ncells) {
        stop(
            "The number of cells in the `response_matrix`, `grna_matrix`, and `extra_covariates` (if supplied) must coincide."
        )
    }

    # if cell barcodes are provided for at least two of `response_matrix`,
    # `grna_matrix`, and `extra_covariates`,
    # then they must be identical and in the same order
    barcode_list <- list(
        response_matrix = colnames(response_matrix),
        grna_matrix = colnames(grna_matrix),
        extra_covariates = rownames(extra_covariates)
    )
    # for matrices we can just check if the names are not NULL, but
    # `extra_covariates` is a data.frame so
    # non-default names were provided if (1) it is not just `data.frame()` and
    # (2) the names are not
    # just "1", "2", ...
    were_names_provided <- c(
        !is.null(barcode_list$response_matrix),
        !is.null(barcode_list$grna_matrix),
        nrow(extra_covariates) > 0 &&
            !identical(
                barcode_list$extra_covariates,
                as.character(seq_len(nrow(extra_covariates)))
            )
    )
    # If at least 2 non-default barcode names were provided, they must all be
    # identical and in the same order.
    # This is done by looping over all pairs of non-default names
    barcodes_with_names <- barcode_list[were_names_provided]
    if (length(barcodes_with_names) >= 2) {
        for (i in seq_len(length(barcodes_with_names) - 1)) {
            for (j in seq(i + 1, length(barcodes_with_names))) {
                if (
                    !identical(
                        barcodes_with_names[[i]],
                        barcodes_with_names[[j]]
                    )
                ) {
                    stop(
                        "You have provided cell barcodes in the `",
                        names(barcodes_with_names)[i],
                        "` and `",
                        names(barcodes_with_names)[j],
                        "`. These cell barcodes must be identical and in the same order across objects."
                    )
                }
            }
        }
    }

    # check that column names of extra_covariates are not already taken
    reserved_covariate_names <- c(
        "response_n_nonzero",
        "response_n_umis",
        "response_p_mito",
        "grna_n_nonzero",
        "grna_n_umis"
    )
    extra_covariate_names <- colnames(extra_covariates)
    if (any(extra_covariate_names %in% reserved_covariate_names)) {
        stop(
            "The covariate names `response_n_nonzero`, `response_n_umis`, `response_p_mito`, `grna_n_nonzero`, and `grna_n_umis` are reserved. Change the column names of the `extra_covariates` data frame."
        )
    }

    # verify that the types of the extra covariates are acceptable; if
    # factor, check correctness of factor
    for (extra_covariate_name in extra_covariate_names) {
        v <- extra_covariates[, extra_covariate_name]
        accept_type <- methods::is(v, "numeric") ||
            methods::is(v, "character") ||
            methods::is(v, "factor") ||
            methods::is(v, "logical")
        if (!accept_type) {
            stop(
                "The column `",
                extra_covariate_name,
                "` of the `extra_covariates` data frame should be of type numeric, character, or factor."
            )
        }
        if (methods::is(v, "factor")) {
            if (length(levels(v)) != length(unique(v))) {
                stop(
                    "The factor column `",
                    extra_covariate_name,
                    "` is incorrectly formatted. The number of unique elements in this column must equal the number of factors in this column. Consider relabeling the factors or converting this column into a character vector."
                )
            }
        }
    }

    # verify that moi is specified
    if (!(moi %in% c("low", "high"))) {
        stop("`moi` should be either `low` or `high`.")
    }

    # Reject missing or infinite extra covariate values.
    if (any(is.na(extra_covariates))) {
        stop("`extra_covariates` has NA values that need to be removed.")
    }
    if (
        vapply(
            extra_covariates,
            function(vect) any(is.infinite(vect)),
            FUN.VALUE = logical(1)
        ) |>
            any()
    ) {
        stop("`extra_covariates` has infinite values that need to be removed.")
    }

    # Fail if grna_target_data_frame has NA values in grna_id or grna_target
    # columns
    if (any(is.na(grna_target_data_frame[, c("grna_id", "grna_target")]))) {
        stop(
            "The `grna_id` or `grna_target` column of `grna_target_data_frame` contains NA values."
        )
    }

    return(NULL)
}


#' Check proposed analysis settings against the imported data
#' @param sceptre_object Imported sceptre object with count matrices and covariates.
#' @param formula_object Model formula; offsets are not supported.
#' @param response_grna_target_pairs_list Named discovery and positive-control tables.
#' @param control_group,resampling_mechanism,side,low_moi Chosen control cells,
#'   resampling method, test direction and MOI setting.
#' @param grna_integration_strategy,resampling_approximation Chosen guide-combination
#'   method and resampling approximation.
#' @return NULL on success, otherwise the first error in the order below.
#'   Check pair names and duplicate complete rows without changing their order.
#'   The current checks access missing-value columns before checking their names.
#'   ondisc-backed matrices currently support permutations only. This helper
#'   checks settings; building the design matrix happens separately.
#' @noRd
check_set_analysis_parameters <- function(
    sceptre_object,
    formula_object,
    response_grna_target_pairs_list,
    control_group,
    resampling_mechanism,
    side,
    low_moi,
    grna_integration_strategy,
    resampling_approximation
) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    response_matrix <- get_response_matrix(sceptre_object)
    covariate_data_frame <- sceptre_object@data@covariate_data_frame
    grna_target_data_frame <- guide_group_metadata(sceptre_object)

    # if response_grna_target_pairs has been supplied, check its
    # characteristics
    for (idx in seq_along(response_grna_target_pairs_list)) {
        response_grna_target_pairs <- response_grna_target_pairs_list[[idx]]
        if (nrow(response_grna_target_pairs) >= 1L) {
            df_name <- names(response_grna_target_pairs_list)[idx]
            # ensure no NA
            if (
                any(is.na(response_grna_target_pairs[, c(
                    "grna_target",
                    "response_id"
                )]))
            ) {
                stop(
                    "The `grna_target` or `response_id` column of the `",
                    df_name,
                    "` data frame contains NA values."
                )
            }
            # i. verify that `grna_target` and `response_id` are columns
            if (
                !all(
                    c("grna_target", "response_id") %in%
                        colnames(response_grna_target_pairs)
                )
            ) {
                stop(
                    "The data frame `",
                    df_name,
                    "` must contain the columns `grna_target` and `response_id`."
                )
            }
            # ii. check that the response ids in the
            # `response_grna_target_pairs` data frame are a subset of the
            # response ids
            if (
                !all(
                    response_grna_target_pairs$response_id %in%
                        rownames(response_matrix)
                )
            ) {
                stop(
                    "The column `response_id` of the `",
                    df_name,
                    "` data frame must be a subset of the row names of the response expression matrix."
                )
            }
            # iii. check that the `grna_target` column of the
            # `response_grna_target_pairs` data frame is a subset of the
            # `grna_target` column of the `grna_target_data_frame`
            if (
                !all(
                    response_grna_target_pairs$grna_target %in%
                        grna_target_data_frame$grna_target
                )
            ) {
                stop(
                    "The column `grna_target` of the `",
                    df_name,
                    "` data frame must be a subset of the colummn `grna_target` of the `grna_target_data_frame`."
                )
            }
            # iv. ensure that "non-targeting" is not a group in the pairs to
            # analyze data frame
            if (
                "non-targeting" %in%
                    unique(response_grna_target_pairs$grna_target)
            ) {
                stop(
                    "The `",
                    df_name,
                    "` data frame cannot contain the gRNA group `non-targeting`. To test non-targeting gRNAs against responses, run the calibration check."
                )
            }
            # v. ensure that rows are not duplicated
            if (
                nrow(response_grna_target_pairs) !=
                    nrow(dplyr::distinct(response_grna_target_pairs))
            ) {
                stop("The `", df_name, "` data frame has duplicate rows.")
            }
        }
    }

    # check that there are no offsets in the formula object
    if (formula_contains_offset_call(formula_object)) {
        stop("Offsets are not currently supported in formula objects.")
    }

    # check that the variables in the formula object are a subset of the
    # column names of the covariate data frame
    formula_object_vars <- all.vars(formula_object)
    check_var <- formula_object_vars %in% colnames(covariate_data_frame)
    if (!all(check_var)) {
        stop(paste0(
            "The variables in the `formula_object` must be a subset of the columns of the `covariate_data_frame`. Check the following variables: ",
            paste0(formula_object_vars[!check_var], collapse = ", ")
        ))
    }

    # verify that control group is either nt_cells, complement, or default
    if (!control_group %in% c("nt_cells", "complement", "default")) {
        stop(
            "`control_group` should set to `nt_cells`, `complement`, or `default`."
        )
    }

    # verify that resampling_mechanism is one of "permutations", "crt", or
    # "default"
    if (!(resampling_mechanism %in% c("permutations", "crt", "default"))) {
        stop(
            "`resampling_mechanism` should set to `permutations`, `crt`, or `default`."
        )
    }

    # verify that the moi is consistent with the control group
    if (!low_moi && control_group == "nt_cells") {
        stop("The control group cannot be the NT cells in high MOI.")
    }

    # verify that, if the control_group is NT cells, there are NT gRNAs
    # present
    if (control_group == "nt_cells") {
        nt_present <- "non-targeting" %in% grna_target_data_frame$grna_target
        if (!nt_present) {
            stop(
                "The string 'non-targeting' must be present in the `grna_target` column of the `grna_target_data_frame` is `control_group` is set to 'nt_cells'."
            )
        }
    }

    # verify that "side" is among "both", "left", or "right"
    if (!(side %in% c("both", "left", "right"))) {
        stop("`side` must be one of 'both', 'left', or 'right'.")
    }

    # Check that the guide integration strategy is supported.
    if (
        !(grna_integration_strategy %in% c("union", "singleton", "bonferroni"))
    ) {
        stop(
            "`grna_integration_strategy` must be either 'union', 'singleton', or 'bonferroni'."
        )
    }

    # if using a backing .odm file, verify that resampling mechanism is
    # permutations
    if (
        methods::is(get_response_matrix(sceptre_object), "odm") &&
            (resampling_mechanism != "permutations")
    ) {
        stop(
            "`resampling_mechanism` must be set to 'permutations' when using an ondisc-backed sceptre_object."
        )
    }

    # verify resampling_approximation acceptable
    if (!(resampling_approximation %in% c("skew_normal", "no_approximation"))) {
        stop(
            "`resampling_approximation` must be set to 'skew_normal' or 'no_approximation'."
        )
    }

    return(NULL)
}


#' Check the parallel option and whether the platform supports forked workers
#' @param parallel One nonmissing logical value.
#' @param os_type Operating-system type; supplied explicitly in platform tests.
#' @return NULL, or an error for an invalid option or parallel work on Windows.
#'   This checks support without starting workers or changing platform settings.
#' @noRd
check_parallel_supported <- function(parallel, os_type = .Platform$OS.type) {
    if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
        stop("`parallel` should be a single logical value.")
    }
    if (parallel && identical(os_type, "windows")) {
        stop(
            "`parallel = TRUE` is not supported on Windows because sceptre's parallel mode uses fork-based processing. ",
            "Please set `parallel = FALSE`.",
            call. = FALSE
        )
    }
    return(NULL)
}


#' Check the requested worker count using the current validation rule
#' @param n_processors Exactly `"auto"` or a numeric value at least 2.
#' @return NULL on success, otherwise an error. The existing comparison also accepts
#'   fractional values and Inf; this helper does not round them or choose workers.
#'   Call it after the other argument checks, as in the public workflows.
#' @noRd
check_processor_count <- function(n_processors) {
    if (
        !(identical(n_processors, "auto") ||
            (is.numeric(n_processors) && n_processors >= 2))
    ) {
        stop(
            "`n_processors` should be set to the string 'auto' or an integer greater than or equal to 2."
        )
    }
    NULL
}


#' Check guide-assignment method and settings before assigning cells
#' @param sceptre_object Object supplying the low/high-MOI setting.
#' @param assignment_method Chosen maximum, thresholding or mixture method.
#' @param hyperparameters Named method settings after merging user values and defaults.
#' @param n_processors Requested worker count, checked after method settings.
#' @return NULL on success, otherwise the first method, range or setting-name error
#'   below. The perturbation parameter uses its documented `[0, 20]` range. Other
#'   type, missing-value and scalar checks retain their current rules. No models
#'   are fitted and no assignments are changed here.
#' @noRd
check_assign_grna_inputs <- function(
    sceptre_object,
    assignment_method,
    hyperparameters,
    n_processors
) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    if (!(assignment_method %in% c("maximum", "mixture", "thresholding"))) {
        stop(
            "`assignment_method` must be `mixture`, `maximum`, or `thresholding`."
        )
    }

    # check that assignment method is OK for high MOI
    if (!sceptre_object@data@low_moi && assignment_method == "maximum") {
        stop("`assignment_method` cannot be `maximum` in high-MOI.")
    }

    # check the hyperparameters
    hyperparam_names <- names(hyperparameters)
    # i. maximum option
    if (assignment_method == "maximum") {
        if (
            !setequal(
                hyperparam_names,
                c("umi_fraction_threshold", "min_grna_n_umis_threshold")
            )
        ) {
            stop(
                "The hyperparameter list must contain the fields `umi_fraction_threshold` and `min_grna_n_umis_threshold`."
            )
        }
        if (
            !(hyperparameters[["umi_fraction_threshold"]] > 0.0 &&
                hyperparameters[["umi_fraction_threshold"]] < 1.0)
        ) {
            stop(
                "`umi_fraction_threshold` should be a numeric greater than 0.0 and be less than 1.0."
            )
        }
        if (hyperparameters[["min_grna_n_umis_threshold"]] < 0L) {
            stop(
                "`min_grna_n_umis_threshold` should be an integer greater than 0."
            )
        }
    }
    # ii. thresholding operation
    if (assignment_method == "thresholding") {
        if (!setequal(hyperparam_names, "threshold")) {
            stop(
                "The hyperparameter list must contain the single entry `threshold`."
            )
        }
        if (hyperparameters[["threshold"]] < 1) {
            stop("`threshold` should be an integer greater than or equal to 1.")
        }
    }
    # iii. mixture model
    if (assignment_method == "mixture") {
        if (
            !setequal(
                hyperparam_names,
                c(
                    "n_em_rep",
                    "pi_guess_range",
                    "g_pert_guess_range",
                    "n_nonzero_cells_cutoff",
                    "backup_threshold",
                    "probability_threshold",
                    "formula_object"
                )
            )
        ) {
            stop(
                "The hyperparameter list must contain the fields `n_em_rep`, `pi_guess_range`, `g_pert_guess_range`, `n_nonzero_cells_cutoff`, `backup_threshold`, and `formula_object`."
            )
        }
        if (
            hyperparameters[["n_em_rep"]] <= 0 ||
                hyperparameters[["n_em_rep"]] >= 100
        ) {
            stop(
                "`n_em_rep` should be an integer greater than zero and less than 100."
            )
        }
        pi_guess_range <- hyperparameters[["pi_guess_range"]]
        if (
            length(pi_guess_range) != 2 ||
                any(pi_guess_range < 0) ||
                any(pi_guess_range > 1) ||
                pi_guess_range[2] <= pi_guess_range[1]
        ) {
            stop(
                "`pi_guess_range` should be a numeric vector of length 2, where the first entry is less than the second and both entries are in the interval [0,1]."
            )
        }
        g_pert_guess_range <- hyperparameters[["g_pert_guess_range"]]
        if (
            length(g_pert_guess_range) != 2 ||
                any(g_pert_guess_range < 0) ||
                any(g_pert_guess_range > 20) ||
                g_pert_guess_range[2] <= g_pert_guess_range[1]
        ) {
            stop(
                "`g_pert_guess_range` should be a numeric vector of length 2, where the first entry is less than the second and both entries are in the interval [0,20]."
            )
        }
        if (hyperparameters[["n_nonzero_cells_cutoff"]] <= 0) {
            stop(
                "`n_nonzero_cells_cutoff` should be an integer greater than or equal to 1."
            )
        }
        if (hyperparameters[["backup_threshold"]] <= 0) {
            stop(
                "`backup_threshold` should be an integer greater than or equal to 1."
            )
        }
        if (
            hyperparameters[["probability_threshold"]] <= 0 ||
                hyperparameters[["probability_threshold"]] >= 1
        ) {
            stop(
                "`probability_threshold` should be a numeric in the interval (0,1)."
            )
        }
    }

    # check n_processors argument
    check_processor_count(n_processors)

    return(NULL)
}


#' Check pair-count thresholds, QC quantiles and whether any guide was assigned
#' @param n_nonzero_trt_thresh,n_nonzero_cntrl_thresh Nonnegative pair-QC cutoffs.
#' @param response_n_umis_range,response_n_nonzero_range Ordered quantile pairs.
#' @param initial_grna_assignment_list Named lists of original-cell assignments.
#' @return NULL on success, otherwise the first error below. These are the checks
#'   performed here, rather than validation of every public QC argument. This
#'   function does not remove cells or count gene-guide pairs.
#' @noRd
check_run_qc_inputs <- function(
    n_nonzero_trt_thresh,
    n_nonzero_cntrl_thresh,
    response_n_umis_range,
    response_n_nonzero_range,
    initial_grna_assignment_list
) {
    if (n_nonzero_trt_thresh < 0 || n_nonzero_cntrl_thresh < 0) {
        stop(
            "`n_nonzero_trt_thresh` and `n_nonzero_cntrl_thresh` must be greater than or equal to zero."
        )
    }
    if (
        min(response_n_umis_range) < 0 ||
            max(response_n_umis_range) > 1 ||
            length(response_n_umis_range) != 2L ||
            response_n_umis_range[1] > response_n_umis_range[2]
    ) {
        stop("`response_n_umis_range` must an interval in the range [0,1].")
    }
    if (
        min(response_n_nonzero_range) < 0 ||
            max(response_n_nonzero_range) > 1 ||
            length(response_n_nonzero_range) != 2L ||
            response_n_nonzero_range[1] > response_n_nonzero_range[2]
    ) {
        stop("`response_n_nonzero_range` must an interval in the range [0,1].")
    }
    if (
        all(
            vapply(
                initial_grna_assignment_list,
                length,
                FUN.VALUE = integer(1)
            ) ==
                0L
        )
    ) {
        stop(
            "At least one gRNA must be assigned to at least one cell to call `run_qc()`."
        )
    }
    return(NULL)
}


#' Check that a single count can safely be converted to a C++ integer
#'
#' Require one finite, whole numeric value between `minimum` and INT_MAX, the
#' largest C++ integer. Reject invalid values rather than rounding or converting
#' them first. Return the input invisibly. Zero is allowed for unused resampling
#' stages when minimum is zero.
#' @noRd
check_native_scalar_count <- function(value, name, minimum = 0) {
    if (
        !is.numeric(value) ||
            length(value) != 1L ||
            is.na(value) ||
            !is.finite(value) ||
            value != floor(value) ||
            value < minimum ||
            value > .Machine$integer.max
    ) {
        stop(
            "`",
            name,
            "` must be a finite integer between ",
            minimum,
            " and ",
            .Machine$integer.max,
            "."
        )
    }
    invisible(value)
}

#' Check B1, B2 and B3 and return the total number of resamples
#'
#' Each stage size must be a nonnegative whole number. Add using doubles first,
#' then check that the total also fits a C++ integer before converting it. Report
#' an error if it is too large; do not silently reduce the number of resamples.
#' @noRd
checked_resampling_total <- function(B1, B2, B3) {
    check_native_scalar_count(B1, "B1")
    check_native_scalar_count(B2, "B2")
    check_native_scalar_count(B3, "B3")
    total <- sum(as.double(B1), as.double(B2), as.double(B3))
    check_native_scalar_count(total, "B1 + B2 + B3")
    as.integer(total)
}

#' Calculate B3 for empirical p-values and check the total resampling size
#'
#' For a nonnegative whole pair count, use ceiling(10 * n_pairs / alpha) for a
#' two-sided test or ceiling(5 * n_pairs / alpha) for a one-sided test. Check both
#' B3 and B1+B2+B3 before converting B3 to an integer. Large requests error
#' rather than being silently capped.
#' @noRd
compute_empirical_budget <- function(n_pairs, alpha, side_code, B1, B2) {
    check_native_scalar_count(n_pairs, "n_pairs")
    multiplier <- if (side_code == 0L) 10 else 5
    budget <- ceiling(multiplier * as.double(n_pairs) / alpha)
    check_native_scalar_count(budget, "Derived B3")
    checked_resampling_total(B1, B2, budget)
    as.integer(budget)
}

#' Check calibration sizes against NT guides that still have cells after QC
#' @param sceptre_object Post-QC object. Count active NT guides from assignments,
#'   since metadata can still list guides whose last cell was removed by QC.
#' @param n_calibration_pairs Positive whole number fitting a C++ integer.
#' @param n_processors Requested worker count, checked last.
#' @param calibration_group_size Positive number of NT guides per group, no larger
#'   than the active NT-guide count. NT controls require at least two active NT
#'   guides; complement controls require at least one.
#' @return NULL on success, otherwise the first error below. Sizes are checked
#'   without rounding them, constructing pairs or drawing random samples.
#' @noRd
check_calibration_check_inputs <- function(
    sceptre_object,
    n_calibration_pairs,
    n_processors,
    calibration_group_size
) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    control_group_complement <- sceptre_object@spec@control_group_complement
    # QC can remove every cell assigned to a guide while its metadata remains.
    n_nt_grnas <- sum(
        lengths(sceptre_object@preparation@analysis_groups$indiv_nt_grna_idxs) > 0L
    )
    # check number of gRNAS
    if (!control_group_complement) {
        if (n_nt_grnas < 2) {
            stop(
                "Two or more non-targeting gRNAs must be present to run the calibration check when using the NT cells as the control group. gRNAs that are non-targeting should be assigned a gRNA group label of 'non-targeting' in the `grna_target_data_frame`."
            )
        }
    } else {
        if (n_nt_grnas < 1) {
            stop(
                "At least one non-targeting gRNA must be present to run the calibration check. gRNAs that are non-targeting should be assigned a gRNA group label of 'non-targeting' in the `grna_target_data_frame`."
            )
        }
    }
    # check number of calibration pairs
    check_native_scalar_count(n_calibration_pairs, "n_calibration_pairs")
    if (n_calibration_pairs == 0L) {
        stop("Cannot run a calibration check on zero negative control pairs.")
    }
    check_native_scalar_count(
        calibration_group_size,
        "calibration_group_size",
        1
    )
    if (calibration_group_size > n_nt_grnas) {
        stop(
            "`calibration_group_size` exceeds the active non-targeting gRNA count after QC (",
            n_nt_grnas,
            ")."
        )
    }
    # check n_processors
    check_processor_count(n_processors)
    return(NULL)
}


#' Check that discovery or power analysis has pairs and usable control cells
#' @param response_grna_group_pairs Requested pair table with QC information.
#' @param control_group_complement Whether controls are all other kept cells.
#' @param grna_target_data_frame Metadata used to check for NT guides when required.
#' @param pc_analysis TRUE for the positive-control power check; FALSE for discovery.
#' @param calibration_result Earlier calibration results; an empty table prompts
#'   the advisory message to run calibration first.
#' @param n_ok_pairs Number of pairs passing QC; zero stops the analysis.
#' @param n_processors Requested worker count, checked last.
#' @return NULL on success, otherwise messages/errors in the order below: pairs,
#'   NT guides, calibration advisory, QC-passing count, then worker count.
#'   This function checks readiness without filtering or running the pairs.
#' @noRd
check_discovery_analysis_inputs <- function(
    response_grna_group_pairs,
    control_group_complement,
    grna_target_data_frame,
    pc_analysis,
    calibration_result,
    n_ok_pairs,
    n_processors
) {
    # check that positive control pairs are available
    if (nrow(response_grna_group_pairs) == 0L) {
        stop(
            if (pc_analysis) "Positive control" else "Discovery",
            " pairs have not been supplied. Thus, the ",
            if (pc_analysis) "power check" else "discovery analysis",
            " cannot be run. You can supply ",
            if (pc_analysis) "positive control" else "discovery",
            " pairs in the function set_analysis_parameters()."
        )
    }

    # check that negative control gRNAs are present (if the control group is
    # the nt cells)
    if (!control_group_complement) {
        nt_present <- "non-targeting" %in% grna_target_data_frame$grna_group
        if (!nt_present) {
            stop(
                "At least one non-targeting gRNA must be present to run a ",
                if (pc_analysis) "power check" else "discovery analysis",
                " when the control group is the NT set."
            )
        }
    }

    # check for the presence of a calibration result
    if (nrow(calibration_result) == 0L) {
        warn_text <- paste0(
            "The calibration check (`run_calibration_check()`) should be run before the ",
            ifelse(pc_analysis, "power check", "discovery analysis"),
            "."
        )
        message(crayon::red(warn_text))
    }

    # check that at least one pair passes qc
    if (n_ok_pairs == 0L) {
        stop("Zero pairs pass pairwise QC. Cannot run analysis.")
    }

    # check n_processors
    check_processor_count(n_processors)

    return(NULL)
}

# Detect actual offset calls rather than covariate labels containing the word.
#' @noRd
formula_contains_offset_call <- function(expression) {
    if (!is.call(expression)) return(FALSE)
    head <- expression[[1L]]
    offset <- identical(head, as.name("offset")) ||
        (is.call(head) && length(head) == 3L &&
            as.character(head[[1L]]) %in% c("::", ":::") &&
            identical(head[[3L]], as.name("offset")))
    offset || any(vapply(as.list(expression)[-1L], formula_contains_offset_call, logical(1)))
}
