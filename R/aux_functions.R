#' Choose the existing default formula from source-cell covariates
#'
#' Retains column order; excludes response mito and optionally guide count
#' covariates. Names matching n_umis/n_nonzero get log(x), or log(x+1) if any
#' zero occurs. Other columns with at least 15 distinct values are omitted.
#' Returns a formula for model.matrix; this is selection, not input validation.
#' Does not modify covariates or draw random numbers.
#' @noRd
auto_construct_formula_object <- function(
    cell_covariates,
    include_grna_covariates
) {
    MAX_N_LEVELS_ALLOWED <- 15L
    cell_covariate_names <- colnames(cell_covariates)
    cell_covariate_names <- cell_covariate_names[
        cell_covariate_names != "response_p_mito"
    ]
    if (!include_grna_covariates) {
        # by default, do not use grna count-based covariates in low moi
        cell_covariate_names <- cell_covariate_names[
            !(cell_covariate_names %in% c("grna_n_umis", "grna_n_nonzero"))
        ]
    }
    form_str <- vapply(
        cell_covariate_names,
        function(curr_name) {
            count_based_covariate <- grepl(
                pattern = "n_umis|n_nonzero",
                x = curr_name
            )
            if (count_based_covariate) {
                if (any(cell_covariates[[curr_name]] == 0)) {
                    out <- paste0("log(", curr_name, "+1)")
                } else {
                    out <- paste0("log(", curr_name, ")")
                }
            } else {
                if (
                    length(unique(cell_covariates[[curr_name]])) >=
                        MAX_N_LEVELS_ALLOWED
                ) {
                    out <- NA_character_
                } else {
                    out <- curr_name
                }
            }
            return(out)
        },
        FUN.VALUE = character(1)
    ) |>
        stats::na.omit() |>
        paste0(collapse = " + ")
    form <- paste0("~ ", form_str) |> stats::as.formula()
    return(form)
}


#' Combine import summaries and extra covariates in source-cell order
#'
#' Aligned feature-by-cell response/guide counts and response labels produce a
#' data frame with prefixed response/guide summaries, followed by supplied extra
#' covariates. Response column IDs become row names when available. Converts
#' sparse orientation for column scans; does not create a backed store or use
#' RNG.
#' @noRd
auto_compute_cell_covariates <- function(
    response_matrix,
    grna_matrix,
    extra_covariates,
    response_names
) {
    # compute the response covariates
    covariate_df <- compute_cell_covariates(
        matrix_in = response_matrix,
        feature_names = response_names,
        compute_p_mito = TRUE,
        compute_max_feature = FALSE
    )
    colnames(covariate_df) <- paste0("response_", colnames(covariate_df))

    # compute the grna covariates
    grna_covariate_df <- compute_cell_covariates(
        matrix_in = grna_matrix,
        feature_names = character(0),
        compute_p_mito = TRUE,
        compute_max_feature = TRUE
    )
    colnames(grna_covariate_df) <- paste0("grna_", colnames(grna_covariate_df))
    covariate_df <- cbind(covariate_df, grna_covariate_df)

    # if extra covariates have been provided, add those as well
    if (nrow(extra_covariates) >= 1L) {
        covariate_df <- cbind(covariate_df, extra_covariates)
    }
    cell_ids <- colnames(response_matrix)
    if (!is.null(cell_ids)) {
        rownames(covariate_df) <- cell_ids
    }
    return(covariate_df)
}


#' Partition ordered labels or canonical row indices for existing fork workers
#'
#' Returns an ordered list of blocks. Parallel work with at least twice the
#' requested worker count uses the existing seed-4 shuffle and cut; smaller or
#' serial work retains one input-order block. "auto" uses the preferred core
#' fraction once. withr restores caller RNG/options; no worker is launched here.
#' @noRd
partition_response_ids <- function(response_ids, parallel, n_processors) {
    groups_set <- FALSE
    if (parallel) {
        if (identical(n_processors, "auto")) {
            # Read the user's preferred fraction (default 0.5 = half cores),
            # then
            # neutralize that option on the parallelly call so parallelly's own
            # option-aware default does not apply the same fraction a second
            # time.
            # parallelly's `fraction` arg defaults to getOption(
            # "parallelly.availableCores.fraction", 1), so omitting it would
            # silently double-scale when the user has set that option. The arg
            # itself was only added in 1.47.0, so passing `fraction = 1`
            # directly
            # would break parallelly < 1.47.0 -- `withr::with_options`
            # neutralizes
            # the option for the call on any parallelly >= 1.23.0.
            frac <- getOption("parallelly.availableCores.fraction", 0.5)
            n_processors <- max(
                1L,
                floor(
                    withr::with_options(
                        list(parallelly.availableCores.fraction = 1),
                        parallelly::availableCores(logical = FALSE)
                    ) *
                        frac
                )
            )
        }
        if (length(response_ids) >= 2 * n_processors) {
            s <- withr::with_seed(4, sample(response_ids))
            out <- split(
                s,
                cut(
                    seq_along(s),
                    n_processors,
                    labels = paste0("group_", seq(1, n_processors))
                )
            )
            groups_set <- TRUE
        }
    }
    if (!groups_set) {
        out <- list(group_1 = response_ids)
    }
    return(out)
}


#' Resolve the existing worker-log subdirectory, creating it if needed
#'
#' A scalar base path yields base/sceptre_logs/ with a trailing slash. Directory
#' creation is recursive; no log file is opened and no random state is changed.
#' @noRd
get_log_dir <- function(log_dir) {
    log_dir <- paste0(log_dir, "/sceptre_logs/")
    if (!dir.exists(log_dir)) {
        dir.create(log_dir, recursive = TRUE)
    }
    return(log_dir)
}
