construct_data_frame_v2 <- function(
    curr_df,
    curr_response_result,
    output_amount
) {
    curr_df$p_value <- vapply(
        X = curr_response_result,
        FUN = function(l) l$p,
        FUN.VALUE = numeric(1)
    )
    fold_change <- vapply(
        curr_response_result,
        FUN = function(l) l$fc,
        FUN.VALUE = numeric(1)
    )
    curr_df$fold_change <- fold_change
    curr_df$se_fold_change <- vapply(
        curr_response_result,
        FUN = function(l) l$se,
        FUN.VALUE = numeric(1)
    )
    curr_df$log_2_fold_change <- log(fold_change) / log(2)
    if (output_amount >= 2L) {
        curr_df$stage <- vapply(
            curr_response_result,
            FUN = function(l) l$stage,
            FUN.VALUE = integer(1)
        )
        curr_df$z_orig <- vapply(
            curr_response_result,
            FUN = function(l) l$z_orig,
            FUN.VALUE = numeric(1)
        )
        curr_df$xi <- vapply(
            curr_response_result,
            FUN = function(l) l$sn_params[1L],
            FUN.VALUE = numeric(1)
        )
        curr_df$omega <- vapply(
            curr_response_result,
            FUN = function(l) l$sn_params[2L],
            FUN.VALUE = numeric(1)
        )
        curr_df$alpha <- vapply(
            curr_response_result,
            FUN = function(l) l$sn_params[3L],
            FUN.VALUE = numeric(1)
        )
    }
    if (output_amount >= 3L) {
        to_append <- lapply(curr_response_result, FUN = function(l) {
            m <- data.table::as.data.table(matrix(l$resampling_dist, nrow = 1))
            colnames(m) <- paste0("z_null_", seq(1L, ncol(m)))
            m
        }) |>
            data.table::rbindlist(fill = TRUE)
        curr_df <- cbind(curr_df, to_append)
    }
    return(curr_df)
}


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
    )
    left_out <- names(form_str)[is.na(form_str)]
    if (length(left_out) >= 1L) {
        message(
            "The default formula leaves out the following covariates, which ",
            "have ",
            MAX_N_LEVELS_ALLOWED,
            " or more distinct values: ",
            paste0("`", left_out, "`", collapse = ", "),
            ". To include them, supply a `formula_object`."
        )
    }
    form_str <- form_str |>
        stats::na.omit() |>
        paste0(collapse = " + ")
    form <- paste0("~ ", form_str) |> stats::as.formula()
    return(form)
}


auto_compute_cell_covariates <- function(
    response_matrix,
    grna_matrix,
    extra_covariates,
    response_names,
    response_covariates = NULL
) {
    # compute the response covariates (unless computed already, as for
    # response counts held on disk)
    covariate_df <- if (is.null(response_covariates)) {
        compute_cell_covariates(
            matrix_in = response_matrix,
            feature_names = response_names,
            compute_p_mito = TRUE,
            compute_max_feature = FALSE
        )
    } else {
        response_covariates
    }
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


# Return the seed of an analysis: `seed` itself, checked, or, when `seed` is
# NULL, a seed drawn from R's random number generator, so that set.seed()
# governs it.
resolve_seed <- function(seed) {
    if (is.null(seed)) {
        return(sample.int(.Machine$integer.max, 1L))
    }
    if (
        !is.numeric(seed) ||
            length(seed) != 1L ||
            !is.finite(seed) ||
            seed != round(seed) ||
            seed < 0 ||
            seed > .Machine$integer.max
    ) {
        stop(
            "`seed` must be NULL or a single whole number between 0 and ",
            .Machine$integer.max,
            ".",
            call. = FALSE
        )
    }
    as.integer(seed)
}


# Return the BiocParallel backend for an analysis function, or NULL to run
# serially. A user-supplied `BPPARAM` takes precedence; otherwise
# `parallel = TRUE` forks `n_processors` workers, as in earlier versions.
get_bpparam <- function(parallel, n_processors, BPPARAM) {
    if (!is.null(BPPARAM)) {
        if (!methods::is(BPPARAM, "BiocParallelParam")) {
            stop(
                "`BPPARAM` should be a `BiocParallelParam` object, such as ",
                "`BiocParallel::MulticoreParam()` or ",
                "`BiocParallel::SnowParam()`."
            )
        }
        if (BiocParallel::bpnworkers(BPPARAM) < 2L) {
            return(NULL)
        }
        return(BPPARAM)
    }
    check_parallel_supported(parallel)
    if (!parallel) {
        return(NULL)
    }
    if (
        !(identical(n_processors, "auto") ||
            (is.numeric(n_processors) && n_processors >= 2))
    ) {
        stop(
            "`n_processors` should be set to the string 'auto' or an integer greater than or equal to 2."
        )
    }
    if (identical(n_processors, "auto")) {
        # Read the user's preferred fraction (default 0.5 = half cores), then
        # neutralize that option on the parallelly call so parallelly's own
        # option-aware default does not apply the same fraction a second
        # time.
        # parallelly's `fraction` arg defaults to getOption(
        # "parallelly.availableCores.fraction", 1), so omitting it would
        # silently double-scale when the user has set that option. The arg
        # itself was only added in 1.47.0, so passing `fraction = 1` directly
        # would break parallelly < 1.47.0 -- `withr::with_options`
        # neutralizes the option for the call on any parallelly >= 1.23.0.
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
    BiocParallel::MulticoreParam(workers = as.integer(n_processors))
}


# Split the IDs into one group per worker of `BPPARAM` (a single group when
# running serially). The IDs are dealt out in turn, so that neighboring IDs,
# which may have similar workloads, go to different workers.
partition_response_ids <- function(response_ids, BPPARAM) {
    n_processors <- if (is.null(BPPARAM)) {
        1L
    } else {
        BiocParallel::bpnworkers(BPPARAM)
    }
    if (n_processors >= 2L && length(response_ids) >= 2 * n_processors) {
        group <- (seq_along(response_ids) - 1L) %% n_processors + 1L
        out <- split(
            response_ids,
            factor(
                paste0("group_", group),
                levels = paste0("group_", seq(1, n_processors))
            )
        )
    } else {
        out <- list(group_1 = response_ids)
    }
    return(out)
}


# Run `FUN` on each element of `X` with the BiocParallel backend `BPPARAM`,
# and stop with a worker's original error if any worker failed.
bplapply_checked <- function(X, FUN, BPPARAM) {
    res <- tryCatch(
        BiocParallel::bptry(BiocParallel::bplapply(X, FUN, BPPARAM = BPPARAM)),
        # BiocParallel itself fails, for example, when a worker is killed
        error = function(cnd) {
            stop(
                "Parallel processing failed; a worker may have exited ",
                "unexpectedly: ",
                conditionMessage(cnd),
                call. = FALSE
            )
        }
    )
    check_worker_results(res)
}


# BiocParallel::bptry() returns an error condition for a worker that failed
# (and for the elements it did not evaluate after that failure), and NULL for a
# worker that exited without returning a result. Stop with the worker's
# original error before its result is combined with the others. A successful
# worker never returns NULL here: every caller's worker function returns a
# list or a data.table, possibly empty.
check_worker_results <- function(res) {
    failed <- which(vapply(
        res,
        inherits,
        what = c("error", "try-error"),
        FUN.VALUE = logical(1)
    ))
    if (length(failed) >= 1L) {
        # report a failure itself rather than an element left unevaluated
        unevaluated <- vapply(
            res[failed],
            inherits,
            what = "unevaluated_error",
            FUN.VALUE = logical(1)
        )
        i <- c(failed[!unevaluated], failed[unevaluated])[1]
        cause <- if (inherits(res[[i]], "try-error")) {
            attr(res[[i]], "condition")
        } else {
            res[[i]]
        }
        detail <- if (inherits(cause, "condition")) {
            conditionMessage(cause)
        } else {
            as.character(res[[i]])
        }
        stop("Parallel worker ", i, " failed: ", detail, call. = FALSE)
    }
    for (i in seq_along(res)) {
        if (is.null(res[[i]])) {
            stop(
                "Parallel worker ",
                i,
                " returned no result; it may have exited unexpectedly.",
                call. = FALSE
            )
        }
    }
    invisible(res)
}


get_log_dir <- function(log_dir) {
    log_dir <- paste0(log_dir, "/sceptre_logs/")
    if (!dir.exists(log_dir)) {
        dir.create(log_dir, recursive = TRUE)
    }
    return(log_dir)
}
