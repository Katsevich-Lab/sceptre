# helper function
get_id_from_idx <- function(
    response_idx,
    print_progress,
    response_ids,
    print_multiple = 5L,
    gc_multiple = 200L,
    feature = "response",
    str = "Analyzing pairs containing ",
    parallel = FALSE,
    f_name = NULL
) {
    if (
        (response_idx == 1 || response_idx %% print_multiple == 0) &&
            print_progress
    ) {
        msg <- paste0(
            str,
            feature,
            " ",
            as.character(response_ids[response_idx]),
            " (",
            response_idx,
            " of ",
            length(response_ids),
            ")\n"
        )
        if (parallel) {
            write(x = msg, file = f_name, append = TRUE)
        } else {
            message(msg, appendLF = FALSE)
        }
    }
    if (response_idx %% gc_multiple == 0) {
        gc() |> invisible()
    }
    response_id <- as.character(response_ids[response_idx])
    return(response_id)
}


# core function 1. run permutation test in memory
run_perm_test_in_memory <- function(
    response_matrix,
    grna_assignments,
    covariate_matrix,
    response_grna_group_pairs,
    synthetic_idxs,
    output_amount,
    resampling_approximation,
    B1,
    B2,
    B3,
    calibration_check,
    control_group_complement,
    n_nonzero_trt_thresh,
    n_nonzero_cntrl_thresh,
    side_code,
    low_moi,
    response_precomputations,
    cells_in_use,
    print_progress,
    parallel,
    n_processors,
    log_dir,
    analysis_type
) {
    # 0. define several variables
    subset_to_nt_cells <- calibration_check && !control_group_complement
    run_outer_regression <- calibration_check || control_group_complement
    all_nt_idxs <- if (!is.null(grna_assignments$all_nt_idxs)) {
        grna_assignments$all_nt_idxs
    } else {
        NA
    }
    grna_group_idxs <- grna_assignments$grna_group_idxs
    index <- prepare_index_context(cells_in_use, ncol(response_matrix),
        grna_assignments, response_grna_group_pairs, rownames(response_matrix),
        calibration_check, control_group_complement, low_moi)
    get_idx_f <- prepared_index_factory(index)
    common_view <- make_fit_view(index)
    response_ids <- unique(response_grna_group_pairs$response_id)
    fit_parametric_curve <- (resampling_approximation == "skew_normal")

    # 1. subset covariate matrix to cells_in_use and then to nt cells (if
    # applicable)
    covariate_matrix <- covariate_matrix[common_view$source_cells, , drop = FALSE]
    # Preserve the legacy NT-only one-column behavior until its separately
    # recorded compatibility decision; the index itself never drops dimensions.
    if (subset_to_nt_cells && ncol(covariate_matrix) == 1L) {
        covariate_matrix <- covariate_matrix[, 1L]
    }

    # 2. define function to loop subset of response IDs
    analyze_given_response_ids <- function(curr_response_ids, proc_id = NULL) {
        if (parallel && print_progress) {
            f_name <- paste0(
                get_log_dir(log_dir),
                analysis_type,
                "_",
                proc_id,
                ".out"
            )
            file.create(f_name) |> invisible()
        } else {
            f_name <- NULL
        }
        precomp_out_list <- list()
        result_out_list <- vector(
            mode = "list",
            length = length(curr_response_ids)
        )

        curr_response_idxs <- match(curr_response_ids, index$response_ids)
        for (response_idx in seq_along(curr_response_ids)) {
            response_id <- get_id_from_idx(
                response_idx,
                print_progress,
                curr_response_ids,
                parallel = parallel,
                f_name = f_name
            )

            # 3. load the expressions of the current response
            store_response_idx <- curr_response_idxs[response_idx]
            expression_vector <- storage_read_response(response_matrix,
                store_response_idx, common_view$source_cells)

            # 4. obtain the gRNA groups to analyze
            request_rows <- index$rows_by_response[[as.character(store_response_idx)]]
            curr_df <- response_grna_group_pairs[request_rows, ]
            grna_groups <- as.character(curr_df$grna_group)

            # 6. perform outer regression (if applicable)
            if (run_outer_regression) {
                # if precomp is available, load
                if (!is.null(lookup_response_fit(response_precomputations, response_id))) {
                    response_precomp <- lookup_response_fit(response_precomputations, response_id)
                } else {
                    # perform the regression to get the coefficients
                    response_precomp <- fit_response_current(
                        expressions = expression_vector,
                        covariate_matrix = covariate_matrix
                    )
                    # save precomputation
                    precomp_out_list[[response_id]] <- response_precomp
                }
                pieces_precomp <- prepare_score_current(
                    expression_vector,
                    covariate_matrix,
                    response_precomp,
                    full_test_stat = TRUE
                )
                curr_response_result <- perm_test_glm_factored_out(
                    synthetic_idxs,
                    B1,
                    B2,
                    B3,
                    fit_parametric_curve,
                    output_amount,
                    grna_groups,
                    expression_vector,
                    pieces_precomp,
                    get_idx_f,
                    side_code
                )
            } else {
                curr_response_result <- discovery_ntcells_perm_test(
                    synthetic_idxs,
                    B1,
                    B2,
                    B3,
                    fit_parametric_curve,
                    output_amount,
                    covariate_matrix,
                    all_nt_idxs,
                    grna_group_idxs,
                    grna_groups,
                    expression_vector,
                    side_code,
                    index = index
                )
            }

            # 9. combine the response-wise results into a data table; insert
            # into list
            result_out_list[[response_idx]] <- construct_data_frame_v2(
                curr_df,
                curr_response_result,
                output_amount
            )
        }

        # prepare output
        ret_pass_qc <- data.table::rbindlist(result_out_list, fill = TRUE)
        return(list(
            ret_pass_qc = ret_pass_qc,
            precomp_out_list = precomp_out_list
        ))
    }

    # 3. partition the response IDs
    partitioned_response_ids <- partition_response_ids(
        response_ids = response_ids,
        parallel = parallel,
        n_processors = n_processors
    )

    # 4. run the analysis
    if (!parallel) {
        res <- lapply(partitioned_response_ids, analyze_given_response_ids)
    } else {
        message("Running ", analysis_type, " in parallel.")
        log_file_pattern <- paste0(analysis_type, "_*.out")
        if (print_progress) {
            message(
                "Change directories to ",
                crayon::blue(get_log_dir(log_dir)),
                " and view the files ",
                crayon::blue(log_file_pattern),
                " for progress updates."
            )
        }
        res <- parallel::mclapply(
            seq_along(partitioned_response_ids),
            function(proc_id) {
                analyze_given_response_ids(
                    partitioned_response_ids[[proc_id]],
                    proc_id
                )
            },
            mc.cores = length(partitioned_response_ids)
        )
        message(crayon::green(" \u2713"))
    }

    # 5. combine and sort result
    result <- lapply(res, function(chunk) chunk[["ret_pass_qc"]]) |>
        data.table::rbindlist(fill = TRUE)
    precomp_out_list <- lapply(res, function(chunk) {
        chunk[["precomp_out_list"]]
    }) |>
        purrr::flatten()
    response_precomputations <- merge_response_fits(response_precomputations, precomp_out_list)

    return(list(
        result = result,
        response_precomputations = response_precomputations
    ))
}


# core function 2: run crt in memory
run_crt_in_memory_v2 <- function(
    response_matrix,
    grna_assignments,
    covariate_matrix,
    response_grna_group_pairs,
    output_amount,
    resampling_approximation,
    B1,
    B2,
    B3,
    calibration_check,
    control_group_complement,
    n_nonzero_trt_thresh,
    n_nonzero_cntrl_thresh,
    side_code,
    low_moi,
    response_precomputations,
    cells_in_use,
    print_progress,
    parallel,
    n_processors,
    log_dir,
    analysis_type
) {
    # 0. define several variables
    subset_to_nt_cells <- calibration_check && !control_group_complement
    run_outer_regression <- calibration_check || control_group_complement
    all_nt_idxs <- if (!is.null(grna_assignments$all_nt_idxs)) {
        grna_assignments$all_nt_idxs
    } else {
        NA
    }
    grna_group_idxs <- grna_assignments$grna_group_idxs
    index <- prepare_index_context(cells_in_use, ncol(response_matrix),
        grna_assignments, response_grna_group_pairs, rownames(response_matrix),
        calibration_check, control_group_complement, low_moi)
    get_idx_f <- prepared_index_factory(index)
    common_view <- make_fit_view(index)
    response_ids <- unique(response_grna_group_pairs$response_id)
    fit_parametric_curve <- (resampling_approximation == "skew_normal")

    # 1. subset covariate matrix to cells_in_use and then to nt cells (if
    # applicable)
    covariate_matrix <- covariate_matrix[common_view$source_cells, , drop = FALSE]
    # Preserve the legacy NT-only one-column behavior until its separately
    # recorded compatibility decision; the index itself never drops dimensions.
    if (subset_to_nt_cells && ncol(covariate_matrix) == 1L) {
        covariate_matrix <- covariate_matrix[, 1L]
    }

    # 2. define the run precomputation function
    run_precomp_on_given_responses <- function(
        curr_response_ids,
        proc_id = NULL
    ) {
        if (parallel && print_progress) {
            f_name <- paste0(
                get_log_dir(log_dir),
                analysis_type,
                "_",
                proc_id,
                ".out"
            )
            file.create(f_name) |> invisible()
        } else {
            f_name <- NULL
        }
        precomp_out_list <- list()

        curr_response_idxs <- match(curr_response_ids, index$response_ids)
        for (response_idx in seq_along(curr_response_ids)) {
            response_id <- get_id_from_idx(
                response_idx,
                print_progress,
                curr_response_ids,
                str = "Running precomputation on ",
                parallel = parallel,
                f_name = f_name
            )

            # 3. load the expressions of the current response
            expression_vector <- storage_read_response(response_matrix,
                curr_response_idxs[response_idx], common_view$source_cells)

            # 4. # if precomp is available, load
            if (!is.null(lookup_response_fit(response_precomputations, response_id))) {
                response_precomp <- lookup_response_fit(response_precomputations, response_id)
            } else {
                # perform the regression to get the coefficients
                response_precomp <- fit_response_current(
                    expressions = expression_vector,
                    covariate_matrix = covariate_matrix
                )
                # save precomputation
                precomp_out_list[[response_id]] <- response_precomp
            }
        }
        return(precomp_out_list)
    }

    # 5. update the response precomputations (if applicable)
    if (run_outer_regression) {
        partitioned_response_ids <- partition_response_ids(
            response_ids = response_ids,
            parallel = parallel,
            n_processors = n_processors
        )
        if (!parallel) {
            res <- lapply(
                partitioned_response_ids,
                run_precomp_on_given_responses
            )
        } else {
            message("Running ", analysis_type, " in parallel.")
            log_file_pattern <- paste0(analysis_type, "_*.out")
            if (print_progress) {
                message(
                    "Change directories to ",
                    crayon::blue(get_log_dir(log_dir)),
                    " and view the files ",
                    crayon::blue(log_file_pattern),
                    " for progress updates."
                )
            }
            res <- parallel::mclapply(
                seq_along(partitioned_response_ids),
                function(proc_id) {
                    run_precomp_on_given_responses(
                        partitioned_response_ids[[proc_id]],
                        proc_id
                    )
                },
                mc.cores = length(partitioned_response_ids)
            )
        }
        precomp_out_list <- purrr::flatten(res)
        response_precomputations <- merge_response_fits(response_precomputations, precomp_out_list)
    }

    # 6. define the function to analyze the pairs
    perform_association_analysis <- function(curr_grna_groups, proc_id) {
        result_out_list <- vector(
            mode = "list",
            length = length(curr_grna_groups)
        )
        f_name <- if (parallel && print_progress) {
            paste0(get_log_dir(log_dir), analysis_type, "_", proc_id, ".out")
        } else {
            NULL
        }

        for (grna_group_idx in seq_along(curr_grna_groups)) {
            curr_grna_group <- get_id_from_idx(
                grna_group_idx,
                print_progress,
                curr_grna_groups,
                feature = "gRNA group",
                parallel = parallel,
                f_name = f_name
            )

            # 7. obtain the genes to analyze
            group_idx <- match(curr_grna_group, index$group_ids)
            request_rows <- index$rows_by_group[[as.character(group_idx)]]
            curr_df <- response_grna_group_pairs[request_rows, ]
            response_ids <- as.character(curr_df$response_id)
            response_idxs <- index$pair_rows$response_idx[request_rows]

            # 8. call the low-level analysis function
            if (run_outer_regression) {
                curr_response_result <- crt_glm_factored_out(
                    B1,
                    B2,
                    B3,
                    fit_parametric_curve,
                    output_amount,
                    response_ids,
                    response_precomputations,
                    covariate_matrix,
                    get_idx_f,
                    curr_grna_group,
                    subset_to_nt_cells,
                    all_nt_idxs,
                    response_matrix,
                    side_code,
                    cells_in_use,
                    response_idxs = response_idxs
                )
            } else {
                curr_response_result <- discovery_ntcells_crt(
                    B1,
                    B2,
                    B3,
                    fit_parametric_curve,
                    output_amount,
                    get_idx_f,
                    response_ids,
                    covariate_matrix,
                    curr_grna_group,
                    all_nt_idxs,
                    response_matrix,
                    side_code,
                    cells_in_use,
                    response_idxs = response_idxs,
                    index = index
                )
            }
            result_out_list[[grna_group_idx]] <- construct_data_frame_v2(
                curr_df,
                curr_response_result,
                output_amount
            )
        }
        # 9. prepare output
        ret_pass_qc <- data.table::rbindlist(result_out_list, fill = TRUE)
    }

    # 10. perform the association analyses
    grna_groups <- unique(response_grna_group_pairs$grna_group)
    partitioned_grna_group_ids <- partition_response_ids(
        response_ids = grna_groups,
        parallel = parallel,
        n_processors = n_processors
    )
    if (!parallel) {
        res <- lapply(partitioned_grna_group_ids, perform_association_analysis)
    } else {
        res <- parallel::mclapply(
            seq_along(partitioned_grna_group_ids),
            function(proc_id) {
                perform_association_analysis(
                    partitioned_grna_group_ids[[proc_id]],
                    proc_id
                )
            },
            mc.cores = length(partitioned_grna_group_ids)
        )
        message(crayon::green(" \u2713"))
    }
    result <- res |> data.table::rbindlist(fill = TRUE)

    return(list(
        result = result,
        response_precomputations = response_precomputations
    ))
}
