#' Run permutation tests one response gene at a time
#' @param context Data and settings from `prepare_analysis_context()`: selected
#'   cells, gene-guide pairs that passed QC, aligned design and reusable fits.
#' @param execution Output level, progress, parallel-worker and log settings.
#' @param synthetic_idxs Complete set of permutation draws, shared across this run.
#' @return Result table and the old plus newly computed reusable response fits.
#' @details Read each gene once per outer loop, even if its fit is saved. Reuse its
#'   fit and score arrays across guide groups when they use the same cells.
#'   NT-control discovery needs a separate treatment-plus-NT fit for each group.
#'   Visit genes in first-request order, or the existing seeded block order when
#'   partitioned. Keep requests within each gene in order. Progress reporting and
#'   periodic garbage collection happen before each gene is read.
#' @noRd
run_permutation_analysis <- function(context, execution, synthetic_idxs) {
    index <- context$index
    treatment_lookup <- make_treatment_index_lookup(index)
    response_rows <- if (is.null(execution$work_keys)) {
        unique(index$pair_rows$response_idx)
    } else {
        execution$work_keys
    }

    # Analyze a block of response genes assigned to one worker
    # block_rows: Original count-matrix gene rows, in the order to process.
    # proc_id: Worker number used in its log filename, or NULL for serial work.
    # Returns: This block's results and newly computed shared response fits.
    #   Each worker builds its own output lists and reads the prepared inputs.
    analyze_response_block <- function(block_rows, proc_id = NULL) {
        log_file <- analysis_log_file(
            context$analysis,
            proc_id,
            execution$parallel,
            execution$print_progress,
            execution$log_dir
        )
        updates <- list()
        results <- vector("list", length(block_rows))
        for (position in seq_along(block_rows)) {
            row <- block_rows[position]
            id <- as.character(index$response_ids[row])
            report_worker_progress(
                position,
                id,
                length(block_rows),
                execution$print_progress,
                execution$parallel,
                log_file
            )
            collect_worker_garbage(position)
            expression <- storage_read_response(
                context$response_store,
                row,
                context$common_view$source_cells
            )
            request_rows <- index$rows_by_response[[as.character(row)]]
            pairs <- context$pairs[request_rows, , drop = FALSE]
            groups <- as.character(pairs$grna_group)
            if (context$shared_response_fit) {
                fit <- lookup_response_fit(context$response_fits, id)
                if (is.null(fit)) {
                    fit <- fit_response(expression, context$design)
                    updates[[id]] <- fit
                }
                score <- prepare_score_components(
                    expression,
                    context$design,
                    fit
                )
                native_results <- test_permutation_shared_fit(
                    expression,
                    score,
                    treatment_lookup,
                    groups,
                    synthetic_idxs,
                    context$test_options,
                    execution$output_amount
                )
            } else {
                native_results <- test_permutation_group_fit(
                    expression,
                    context$design,
                    index,
                    groups,
                    synthetic_idxs,
                    context$test_options,
                    execution$output_amount
                )
            }
            results[[position]] <- pair_results_from_native(
                pairs,
                native_results,
                execution$output_amount
            )
        }
        list(result = bind_analysis_results(results), response_fits = updates)
    }

    blocks <- partition_response_ids(
        response_rows,
        execution$parallel,
        execution$n_processors
    )
    if (execution$parallel) {
        announce_parallel_analysis(
            context$analysis,
            execution$print_progress,
            execution$log_dir
        )
    }
    chunks <- execute_analysis_blocks(
        blocks,
        analyze_response_block,
        execution$parallel,
        analysis = context$analysis,
        phase = "response_tests"
    )
    if (execution$parallel) {
        message(crayon::green(" \u2713"))
    }
    result <- lapply(chunks, function(chunk) chunk$result) |>
        bind_analysis_results()
    updates <- lapply(chunks, function(chunk) chunk$response_fits) |>
        purrr::flatten()
    list(
        result = result,
        response_precomputations = merge_response_fits(
            context$response_fits,
            updates
        )
    )
}

#' Run CRT tests one guide group at a time
#' @param context Data and settings from `prepare_analysis_context()`.
#' @param execution Output level, progress, parallel-worker and log settings.
#' @return Result table and the old plus newly computed reusable response fits.
#' @details First prepare response fits that can be shared across groups, reading
#'   each gene before looking up its saved fit. Then, for each guide group, fit
#'   its assignment probabilities and generate B1+B2+B3 draws once, reusing those
#'   draws across the group's requested genes. With NT-control discovery, response
#'   fits use group-specific treatment-plus-NT cells and are not saved for reuse.
#'   Preserve block/gene order, saved-fit order and progress timing. The guide
#'   phase appends to the log started during response fitting.
#' @noRd
run_crt_analysis <- function(context, execution) {
    index <- context$index
    treatment_lookup <- make_treatment_index_lookup(index)
    response_fits <- context$response_fits

    # Prepare missing shared response fits for one worker's genes
    # block_rows: Original count-matrix gene rows, in processing order.
    # proc_id: Worker number used in the progress log, or NULL for serial work.
    # Returns: Newly fitted genes only, named and ordered as processed. Each gene is
    #   read before looking up its saved fit, including when that fit already exists.
    prepare_response_block <- function(block_rows, proc_id = NULL) {
        log_file <- analysis_log_file(
            context$analysis,
            proc_id,
            execution$parallel,
            execution$print_progress,
            execution$log_dir
        )
        updates <- list()
        for (position in seq_along(block_rows)) {
            row <- block_rows[position]
            id <- as.character(index$response_ids[row])
            report_worker_progress(
                position,
                id,
                length(block_rows),
                execution$print_progress,
                execution$parallel,
                log_file,
                prefix = "Running precomputation on "
            )
            collect_worker_garbage(position)
            expression <- storage_read_response(
                context$response_store,
                row,
                context$common_view$source_cells
            )
            fit <- lookup_response_fit(response_fits, id)
            if (is.null(fit)) {
                updates[[id]] <- fit_response(expression, context$design)
            }
        }
        updates
    }

    work_groups <- if (is.null(execution$work_keys)) {
        unique(index$pair_rows$group_idx)
    } else {
        execution$work_keys
    }
    if (context$shared_response_fit) {
        response_rows <- unique(index$pair_rows$response_idx[
            index$pair_rows$group_idx %in% work_groups
        ])
        blocks <- partition_response_ids(
            response_rows,
            execution$parallel,
            execution$n_processors
        )
        if (execution$parallel) {
            announce_parallel_analysis(
                context$analysis,
                execution$print_progress,
                execution$log_dir
            )
        }
        chunks <- execute_analysis_blocks(
            blocks,
            prepare_response_block,
            execution$parallel,
            analysis = context$analysis,
            phase = "response_fits"
        )
        response_fits <- merge_response_fits(
            response_fits,
            purrr::flatten(chunks)
        )
    }

    # Analyze a block of guide groups assigned to one worker
    # block_groups: Positions in the prepared guide-group list, starting at 1.
    # proc_id: Worker number used to append to its log, or NULL for serial work.
    # Returns: Results in block, group and within-group request order. For NT controls,
    #   construct only the current group's treatment-plus-NT selection as needed.
    analyze_guide_block <- function(block_groups, proc_id = NULL) {
        log_file <- analysis_log_file(
            context$analysis,
            proc_id,
            execution$parallel,
            execution$print_progress,
            execution$log_dir,
            initialize = FALSE
        )
        results <- vector("list", length(block_groups))
        perturbation_updates <- list()
        for (position in seq_along(block_groups)) {
            group_idx <- block_groups[position]
            id <- as.character(index$group_ids[group_idx])
            report_worker_progress(
                position,
                id,
                length(block_groups),
                execution$print_progress,
                execution$parallel,
                log_file,
                feature = "gRNA group"
            )
            collect_worker_garbage(position)
            request_rows <- index$rows_by_group[[as.character(group_idx)]]
            pairs <- context$pairs[request_rows, , drop = FALSE]
            response_rows <- index$pair_rows$response_idx[request_rows]
            treatment <- if (context$shared_response_fit) {
                treatment_lookup(id)
            } else {
                NULL
            }
            view <- if (!context$shared_response_fit) {
                make_fit_view(index, id)
            } else {
                NULL
            }
            probabilities <- NULL
            if (!index$calibration_check) {
                fit_design <- if (context$shared_response_fit) {
                    context$design
                } else {
                    context$design[view$fit_to_analysis, , drop = FALSE]
                }
                fit_cells <- if (context$shared_response_fit) {
                    context$common_view$source_cells
                } else {
                    view$source_cells
                }
                trt <- if (context$shared_response_fit) {
                    treatment$trt_idxs
                } else {
                    view$trt_fit_pos
                }
                fit_identity <- perturbation_fit_identity(
                    fit_cells,
                    fit_design,
                    trt
                )
                saved <- lookup_perturbation_fit(
                    context$perturbation_fits,
                    id,
                    fit_identity
                )
                if (is.null(saved)) {
                    proposed <- fit_perturbation_model(trt, fit_design)
                    probabilities <- proposed$probabilities
                    perturbation_updates[id] <- list(
                        if (is.null(proposed$fit)) {
                            NULL
                        } else {
                            c(proposed$fit, list(identity = fit_identity))
                        }
                    )
                } else {
                    probabilities <- reconstruct_perturbation_probabilities(
                        saved,
                        fit_design
                    )
                }
            }
            if (context$shared_response_fit) {
                native_results <- test_crt_shared_fit(
                    context$response_store,
                    response_rows,
                    as.character(pairs$response_id),
                    context$design,
                    context$common_view$source_cells,
                    treatment,
                    response_fits,
                    context$test_options,
                    execution$output_amount,
                    probabilities = probabilities
                )
            } else {
                native_results <- test_crt_group_fit(
                    context$response_store,
                    response_rows,
                    context$design,
                    view,
                    context$test_options,
                    execution$output_amount,
                    probabilities = probabilities,
                    group_design = fit_design
                )
            }
            results[[position]] <- pair_results_from_native(
                pairs,
                native_results,
                execution$output_amount
            )
            if (!index$calibration_check) rm(fit_design)
        }
        list(
            result = bind_analysis_results(results),
            perturbation_updates = perturbation_updates
        )
    }

    blocks <- partition_response_ids(
        work_groups,
        execution$parallel,
        execution$n_processors
    )
    chunks <- execute_analysis_blocks(
        blocks,
        analyze_guide_block,
        execution$parallel,
        analysis = context$analysis,
        phase = "guide_tests"
    )
    if (execution$parallel) {
        message(crayon::green(" \u2713"))
    }
    perturbation_fits <- if (index$calibration_check) {
        NULL
    } else {
        updates <- do.call(
            c,
            unname(lapply(chunks, function(chunk) chunk$perturbation_updates))
        )
        cache <- merge_perturbation_fits(context$perturbation_fits, updates)
        attr(cache, "fit_identity") <- context$perturbation_identity
        cache
    }
    list(
        result = bind_analysis_results(lapply(chunks, function(chunk) {
            chunk$result
        })),
        response_precomputations = response_fits,
        perturbation_fits = perturbation_fits
    )
}
