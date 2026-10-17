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
    response_rows <- unique(index$pair_rows$response_idx)

    #' Analyze a block of response genes assigned to one worker
    #' @param block_rows Original count-matrix gene rows, in the order to process.
    #' @param proc_id Worker number used in its log filename, or NULL for serial work.
    #' @return This block's results and newly computed shared response fits.
    #'   Each worker builds its own output lists and reads the prepared inputs.
    #' @noRd
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
        execution$parallel
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

    #' Prepare missing shared response fits for one worker's genes
    #' @param block_rows Original count-matrix gene rows, in processing order.
    #' @param proc_id Worker number used in the progress log, or NULL for serial work.
    #' @return Newly fitted genes only, named and ordered as processed. Each gene is
    #'   read before looking up its saved fit, including when that fit already exists.
    #' @noRd
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

    if (context$shared_response_fit) {
        blocks <- partition_response_ids(
            unique(index$pair_rows$response_idx),
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
            execution$parallel
        )
        response_fits <- merge_response_fits(
            response_fits,
            purrr::flatten(chunks)
        )
    }

    #' Analyze a block of guide groups assigned to one worker
    #' @param block_groups Positions in the prepared guide-group list, starting at 1.
    #' @param proc_id Worker number used to append to its log, or NULL for serial work.
    #' @return Results in block, group and within-group request order. For NT controls,
    #'   construct only the current group's treatment-plus-NT selection as needed.
    #' @noRd
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
            if (context$shared_response_fit) {
                native_results <- test_crt_shared_fit(
                    context$response_store,
                    response_rows,
                    as.character(pairs$response_id),
                    context$design,
                    context$common_view$source_cells,
                    treatment_lookup(id),
                    response_fits,
                    context$test_options,
                    execution$output_amount
                )
            } else {
                view <- make_fit_view(index, id)
                native_results <- test_crt_group_fit(
                    context$response_store,
                    response_rows,
                    context$design,
                    view,
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
        bind_analysis_results(results)
    }

    blocks <- partition_response_ids(
        unique(index$pair_rows$group_idx),
        execution$parallel,
        execution$n_processors
    )
    chunks <- execute_analysis_blocks(
        blocks,
        analyze_guide_block,
        execution$parallel
    )
    if (execution$parallel) {
        message(crayon::green(" \u2713"))
    }
    list(
        result = bind_analysis_results(chunks),
        response_precomputations = response_fits
    )
}
