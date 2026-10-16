#' Traverse responses for permutation inference using a required prepared context.
#'
#' @param context Prepared store, eligible pairs, canonical index, common matrix
#'   design/view, compatible response fits and checked test settings.
#' @param execution Output, progress, parallel-worker and log-directory options.
#' @param synthetic_idxs One complete permutation bank shared across this run.
#' @return Result table and merged compact response_precomputations. Traverses
#'   first-request response order (or the existing seeded block order); retains
#'   pair order within each response. Shared fits are inserted in traversal order;
#'   NT discovery fits remain guide-specific. Reads one row per outer response,
#'   including cache hits. Progress/logging and periodic GC precede each row read.
#' @noRd
run_permutation_analysis <- function(context, execution, synthetic_idxs) {
    index <- context$index
    treatment_lookup <- make_treatment_index_lookup(index)
    response_rows <- unique(index$pair_rows$response_idx)

    #' Analyze one ordered block of source response rows.
    #' @param block_rows One-based source rows; proc_id is the log-file block ID.
    #' @return Block results and newly fitted shared responses, preserving order.
    #'   Captures immutable prepared inputs; owns only its local result/cache list.
    #' @noRd
    analyze_response_block <- function(block_rows, proc_id = NULL) {
        log_file <- analysis_log_file(context$analysis, proc_id,
            execution$parallel, execution$print_progress, execution$log_dir)
        updates <- list()
        results <- vector("list", length(block_rows))
        for (position in seq_along(block_rows)) {
            row <- block_rows[position]
            id <- as.character(index$response_ids[row])
            report_worker_progress(position, id, length(block_rows),
                execution$print_progress, execution$parallel, log_file)
            collect_worker_garbage(position)
            expression <- storage_read_response(context$response_store, row,
                                                 context$common_view$source_cells)
            request_rows <- index$rows_by_response[[as.character(row)]]
            pairs <- context$pairs[request_rows, , drop = FALSE]
            groups <- as.character(pairs$grna_group)
            if (context$shared_response_fit) {
                fit <- lookup_response_fit(context$response_fits, id)
                if (is.null(fit)) {
                    fit <- fit_response(expression, context$design)
                    updates[[id]] <- fit
                }
                score <- prepare_score_components(expression, context$design, fit)
                native_results <- test_permutation_shared_fit(expression, score,
                    treatment_lookup, groups, synthetic_idxs, context$test_options,
                    execution$output_amount)
            } else {
                native_results <- test_permutation_group_fit(expression,
                    context$design, index, groups, synthetic_idxs,
                    context$test_options, execution$output_amount)
            }
            results[[position]] <- pair_results_from_native(pairs, native_results,
                                                            execution$output_amount)
        }
        list(result = bind_analysis_results(results), response_fits = updates)
    }

    blocks <- partition_response_ids(response_rows, execution$parallel,
                                     execution$n_processors)
    if (execution$parallel) announce_parallel_analysis(context$analysis,
        execution$print_progress, execution$log_dir)
    chunks <- execute_analysis_blocks(blocks, analyze_response_block,
                                      execution$parallel)
    if (execution$parallel) message(crayon::green(" \u2713"))
    result <- lapply(chunks, function(chunk) chunk$result) |> bind_analysis_results()
    updates <- lapply(chunks, function(chunk) chunk$response_fits) |> purrr::flatten()
    list(result = result, response_precomputations =
        merge_response_fits(context$response_fits, updates))
}

#' Traverse guide groups for CRT inference using a required prepared context.
#'
#' @param context Prepared store, eligible pairs, index, common matrix/view,
#'   checked test settings and compatible saved response fits.
#' @param execution Output, progress, worker and log-directory options.
#' @return Result table and compact response_precomputations. Shared-view response
#'   fits are prepared first, including row reads on cache hits; guide-specific NT
#'   discovery fits are not cached. Each guide worker draws exactly one complete
#'   CRT bank before its response loop. Existing block, row and cache insertion
#'   order, two-phase log append behavior, progress cadence and GC are retained.
#' @noRd
run_crt_analysis <- function(context, execution) {
    index <- context$index
    treatment_lookup <- make_treatment_index_lookup(index)
    response_fits <- context$response_fits

    #' Prepare missing shared response fits for an ordered source-row block.
    #' @param block_rows One-based response rows; proc_id selects the progress log.
    #' @return Only new fits, keyed and inserted in response traversal order. Reads
    #'   each row before its single cache lookup, even when the fit is available.
    #' @noRd
    prepare_response_block <- function(block_rows, proc_id = NULL) {
        log_file <- analysis_log_file(context$analysis, proc_id,
            execution$parallel, execution$print_progress, execution$log_dir)
        updates <- list()
        for (position in seq_along(block_rows)) {
            row <- block_rows[position]
            id <- as.character(index$response_ids[row])
            report_worker_progress(position, id, length(block_rows),
                execution$print_progress, execution$parallel, log_file,
                prefix = "Running precomputation on ")
            collect_worker_garbage(position)
            expression <- storage_read_response(context$response_store, row,
                                                 context$common_view$source_cells)
            fit <- lookup_response_fit(response_fits, id)
            if (is.null(fit)) updates[[id]] <- fit_response(expression, context$design)
        }
        updates
    }

    if (context$shared_response_fit) {
        blocks <- partition_response_ids(unique(index$pair_rows$response_idx),
            execution$parallel, execution$n_processors)
        if (execution$parallel) announce_parallel_analysis(context$analysis,
            execution$print_progress, execution$log_dir)
        chunks <- execute_analysis_blocks(blocks, prepare_response_block,
                                          execution$parallel)
        response_fits <- merge_response_fits(response_fits, purrr::flatten(chunks))
    }

    #' Analyze one block of canonical guide-group indices.
    #' @param block_groups One-based group positions; proc_id selects append log.
    #' @return Bound requested results in block/group/response order. Materializes
    #'   only the current guide's NT fit view, never an all-guide view bank.
    #' @noRd
    analyze_guide_block <- function(block_groups, proc_id = NULL) {
        log_file <- analysis_log_file(context$analysis, proc_id,
            execution$parallel, execution$print_progress, execution$log_dir,
            initialize = FALSE)
        results <- vector("list", length(block_groups))
        for (position in seq_along(block_groups)) {
            group_idx <- block_groups[position]
            id <- as.character(index$group_ids[group_idx])
            report_worker_progress(position, id, length(block_groups),
                execution$print_progress, execution$parallel, log_file,
                feature = "gRNA group")
            collect_worker_garbage(position)
            request_rows <- index$rows_by_group[[as.character(group_idx)]]
            pairs <- context$pairs[request_rows, , drop = FALSE]
            response_rows <- index$pair_rows$response_idx[request_rows]
            if (context$shared_response_fit) {
                native_results <- test_crt_shared_fit(context$response_store,
                    response_rows, as.character(pairs$response_id), context$design,
                    context$common_view$source_cells, treatment_lookup(id),
                    response_fits, context$test_options, execution$output_amount)
            } else {
                view <- make_fit_view(index, id)
                native_results <- test_crt_group_fit(context$response_store,
                    response_rows, context$design, view, context$test_options,
                    execution$output_amount)
            }
            results[[position]] <- pair_results_from_native(pairs, native_results,
                                                            execution$output_amount)
        }
        bind_analysis_results(results)
    }

    blocks <- partition_response_ids(unique(index$pair_rows$group_idx),
        execution$parallel, execution$n_processors)
    chunks <- execute_analysis_blocks(blocks, analyze_guide_block,
                                      execution$parallel)
    if (execution$parallel) message(crayon::green(" \u2713"))
    list(result = bind_analysis_results(chunks),
         response_precomputations = response_fits)
}
