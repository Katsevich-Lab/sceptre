# Helpers for starting workers, recording progress and joining their results.
# Both permutation and CRT use these mechanics, while their gene/guide loops
# and statistical calculations remain in the separate analysis functions.
#' Choose a worker's progress-log path and optionally create the file
#' @param analysis_type,proc_id Analysis name and worker number for the filename.
#' @param parallel,print_progress Return NULL unless both are TRUE.
#' @param log_dir Directory passed to `get_log_dir()`, which creates it if needed.
#' @param initialize TRUE creates or empties the file. FALSE lets the CRT guide
#'   phase append to the file already used during response fitting.
#' @return Log-file path, or NULL when file-based reporting is disabled.
#' @noRd
analysis_log_file <- function(
    analysis_type,
    proc_id,
    parallel,
    print_progress,
    log_dir,
    initialize = TRUE
) {
    if (!parallel || !print_progress) {
        return(NULL)
    }
    path <- paste0(get_log_dir(log_dir), analysis_type, "_", proc_id, ".out")
    if (initialize) {
        file.create(path) |> invisible()
    }
    path
}

#' Print the parallel-analysis message and where to find progress logs
#' @param analysis_type,print_progress,log_dir Analysis name and reporting settings.
#' @return The invisible result of the message call. Looking up the log location
#'   can create its directory when progress reporting is enabled.
#' @noRd
announce_parallel_analysis <- function(analysis_type, print_progress, log_dir) {
    message("Running ", analysis_type, " in parallel.")
    if (print_progress) {
        message(
            "Change directories to ",
            crayon::blue(get_log_dir(log_dir)),
            " and view the files ",
            crayon::blue(paste0(analysis_type, "_*.out")),
            " for progress updates."
        )
    }
}

#' Run each work block serially or in a separate forked worker
#' @param blocks Input blocks, already divided and ordered by the caller.
#' @param worker Function accepting a block and a worker number; serial work
#'   receives NULL instead of a number.
#' @param parallel Whether to use `mclapply()`, with one worker per block.
#' @return Worker results in input-block order, including legitimate NULL results.
#'   Stop if a worker fails
#'   or exits without returning a result. Include its block number and original
#'   error message, before callers combine results or report success. Blocks and
#'   random seeds are not rearranged; `mclapply()` keeps its default behavior.
#' @noRd
execute_analysis_blocks <- function(blocks, worker, parallel) {
    # A successful worker can return NULL for an empty block. Put its result
    # inside a list so that this differs from a child exiting without replying.
    replies <- if (!parallel) {
        lapply(blocks, function(block) list(value = worker(block, NULL)))
    } else {
        parallel::mclapply(
            seq_along(blocks),
            function(proc_id) {
                list(value = worker(blocks[[proc_id]], proc_id))
            },
            mc.cores = length(blocks)
        )
    }
    for (i in seq_along(replies)) {
        result <- replies[[i]]
        if (inherits(result, "try-error")) {
            cause <- attr(result, "condition")
            detail <- if (inherits(cause, "condition")) {
                conditionMessage(cause)
            } else {
                as.character(result)
            }
            stop("Analysis worker ", i, " failed: ", detail, call. = FALSE)
        }
        if (is.null(result)) {
            stop(
                "Analysis worker ",
                i,
                " returned no result; it may have exited unexpectedly.",
                call. = FALSE
            )
        }
    }
    lapply(replies, function(reply) reply$value)
}

#' Join worker result tables in block order
#' @param blocks Ordered list of compatible result tables.
#' @return A data.table from `rbindlist(fill = TRUE)`. Rows keep block order and
#'   columns absent from a block are filled with missing values.
#' @noRd
bind_analysis_results <- function(blocks) {
    data.table::rbindlist(blocks, fill = TRUE)
}


#' Report which gene or guide group a worker is processing
#' @param position Item number within this block and phase, starting at 1.
#' @param id Gene or guide-group name to display.
#' @param total Number of items assigned to this block.
#' @param print_progress,parallel Whether to report progress and use a worker log.
#' @param log_file Worker log path for parallel runs.
#' @param feature,prefix Text supplied by the calling analysis for its message.
#' @return Invisible NULL. Report the first item and every fifth item. Parallel
#'   runs append to the log, retaining the extra newline from `write()`; serial
#'   runs print a message. Disabled reporting does not write files.
#' @noRd
report_worker_progress <- function(
    position,
    id,
    total,
    print_progress,
    parallel,
    log_file,
    feature = "response",
    prefix = "Analyzing pairs containing "
) {
    if ((position == 1L || position %% 5L == 0L) && print_progress) {
        msg <- paste0(
            prefix,
            feature,
            " ",
            id,
            " (",
            position,
            " of ",
            total,
            ")\n"
        )
        if (parallel) {
            write(x = msg, file = log_file, append = TRUE)
        } else {
            message(msg, appendLF = FALSE)
        }
    }
    invisible(NULL)
}

#' Ask R to reclaim unused memory every 200 items
#' @param position Item number within this worker's block, starting at 1.
#' @return Invisible NULL. Run garbage collection after the progress message and
#'   before reading or fitting the current item, even if progress is hidden.
#' @noRd
collect_worker_garbage <- function(position) {
    if (position %% 200L == 0L) {
        gc() |> invisible()
    }
    invisible(NULL)
}
