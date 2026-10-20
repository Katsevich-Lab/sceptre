# Helpers for starting workers, recording progress and joining their results.
# Both permutation and CRT use these mechanics, while their gene/guide loops
# and statistical calculations remain in the separate analysis functions.
#' Choose a worker's progress-log path and optionally create the file
#' @param analysis_type,proc_id Analysis name and worker number for the
#' filename.
#' @param parallel,print_progress Return NULL unless both are TRUE.
#' @param log_dir Directory passed to `get_log_dir()`, which creates it if
#' needed.
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
#' @param analysis_type,print_progress,log_dir Analysis name and reporting
#' settings.
#' @return The invisible result of the message call. Looking up the log location
#'   can create its directory when progress reporting is enabled.
#' @noRd
announce_parallel_analysis <- function(analysis_type, print_progress, log_dir) {
    message("Running ", analysis_type, " in parallel.")
    if (print_progress) {
        log_pattern <- paste0(analysis_type, "_*.out")
        message(
            "Change directories to ",
            crayon::blue(get_log_dir(log_dir)),
            " and view the files ",
            crayon::blue(log_pattern),
            " for progress updates."
        )
    }
}

# Base conditions retain the original cause and useful execution context.
analysis_worker_error <- function(
    block,
    cause = NULL,
    analysis = NULL,
    phase = NULL,
    reason = "error",
    detail = NULL
) {
    if (is.null(detail) && inherits(cause, "condition")) {
        detail <- conditionMessage(cause)
    }
    message <- if (identical(reason, "lost_reply")) {
        paste0(
            "Analysis worker ",
            block,
            " returned no result; it may have exited unexpectedly."
        )
    } else {
        paste0("Analysis worker ", block, " failed: ", detail)
    }
    structure(
        list(
            message = message,
            call = NULL,
            block = as.integer(block),
            analysis = analysis,
            phase = phase,
            reason = reason,
            cause = cause,
            parent = cause
        ),
        class = c("sceptre_worker_error", "error", "condition")
    )
}

invoke_analysis_worker <- function(
    blocks,
    worker,
    i,
    proc_id,
    parallel,
    analysis,
    phase
) {
    # Errors acquire execution context. A child interrupt is transported to
    # the parent and re-signalled unchanged; a parent interrupt escapes.
    tryCatch(
        list(value = worker(blocks[[i]], proc_id)),
        error = function(cause) {
            stop(analysis_worker_error(i, cause, analysis, phase))
        },
        interrupt = function(cause) {
            if (!parallel) {
                stop(cause)
            }
            list(interrupt = cause)
        }
    )
}

validate_analysis_replies <- function(replies, analysis, phase) {
    for (i in seq_along(replies)) {
        result <- replies[[i]]
        if (inherits(result, "try-error")) {
            cause <- attr(result, "condition")
            if (inherits(cause, "sceptre_worker_error")) {
                stop(cause)
            }
            detail <- if (inherits(cause, "condition")) {
                conditionMessage(cause)
            } else {
                as.character(result)
            }
            stop(analysis_worker_error(
                i,
                cause,
                analysis,
                phase,
                detail = detail
            ))
        }
        if (is.null(result)) {
            stop(analysis_worker_error(
                i,
                analysis = analysis,
                phase = phase,
                reason = "lost_reply"
            ))
        }
        if (inherits(result$interrupt, "interrupt")) stop(result$interrupt)
    }
    invisible(NULL)
}

#' Run each work block serially or in a separate forked worker
#' @param blocks Input blocks, already divided and ordered by the caller.
#' @param worker Function accepting a block and a worker number; serial work
#'   receives NULL instead of a number.
#' @param parallel Whether to use `mclapply()`, with one worker per block.
#' @param analysis,phase Optional scientific analysis and execution phase
#' labels.
#' @return Worker results in input-block order, including legitimate NULL
#' results.
#'   Stop if a worker fails
#'   or exits without returning a result. Include its block number and original
#'   error message, before callers combine results or report success. Blocks and
#'   random seeds are not rearranged; `mclapply()` keeps its default behavior.
#' @noRd
execute_analysis_blocks <- function(
    blocks,
    worker,
    parallel,
    analysis = NULL,
    phase = NULL
) {
    # A successful worker can return NULL for an empty block. Put its result
    # inside a list so that this differs from a child exiting without replying.
    invoke_worker <- function(i, proc_id) {
        invoke_analysis_worker(
            blocks,
            worker,
            i,
            proc_id,
            parallel,
            analysis,
            phase
        )
    }
    replies <- if (!parallel) {
        lapply(seq_along(blocks), function(i) invoke_worker(i, NULL))
    } else {
        parallel::mclapply(
            seq_along(blocks),
            function(proc_id) {
                invoke_worker(proc_id, proc_id)
            },
            mc.cores = length(blocks)
        )
    }
    validate_analysis_replies(replies, analysis, phase)
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
#' @param print_progress,parallel Whether to report progress and use a worker
#' log.
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
