# Execution concerns shared by the response-oriented permutation plan and the
# guide-oriented CRT plan. Numerical preparation and traversal stay in those
# separate plans; these helpers only handle worker dispatch, logs and binding.
#' Resolve and optionally initialize a parallel progress file.
#' @param analysis_type,proc_id Existing role and block ID used in its filename.
#' @param parallel,print_progress Logical flags; either false returns NULL.
#' @param log_dir Existing logging root, resolved by get_log_dir.
#' @param initialize Whether to truncate/create the file; CRT guide phase
#'   appends.
#' @return File path or NULL. May create directories and a file; no RNG changes.
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

#' Announce parallel association work and its existing optional log location.
#' @param analysis_type,print_progress,log_dir Role and progress settings.
#' @return Invisible message result. Reporting may create the logging directory.
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

#' Execute ordered work blocks using the existing serial or fork dispatch.
#' @param blocks Ordered list of inputs, already partitioned by the caller.
#' @param worker Function accepting a block and NULL or one-based process ID.
#' @param parallel Logical switch; fork count equals the number of blocks.
#' @return Block results in input-block order. Retains mclapply defaults and
#'   caller
#'   worker effects; does not repartition or reseed work.
#' @noRd
execute_analysis_blocks <- function(blocks, worker, parallel) {
    if (!parallel) {
        lapply(blocks, function(block) worker(block, NULL))
    } else {
        parallel::mclapply(
            seq_along(blocks),
            function(proc_id) {
                worker(blocks[[proc_id]], proc_id)
            },
            mc.cores = length(blocks)
        )
    }
}

#' Bind worker tables with the existing schema and row order.
#' @param blocks Ordered list of compatible result tables.
#' @return data.table using rbindlist(fill=TRUE), including its factor handling.
#' @noRd
bind_analysis_results <- function(blocks) {
    data.table::rbindlist(blocks, fill = TRUE)
}


#' Report progress for one item in a worker block.
#'
#' @param position One-based item position, restarted for each block and phase.
#' @param id Character display ID; lookup belongs to the caller.
#' @param total Number of items in this block.
#' @param print_progress,parallel Logical reporting and file-output switches.
#' @param log_file Existing worker file path for parallel reporting.
#' @param feature,prefix Existing display text, preserved by each caller.
#' @return Invisible NULL. Emits at item 1 and multiples of five. Parallel write
#'   retains its extra terminating newline; disabled reporting does no file I/O.
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

#' Run periodic worker maintenance independently of progress reporting.
#'
#' @param position One-based block-local position.
#' @return Invisible NULL. Calls GC at multiples of 200, after progress and
#'   before
#'   the caller reads or fits the current item. Does not change RNG state.
#' @noRd
collect_worker_garbage <- function(position) {
    if (position %% 200L == 0L) {
        gc() |> invisible()
    }
    invisible(NULL)
}
