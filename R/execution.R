# Execution concerns shared by the response-oriented permutation plan and the
# guide-oriented CRT plan. Numerical preparation and traversal stay in those
# separate plans; these helpers only handle worker dispatch, logs and binding.
analysis_log_file <- function(analysis_type, proc_id, parallel, print_progress,
                              log_dir, initialize = TRUE) {
    if (!parallel || !print_progress) return(NULL)
    path <- paste0(get_log_dir(log_dir), analysis_type, "_", proc_id, ".out")
    if (initialize) file.create(path) |> invisible()
    path
}

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

execute_analysis_blocks <- function(blocks, worker, parallel) {
    if (!parallel) {
        lapply(blocks, function(block) worker(block, NULL))
    } else {
        parallel::mclapply(seq_along(blocks), function(proc_id) {
            worker(blocks[[proc_id]], proc_id)
        }, mc.cores = length(blocks))
    }
}

bind_analysis_results <- function(blocks) data.table::rbindlist(blocks, fill = TRUE)
