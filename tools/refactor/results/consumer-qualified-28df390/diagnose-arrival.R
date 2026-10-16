#!/usr/bin/env Rscript
root <- normalizePath(commandArgs(TRUE)[1L], mustWork = TRUE)
plain <- function(x) {
    if (is.data.frame(x)) return(list(class = class(x), names = names(x), row_names = rownames(x), columns = lapply(x, plain)))
    if (is.list(x)) return(lapply(x, plain))
    x
}
inspect_arrival <- function(side, lane) {
    work <- file.path(root, side, lane, "standard/work")
    commands <- list.files(work, pattern = "^[.]command[.]sh$", all.files = TRUE, recursive = TRUE, full.names = TRUE)
    roles <- list()
    raw <- list()
    for (script in commands) {
        source <- readLines(script)
        if (!any(grepl("process_association_analysis_results.R", source, fixed = TRUE))) next
        role <- c("calibration_check", "power_check", "discovery_analysis")
        role <- role[vapply(role, function(x) any(grepl(paste0("run_", x), source, fixed = TRUE)), logical(1L))]
        stopifnot(length(role) == 1L)
        files <- list.files(dirname(script), pattern = "^results[0-9]*$", full.names = TRUE)
        rows <- lapply(files, function(path) {
            value <- readRDS(path)
            pods <- unique(value$pod)
            stopifnot(length(pods) == 1L)
            raw[[paste0(role, "/pod-", pods)]] <<- plain(value)
            list(staged_name = basename(path), worker_source = normalizePath(path), pod = pods,
                response_levels = levels(value$response_id), group_levels = levels(value$grna_group),
                response_ids = as.character(value$response_id), group_ids = as.character(value$grna_group))
        })
        result <- readRDS(file.path(dirname(script), paste0("results_run_", role, ".rds")))
        roles[[role]] <- list(staged_arrival = rows, final_response_levels = levels(result$response_id),
            final_target_levels = levels(result$grna_target), final_rows = as.data.frame(result))
    }
    list(roles = roles, raw = raw)
}
parallel <- lapply(c("reference", "candidate"), function(x) inspect_arrival(x, "nextflow"))
names(parallel) <- c("reference", "candidate")
keys <- names(parallel$reference$raw)
raw_equal <- stats::setNames(lapply(keys, function(k) all.equal(parallel$reference$raw[[k]], parallel$candidate$raw[[k]], tolerance = 0)), keys)
result <- list(schema_version = 1L, evidence = "inspected executed workflow artifacts", diagnostic_only = TRUE,
    original_parallel_comparison_remains_failed = TRUE, no_published_output_reordered_or_modified = TRUE,
    cause = "Nextflow collection stages results1/results2 in completion order. Their factor response levels are merged in that order; the original gather sorts tied p-values by factor codes, so final levels and tied rows depend on arrival.",
    parallel = lapply(parallel, `[[`, "roles"), raw_worker_payloads_compared_by_existing_pod_id = raw_equal)
serialized_files <- file.path(root, c("reference", "candidate"), "nextflow-serialized/nextflow-report.json")
if (all(file.exists(serialized_files))) {
    serial <- lapply(c("reference", "candidate"), function(x) inspect_arrival(x, "nextflow-serialized"))
    names(serial) <- c("reference", "candidate")
    result$serialized <- lapply(serial, `[[`, "roles")
    arrival_plan <- function(x) {
        roles <- c("calibration_check", "power_check", "discovery_analysis")
        stats::setNames(lapply(roles, function(name) vapply(x$roles[[name]]$staged_arrival,
            function(row) row$pod, integer(1L))), roles)
    }
    result$serialized_arrival_plans_exact <- identical(arrival_plan(serial$reference), arrival_plan(serial$candidate))
}
controlled_files <- file.path(root, c("reference", "candidate"), "nextflow-arrival-controlled/nextflow-report.json")
if (all(file.exists(controlled_files))) {
    controlled <- lapply(c("reference", "candidate"), function(x) inspect_arrival(x, "nextflow-arrival-controlled"))
    names(controlled) <- c("reference", "candidate")
    result$arrival_controlled <- lapply(controlled, `[[`, "roles")
    result$controlled_arrival_plans_exact <- identical(arrival_plan(controlled$reference), arrival_plan(controlled$candidate))
}
jsonlite::write_json(result, file.path(root, "arrival-order-diagnostic.json"), pretty = TRUE, auto_unbox = TRUE,
    null = "null", na = "null", digits = NA)
cat(jsonlite::toJSON(list(raw_worker_pods_exact = vapply(raw_equal, isTRUE, logical(1L)),
    serialized_arrival_plans_exact = result$serialized_arrival_plans_exact,
    controlled_arrival_plans_exact = result$controlled_arrival_plans_exact), auto_unbox = TRUE), "\n")
