#!/usr/bin/env Rscript
# Real local-executor qualification, with an untouched archive of the pinned
# consumer source. Containers and nf-test remain separate qualification lanes.
args <- commandArgs(TRUE)
if (length(args) != 4L) stop("Usage: pipeline-nextflow.R LIB PIPELINE_CHECKOUT FIXTURES OUT")
library_dir <- normalizePath(args[1L], mustWork = TRUE)
pipeline <- normalizePath(args[2L], mustWork = TRUE)
fixtures <- normalizePath(args[3L], mustWork = TRUE)
dir.create(args[4L], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[4L], mustWork = TRUE)
if (length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("OUT must be empty.")
stopifnot(requireNamespace("jsonlite", quietly = TRUE), requireNamespace("digest", quietly = TRUE))
.libPaths(c(library_dir, .libPaths()))
expected_package <- normalizePath(file.path(library_dir, "sceptre"), mustWork = TRUE)
suppressPackageStartupMessages(library("sceptre", lib.loc = library_dir, character.only = TRUE))
stopifnot(identical(normalizePath(find.package("sceptre")), expected_package))
hash <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
plain <- function(x) {
    if (is.data.frame(x)) return(list(class = class(x), names = names(x), row_names = rownames(x), columns = lapply(x, plain)))
    if (is.list(x)) return(lapply(x, plain))
    x
}
report <- list(schema_version = 1L, evidence = "executed", observed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package = list(path = expected_package, version = as.character(packageVersion("sceptre"))),
    ondisc = list(path = find.package("ondisc"), version = as.character(packageVersion("ondisc"))),
    R = R.version.string, fixtures = fixtures, runs = list(),
    scope = "Unmodified pinned pipeline; actual local Nextflow tasks; explicit source package/dependency libraries; no containers or nf-test.",
    runtime_overrides = c("Local executor, one CPU/task, one concurrent task (explicit arrival-order control), containers/conda disabled.",
        "Job arrays disabled with array=null for exact archived array process selector; default closure is rejected by Nextflow25.10.4 and local does not support arrays.",
        "NXF_OFFLINE=true uses the installed runtime and cached nf-schema plugin.",
        "Each R task starts with seed271828, freezing the original unseeded response-pod sampler for a shared reference/candidate work plan.",
        "Fresh R startup profile asserts package path and records package/dependency provenance in each task."),
    outstanding_gates = c("Historical sceptre0.10.2 containers", "nf-test/hosted executor profiles", "Unmodified default array configuration"))
payload <- list(schema_version = 1L, runs = list())
save_reports <- function() {
    jsonlite::write_json(report, file.path(out, "nextflow-report.json"), pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
    saveRDS(payload, file.path(out, "nextflow-comparison.rds"), version = 3L)
}
quote_groovy <- function(x) paste0("'", gsub("'", "\\'", x, fixed = TRUE), "'")

main <- function() {
    executable <- Sys.which("nextflow")
    if (!nzchar(executable)) stop("nextflow is not on PATH.")
    report$nextflow <<- list(executable = executable, version_output = system2(executable, "-version", stdout = TRUE, stderr = TRUE,
        env = c("NXF_OFFLINE=true", "NXF_DISABLE_CHECK_LATEST=true")))
    script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
    inventory <- jsonlite::read_json(file.path(dirname(normalizePath(script_file)), "consumer-inventory.json"))
    pin <- Filter(function(x) identical(x$name, "sceptre-pipeline"), inventory$repositories)[[1L]]$revision
    archive_tar <- file.path(out, "pipeline-source.tar")
    status <- system2("git", c("-C", shQuote(pipeline), "archive", "--format=tar", "-o", shQuote(archive_tar), pin))
    if (status != 0L) stop("git archive failed.")
    archived <- file.path(out, "pipeline")
    dir.create(archived)
    utils::untar(archive_tar, exdir = archived)
    for (row in inventory$scripts) stopifnot(identical(hash(file.path(archived, row$path)), row$sha256))
    files <- list.files(archived, recursive = TRUE, full.names = TRUE)
    report$source <<- list(revision = pin, archive_sha256 = hash(archive_tar),
        file_hashes = stats::setNames(lapply(files, hash), substring(files, nchar(archived) + 2L)))
    source(file.path(fixtures, "..", "..", "helper-refactor.R"), local = environment())
    input_path <- file.path(fixtures, "backed-inputs.rds")
    case_path <- file.path(fixtures, "cases.rds")
    inputs <- readRDS(input_path)
    case <- readRDS(case_path)$high_perm
    data_dir <- file.path(out, "inputs")
    dir.create(data_dir)
    set.seed(271828)
    object <- refactor_import(inputs, case, backed = TRUE, directory = data_dir)
    object <- sceptre::set_analysis_parameters(object,
        discovery_pairs = inputs$discovery_pairs, positive_control_pairs = inputs$positive_control_pairs,
        formula_object = stats::as.formula("~ batch + z", env = baseenv()),
        side = case$side, grna_integration_strategy = case$integration,
        control_group = case$control, resampling_mechanism = case$mechanism, resampling_approximation = case$approximation)
    sceptre::write_ondisc_backed_sceptre_object(object, data_dir)
    report$inputs <<- list(backed_fixture_sha256 = hash(input_path), cases_sha256 = hash(case_path),
        response_dimensions = dim(inputs$response_matrix), grna_dimensions = dim(inputs$high_grna_matrix), integer_id = object@integer_id)
    profile <- file.path(out, "runtime.Rprofile")
    rquote <- function(x) encodeString(x, quote = '"')
    writeLines(c(paste0(".libPaths(c(", paste(vapply(.libPaths(), rquote, character(1L)), collapse = ","), "))"),
        paste0("suppressPackageStartupMessages(library(sceptre, lib.loc=", rquote(library_dir), "))"),
        paste0("stopifnot(identical(normalizePath(find.package('sceptre')), ", rquote(expected_package), "))"),
        "set.seed(271828)",
        "saveRDS(list(package_path=normalizePath(find.package('sceptre')), version=as.character(utils::packageVersion('sceptre')), ondisc_version=as.character(utils::packageVersion('ondisc')), initial_seed=271828), paste0('runtime-', Sys.getpid(), '.rds'))"), profile)
    config <- file.path(out, "local.config")
    writeLines(c(
        "process { executor = 'local'; container = null; conda = null; cache = false; maxForks = 1",
        " withName: '.*' { cpus = 1; memory = '1 GB'; time = '2m' }",
        " withName: 'ASSIGN_GRNAS|RUN_ASSOCIATION_ANALYSIS|RUN_DISCOVERY_ANALYSIS_TRANS' { array = null }",
        "}", "executor { queueSize = 1 }",
        "docker.enabled = false", "singularity.enabled = false", "apptainer.enabled = false",
        "conda.enabled = false", "podman.enabled = false",
        paste0("env { R_LIBS_USER = ", quote_groovy(library_dir)),
        paste0(" R_LIBS = ", quote_groovy(paste(.libPaths(), collapse = .Platform$path.sep))),
        paste0(" R_PROFILE_USER = ", quote_groovy(profile), " }"),
        "params { grna_assignment_method = 'thresholding'; threshold = 40",
        " n_nonzero_trt_thresh = 7; n_nonzero_cntrl_thresh = 7",
        " response_n_umis_range_lower = 0; response_n_umis_range_upper = 1",
        " response_n_nonzero_range_lower = 0; response_n_nonzero_range_upper = 1",
        " n_calibration_pairs = 4; calibration_group_size = 2",
        " grna_pod_size = 2; pair_pod_size = 2; trial = false }"), config)
    report$runtime_config <<- list(path = config, sha256 = hash(config), content = readLines(config), profile_sha256 = hash(profile))
    for (kind in c("standard", "nuclear")) {
        directory <- file.path(out, kind)
        dir.create(directory)
        published <- file.path(directory, "published")
        command <- c("run", shQuote(file.path(archived, "main.nf")), "-c", shQuote(config),
            "-work-dir", shQuote(file.path(directory, "work")),
            "--sceptre_object_fp", shQuote(file.path(data_dir, "sceptre_object.rds")),
            "--response_odm_fp", shQuote(file.path(data_dir, "response.odm")),
            "--grna_odm_fp", shQuote(file.path(data_dir, "grna.odm")),
            "--outdir", shQuote(published), "-ansi-log", "false")
        if (kind == "nuclear") command <- c(command, "--discovery_pairs", "trans")
        old <- setwd(directory)
        started <- proc.time()[["elapsed"]]
        status <- tryCatch(system2(executable, command, stdout = "stdout.log", stderr = "stderr.log",
            env = c("NXF_OFFLINE=true", "NXF_DISABLE_CHECK_LATEST=true")), error = function(e) {
                writeLines(conditionMessage(e), "launch-error.log"); 127L
            })
        setwd(old)
        row <- list(evidence = "executed", exit_status = status, elapsed_seconds = proc.time()[["elapsed"]] - started,
            directory = directory, command = command, logs = c("stdout.log", "stderr.log", ".nextflow.log"))
        failure <- tryCatch({
            if (status != 0L) stop("Nextflow exited with status ", status, "; see captured logs.")
            runtime_files <- list.files(file.path(directory, "work"), pattern = "^runtime-.*[.]rds$", recursive = TRUE, full.names = TRUE)
            provenance <- lapply(runtime_files, readRDS)
            stopifnot(length(provenance) > 0L,
                all(vapply(provenance, function(x) identical(x$package_path, expected_package), logical(1L))),
                all(vapply(provenance, function(x) identical(x$ondisc_version, report$ondisc$version), logical(1L))))
            row$task_R_processes <- length(provenance)
            row$task_provenance <- provenance
            if (kind == "standard") {
                results <- lapply(c("calibration_check", "power_check", "discovery_analysis"), function(role) readRDS(file.path(published, paste0("results_run_", role, ".rds"))))
                names(results) <- c("calibration", "power", "discovery")
                stopifnot("significant" %in% names(results$discovery), nrow(results$discovery) == nrow(inputs$discovery_pairs),
                    identical(results$discovery$significant, stats::p.adjust(results$discovery$p_value, method = object@multiple_testing_method) < object@multiple_testing_alpha))
                payload$runs[[kind]] <<- lapply(results, plain)
                row$result_rows <- lapply(results, nrow)
            } else {
                shard_files <- list.files(file.path(published, "trans_results"), pattern = "[.]parquet$", full.names = TRUE)
                stopifnot(length(shard_files) > 0L)
                shards <- lapply(shard_files, function(path) as.data.frame(arrow::read_parquet(path)))
                stopifnot(all(vapply(shards, function(x) !"significant" %in% names(x), logical(1L))))
                payload$runs[[kind]] <<- stats::setNames(lapply(shards, plain), basename(shard_files))
                row$shards <- lapply(shards, function(x) list(rows = nrow(x), columns = names(x), classes = lapply(x, class)))
                dataset <- arrow::open_dataset(file.path(published, "trans_results"))
                collected <- dplyr::collect(dplyr::mutate(dataset, response_id = as.character(response_id), grna_target = as.character(grna_target)))
                stopifnot(nrow(collected) == sum(vapply(shards, nrow, integer(1L))))
                row$arrow_dataset_rows <- nrow(collected)
            }
            NULL
        }, error = function(e) conditionMessage(e))
        row$status <- if (is.null(failure)) "passed" else "failed"
        row$error <- failure
        report$runs[[kind]] <<- row
        save_reports()
    }
}
fatal <- tryCatch({ main(); NULL }, error = function(e) conditionMessage(e))
report$fatal_error <- fatal
report$passed <- is.null(fatal) && length(report$runs) == 2L && all(vapply(report$runs, function(x) identical(x$status, "passed"), logical(1L)))
save_reports()
cat(jsonlite::toJSON(list(passed = report$passed, fatal_error = fatal,
    runs = lapply(report$runs, function(x) list(status = x$status, error = x$error))), auto_unbox = TRUE, null = "null"), "\n")
quit(status = if (report$passed) 0L else 1L)
