#!/usr/bin/env Rscript
# Execute real downstream scripts against an explicitly installed sceptre build.
# Usage: Rscript tools/refactor/pipeline-smoke.R LIB PIPELINE FIXTURES OUT
# Run baseline/candidate into separate directories; compare comparison.rds with
# all.equal(..., tolerance = 0). report.json keeps provenance and outstanding gates.

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 4L) stop("Usage: pipeline-smoke.R <sceptre-library-dir> <pipeline-checkout> <fixture-dir> <output-dir>")
library_dir <- normalizePath(argv[1L], mustWork = TRUE)
pipeline_dir <- normalizePath(argv[2L], mustWork = TRUE)
fixture_dir <- normalizePath(argv[3L], mustWork = TRUE)
dir.create(argv[4L], recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(argv[4L], mustWork = TRUE)
if (length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Use a fresh empty output directory to preserve previous runs.")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required for the smoke report.")
.libPaths(c(library_dir, .libPaths()))
expected_package <- normalizePath(file.path(library_dir, "sceptre"), mustWork = TRUE)
if ("sceptre" %in% loadedNamespaces() &&
    normalizePath(getNamespaceInfo(asNamespace("sceptre"), "path")) != expected_package) {
    stop("A different sceptre build is already loaded; start a fresh R process.")
}
suppressPackageStartupMessages(library("sceptre", lib.loc = library_dir, character.only = TRUE))
stopifnot(identical(normalizePath(find.package("sceptre")), expected_package))

report <- list(schema_version = 1L, observed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    evidence = "executed", package_path = expected_package,
    package_version = as.character(packageVersion("sceptre")), R = R.version.string,
    library_paths = .libPaths(), pipeline_checkout = pipeline_dir,
    fixture_dir = fixture_dir, stages = list(), comparisons = list(),
    outstanding_gates = c("Full Nextflow/nf-test execution", "Historical container execution"),
    limitations = c("Only the backed high-MOI permutation fixture is exercised here.",
        "Nuclear threshold CLI values are built but not passed by the pinned script.",
        "The pinned nuclear writer writes Parquet even when use_parquet=FALSE; this harness uses TRUE.",
        "The zero-passing nuclear branch preserves baseline logical NA p-value/effect columns."))
report$dependencies <- stats::setNames(lapply(c("ondisc", "arrow", "Matrix", "data.table", "dplyr", "ggplot2"), function(package) {
    if (!requireNamespace(package, quietly = TRUE)) return(list(available = FALSE))
    list(available = TRUE, version = as.character(packageVersion(package)),
        path = normalizePath(find.package(package)))
}), c("ondisc", "arrow", "Matrix", "data.table", "dplyr", "ggplot2"))
comparison <- list(schema_version = 1L, stages = list())
write_reports <- function() {
    jsonlite::write_json(report, file.path(output_dir, "report.json"), pretty = TRUE,
        auto_unbox = TRUE, null = "null", na = "null", digits = NA)
    saveRDS(comparison, file.path(output_dir, "comparison.rds"), version = 3L)
}
git_value <- function(args) {
    value <- suppressWarnings(system2("git", c("-C", shQuote(pipeline_dir), args), stdout = TRUE, stderr = FALSE))
    if (!is.null(attr(value, "status"))) stop("Cannot read pipeline Git revision.")
    paste(value, collapse = "\n")
}
hash_file <- function(path) {
    # digest is available with the test tooling and hashes raw bytes, not an R serialization.
    if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required for source hashes.")
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
plain <- function(x) {
    if (is.data.frame(x)) return(list(class = class(x), names = names(x),
        row_names = rownames(x), columns = lapply(x, plain)))
    if (is.list(x)) return(lapply(x, plain))
    x
}
table_summary <- function(x) {
    numerical_names <- names(x)[vapply(x, is.numeric, logical(1L))]
    list(class = class(x), rows = nrow(x), columns = names(x),
        column_classes = lapply(x, class),
        numerical = stats::setNames(lapply(numerical_names, function(name) as.vector(x[[name]])), numerical_names),
        significant_present = "significant" %in% names(x))
}

# Rscript -e preserves these trailing arguments for the actual script's own
# commandArgs(). The source is unmodified; tracing only records its real loader's
# return value, after its backing-file identity checks have succeeded.
bootstrap <- paste(
    '.libPaths(c(Sys.getenv("SCEPTRE_SMOKE_LIBRARY"), .libPaths()))',
    'suppressPackageStartupMessages(library("sceptre", lib.loc=Sys.getenv("SCEPTRE_SMOKE_LIBRARY")))',
    'expected <- normalizePath(file.path(Sys.getenv("SCEPTRE_SMOKE_LIBRARY"), "sceptre"))',
    'stopifnot(identical(normalizePath(find.package("sceptre")), expected))',
    'saveRDS(list(pid=Sys.getpid(), package_path=expected, version=as.character(packageVersion("sceptre")), library_paths=.libPaths(), ondisc_version=if(requireNamespace("ondisc", quietly=TRUE)) as.character(packageVersion("ondisc")) else NA_character_, args=commandArgs(TRUE)), "child-metadata.rds")',
    'trace("read_ondisc_backed_sceptre_object", where=asNamespace("sceptre"), print=FALSE, exit=quote({',
    'o <- returnValue(); saveRDS(list(object_id=o@integer_id, response_id=sceptre::get_response_matrix(o)@integer_id, grna_id=sceptre::get_grna_matrix(o)@integer_id, cells_in_use=o@cells_in_use), "reattachment.rds")',
    '}))',
    'set.seed(271828)',
    'smoke_warnings <- character()',
    'tryCatch(withCallingHandlers(source(Sys.getenv("SCEPTRE_SMOKE_SCRIPT"), local=.GlobalEnv, echo=FALSE), warning=function(w) { smoke_warnings <<- c(smoke_warnings, conditionMessage(w)) }), finally=saveRDS(smoke_warnings, "warnings.rds"))',
    sep = ";\n")

main <- function() {
    report$pipeline_revision <<- git_value("rev-parse HEAD")
    inventory_path <- file.path(dirname(fixture_dir), "..", "..", "..", "tools", "refactor", "consumer-inventory.json")
    # Normal repo location is tools/refactor; alternatively the caller can point
    # at a fixture tree whose repository has the checked-in inventory.
    script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
    local_inventory <- file.path(dirname(normalizePath(script_arg)), "consumer-inventory.json")
    if (file.exists(local_inventory)) inventory_path <- local_inventory
    if (!file.exists(inventory_path)) stop("Cannot locate tools/refactor/consumer-inventory.json beside the harness or fixture repository.")
    inventory <- jsonlite::read_json(inventory_path, simplifyVector = FALSE)
    pinned <- Filter(function(x) identical(x$name, "sceptre-pipeline"), inventory$repositories)[[1L]]$revision
    if (!identical(report$pipeline_revision, pinned)) stop("Pipeline HEAD differs from the consumer inventory revision.")
    script_rows <- inventory$scripts
    report$consumer_inventory <<- list(path = normalizePath(inventory_path), sha256 = hash_file(inventory_path))
    for (row in script_rows) {
        if (!identical(hash_file(file.path(pipeline_dir, row$path)), row$sha256)) {
            stop("Pipeline script differs from inspected baseline: ", row$path)
        }
    }
    helper <- normalizePath(file.path(fixture_dir, "..", "..", "helper-refactor.R"), mustWork = TRUE)
    source(helper, local = environment())
    inputs_path <- file.path(fixture_dir, "backed-inputs.rds")
    cases_path <- file.path(fixture_dir, "cases.rds")
    inputs <- readRDS(inputs_path)
    case <- readRDS(cases_path)$high_perm
    stopifnot(identical(case$mechanism, "permutations"), ncol(inputs$response_matrix) > 1000L)
    report$fixture_hashes <<- list(backed_inputs = hash_file(inputs_path), cases = hash_file(cases_path), helper = hash_file(helper))
    report$fixture_dimensions <<- list(response = dim(inputs$response_matrix), grna = dim(inputs$high_grna_matrix))

    run_stage <- function(label, script, args, expected, setup = NULL, validate = NULL) {
        directory <- file.path(output_dir, label)
        dir.create(directory)
        dir.create(file.path(directory, "tmp"))
        if (!is.null(setup)) setup(directory)
        full_script <- file.path(pipeline_dir, "bin", script)
        old <- setwd(directory)
        started <- proc.time()[["elapsed"]]
        status <- tryCatch(system2(file.path(R.home("bin"), "Rscript"),
            c("--vanilla", "-e", shQuote(bootstrap), vapply(args, shQuote, character(1L))),
            stdout = "stdout.log", stderr = "stderr.log",
            env = c(paste0("R_LIBS_USER=", shQuote(library_dir)),
                paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep))),
                paste0("SCEPTRE_SMOKE_LIBRARY=", shQuote(library_dir)),
                paste0("SCEPTRE_SMOKE_SCRIPT=", shQuote(full_script)),
                paste0("TMPDIR=", shQuote(file.path(directory, "tmp"))))),
            error = function(e) { writeLines(conditionMessage(e), "launch-error.log"); 127L })
        setwd(old)
        present <- file.exists(file.path(directory, expected))
        meta_path <- file.path(directory, "child-metadata.rds")
        attach_path <- file.path(directory, "reattachment.rds")
        meta <- if (file.exists(meta_path)) readRDS(meta_path) else NULL
        attached <- if (file.exists(attach_path)) readRDS(attach_path) else NULL
        checks <- list(child_package_path = !is.null(meta) && identical(meta$package_path, expected_package),
            reattached = !is.null(attached),
            backing_ids_match = !is.null(attached) && identical(attached$object_id, attached$response_id) && identical(attached$object_id, attached$grna_id),
            expected_outputs = all(present))
        ok <- identical(as.integer(status), 0L) && all(unlist(checks))
        report$stages[[label]] <<- list(status = if (ok) "passed" else "failed", evidence = "executed",
            script = paste0("bin/", script), script_sha256 = hash_file(full_script), args = args,
            child_exit_status = as.integer(status), elapsed_seconds = proc.time()[["elapsed"]] - started,
            directory = directory, checks = checks, child = meta, reattachment = attached,
            artifacts = stats::setNames(as.list(present), expected), logs = c("stdout.log", "stderr.log"))
        warning_path <- file.path(directory, "warnings.rds")
        report$stages[[label]]$warnings <<- if (file.exists(warning_path)) readRDS(warning_path) else character()
        comparison$stages[[label]] <<- list(status = if (ok) "passed" else "failed")
        inspection_error <- tryCatch({ if (ok) {
            for (name in expected) {
                path <- file.path(directory, name)
                if (grepl("[.]rds$", name)) {
                    value <- readRDS(path)
                    if (is.data.frame(value)) {
                        comparison$stages[[label]][[name]] <<- plain(value)
                        report$stages[[label]]$result_summary[[name]] <<- table_summary(value)
                    } else if (isS4(value) && is(value, "sceptre_object")) {
                        slots <- c("functs_called", "last_function_called", "response_precomputations", "calibration_result", "power_result", "discovery_result", "nuclear", "nf_pipeline", "n_discovery_pairs", "grna_assignment_method", "initial_grna_assignment_list", "grna_assignments_raw", "grna_assignments", "cells_in_use", "M_matrix", "discovery_pairs")
                        comparison$stages[[label]][[name]] <<- stats::setNames(lapply(slots, function(n) plain(methods::slot(value, n))), slots)
                    } else comparison$stages[[label]][[name]] <<- plain(value)
                } else if (grepl("[.]parquet$", name)) {
                    value <- as.data.frame(arrow::read_parquet(path))
                    comparison$stages[[label]][[name]] <<- plain(value)
                    report$stages[[label]]$result_summary[[name]] <<- table_summary(value)
                    shard_checks <- list(no_significant = !"significant" %in% names(value),
                        result_columns = all(c("response_id", "grna_target", "p_value", "log_2_fold_change") %in% names(value)),
                        p_values = (is.numeric(value$p_value) || (is.logical(value$p_value) && all(is.na(value$p_value)))) &&
                            all(is.na(value$p_value) | (value$p_value >= 0 & value$p_value <= 1)))
                    report$stages[[label]]$contract_checks <<- shard_checks
                    if (!all(unlist(shard_checks))) {
                        report$stages[[label]]$status <<- "failed"
                        report$stages[[label]]$contract_failure <<- "Nuclear raw-shard contract failed."
                        ok <- FALSE
                    }
                }
            }
        }
        if (ok && !is.null(validate)) {
            contract_checks <- validate(directory)
            report$stages[[label]]$contract_checks <<- contract_checks
            ok <- all(unlist(contract_checks))
        }
        NULL }, error = function(e) conditionMessage(e))
        if (!is.null(inspection_error)) {
            report$stages[[label]]$inspection_error <<- inspection_error
            ok <- FALSE
        }
        report$stages[[label]]$status <<- if (ok) "passed" else "failed"
        comparison$stages[[label]]$status <<- if (ok) "passed" else "failed"
        write_reports()
        if (ok) directory else NULL
    }
    blocked <- function(label, reason) {
        report$stages[[label]] <<- list(status = "blocked", evidence = "not executed", reason = reason)
        comparison$stages[[label]] <<- list(status = "blocked")
        write_reports()
        NULL
    }

    # Each family starts from independently prepared real backed inputs. Fixed
    # seeds retain backing identity across baseline/candidate runs; no pointer is
    # used in a child without the original script reattaching its backing files.
    prepare_family <- function(name, nuclear = FALSE) {
        directory <- file.path(output_dir, paste0(name, "-inputs"))
        dir.create(directory)
        set.seed(271828)
        if (!nuclear) object <- refactor_prepare(inputs, case, backed = TRUE, directory = directory)
        else {
            object <- refactor_import(inputs, case, backed = TRUE, directory = directory)
            object <- sceptre::set_analysis_parameters(object,
                formula_object = stats::as.formula("~ batch + z", env = baseenv()),
                side = case$side, grna_integration_strategy = case$integration,
                control_group = case$control, resampling_mechanism = case$mechanism,
                resampling_approximation = case$approximation, multiple_testing_method = "none")
            object@nuclear <- TRUE
            object <- sceptre::assign_grnas(object, method = "thresholding", threshold = 40, print_progress = FALSE)
        }
        sceptre::write_ondisc_backed_sceptre_object(object, directory)
        list(object = file.path(directory, "sceptre_object.rds"), response = file.path(directory, "response.odm"), grna = file.path(directory, "grna.odm"))
    }
    capture_prepare <- function(name, nuclear = FALSE) {
        warnings <- character()
        value <- withCallingHandlers(prepare_family(name, nuclear), warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
        })
        report$input_preparation[[name]] <<- list(evidence = "executed", warnings = warnings)
        value
    }
    # Exercise the original upstream CLI consumers, rather than replacing them
    # with helper calls. Helpers only create the unassigned backed input object.
    upstream <- function(name, selected_case, maximum = FALSE) {
        root <- file.path(output_dir, paste0(name, "-inputs"))
        dir.create(root)
        set.seed(271828)
        object <- refactor_import(inputs, selected_case, backed = TRUE, directory = root)
        object <- sceptre::set_analysis_parameters(object,
            discovery_pairs = inputs$discovery_pairs, positive_control_pairs = inputs$positive_control_pairs,
            formula_object = stats::as.formula("~ batch + z", env = baseenv()),
            side = selected_case$side, grna_integration_strategy = selected_case$integration,
            control_group = selected_case$control, resampling_mechanism = selected_case$mechanism,
            resampling_approximation = selected_case$approximation)
        sceptre::write_ondisc_backed_sceptre_object(object, root)
        backing <- file.path(root, c("response.odm", "grna.odm"))
        parameters <- run_stage(paste0(name, "-parameters"), "set_analysis_parameters.R",
            c(file.path(root, "sceptre_object.rds"), backing, rep("default", 7L), rep("null.rds", 3L), "FALSE"),
            c("sceptre_object.rds", "analysis_summary.txt"), setup = function(directory) saveRDS(NULL, file.path(directory, "null.rds")),
            validate = function(directory) {
                value <- readRDS(file.path(directory, "sceptre_object.rds"))
                list(normal_mode = identical(value@nuclear, FALSE), pairs_retained = identical(value@discovery_pairs, object@discovery_pairs))
            })
        if (is.null(parameters)) return(NULL)
        metadata <- file.path(parameters, "sceptre_object.rds")
        plan <- run_stage(paste0(name, "-guide-plan"), "output_grna_info.R", c(metadata, backing, "2", "FALSE"),
            c("grna_pods.txt", "grna_to_pod_map.rds", "low_moi.txt"), validate = function(directory) {
                map <- readRDS(file.path(directory, "grna_to_pod_map.rds"))
                list(unique_guides = !anyDuplicated(map$grna_id), known_guides = all(map$grna_id %in% inputs$grna_target_data_frame$grna_id),
                    matching_pods = identical(scan(file.path(directory, "grna_pods.txt"), quiet = TRUE), as.numeric(unique(map$pod_id))))
            })
        if (is.null(plan)) return(NULL)
        formula <- file.path(parameters, "null.rds")
        map_path <- file.path(plan, "grna_to_pod_map.rds")
        map <- readRDS(map_path)
        workers <- if (maximum) list() else lapply(sort(unique(map$pod_id)), function(pod) {
            run_stage(paste0(name, "-assignment-pod-", pod), "assign_grnas.R",
                c(metadata, backing, map_path, as.character(pod), "thresholding", "40", rep("default", 4L), formula),
                c("grna_assignments.rds", "grna_assignment_formula.rds"))
        })
        if (any(vapply(workers, is.null, logical(1L)))) return(NULL)
        files <- if (maximum) "dummy.rds" else paste0("grna_assignments_", seq_along(workers), ".rds")
        gathered <- run_stage(paste0(name, "-assignment-gather"), "process_grna_assignments.R",
            c(metadata, backing, if (maximum) "maximum" else "thresholding", "default", "default", formula, files),
            c("sceptre_object.rds", "grna_assignment_matrix.rds", "plot_grna_count_distributions.png", "plot_assign_grnas.png", "analysis_summary.txt"),
            setup = function(directory) {
                if (maximum) saveRDS(NULL, file.path(directory, "dummy.rds"))
                else for (i in seq_along(workers)) stopifnot(file.copy(file.path(workers[[i]], "grna_assignments.rds"), file.path(directory, files[i])))
            }, validate = function(directory) {
                value <- readRDS(file.path(directory, "sceptre_object.rds"))
                list(assignment_completed = isTRUE(value@functs_called[["assign_grnas"]]),
                    method_preserved = identical(value@grna_assignment_method, if (maximum) "maximum" else "thresholding"),
                    assignment_dimensions = identical(dim(readRDS(file.path(directory, "grna_assignment_matrix.rds"))), dim(inputs$high_grna_matrix)))
            })
        if (!is.null(gathered)) run_stage(paste0(name, "-qc"), "run_qc.R",
            c(file.path(gathered, "sceptre_object.rds"), backing, "7", "7", "0", "1", "0", "1", "default"),
            c("sceptre_object.rds", "plot_covariates.png", "plot_run_qc.png", "analysis_summary.txt"))
        # The trans sentinel and its normal parameter path are distinct consumers.
        run_stage(paste0(name, "-trans-parameters"), "set_analysis_parameters.R",
            c(file.path(root, "sceptre_object.rds"), backing, rep("default", 7L), "null.rds", "trans.rds", "null.rds", "FALSE"),
            c("sceptre_object.rds", "analysis_summary.txt"), setup = function(directory) {
                saveRDS(NULL, file.path(directory, "null.rds")); saveRDS("trans", file.path(directory, "trans.rds"))
            }, validate = function(directory) {
                value <- readRDS(file.path(directory, "sceptre_object.rds"))
                list(nuclear = isTRUE(value@nuclear), unadjusted = identical(value@multiple_testing_method, "none"),
                    no_precomputed_pairs = nrow(value@discovery_pairs) == 0L && nrow(value@positive_control_pairs) == 0L)
            })
        invisible(NULL)
    }
    tryCatch(upstream("upstream-high", case), error = function(e) blocked("upstream-high", conditionMessage(e)))
    tryCatch(upstream("upstream-low", readRDS(cases_path)$low_complement_perm, maximum = TRUE),
        error = function(e) blocked("upstream-low", conditionMessage(e)))
    standard <- tryCatch(capture_prepare("standard"), error = function(e) {
        blocked("standard-inputs", conditionMessage(e)); NULL
    })
    if (!is.null(standard)) {
        prepared <- run_stage("standard-prepare", "prepare_association_analyses.R",
            c(standard$object, standard$response, standard$grna, "4", "2", "2"),
            c("sceptre_object.rds", "run_calibration_check", "run_power_check", "run_discovery_analysis"))
        current <- if (!is.null(prepared)) file.path(prepared, "sceptre_object.rds") else NULL
        for (kind in c("calibration_check", "power_check", "discovery_analysis")) {
            if (is.null(current)) { blocked(paste0("standard-", kind), "A preceding preparation/gather stage failed."); next }
            object <- readRDS(current)
            pair_slot <- switch(kind, calibration_check = "negative_control_pairs", power_check = "positive_control_pairs_with_info", discovery_analysis = "discovery_pairs_with_info")
            pods <- sort(unique(methods::slot(object, pair_slot)$pod))
            workers <- lapply(pods, function(pod) run_stage(paste0("standard-", kind, "-pod-", pod),
                "run_association_analysis.R", c(current, standard$response, standard$grna, as.character(pod), paste0("run_", kind)),
                c("result.rds", "precomputations.rds")))
            if (!length(pods) || any(vapply(workers, is.null, logical(1L)))) {
                blocked(paste0("standard-", kind, "-gather"), "No pods or at least one worker failed."); current <- NULL; next
            }
            results <- paste0("results_", seq_along(workers), ".rds")
            precomputations <- paste0("precomputations_", seq_along(workers), ".rds")
            gather <- run_stage(paste0("standard-", kind, "-gather"), "process_association_analysis_results.R",
                c("input.rds", "response.odm", "grna.odm", paste0("run_", kind), results, precomputations),
                c("sceptre_object.rds", paste0("results_run_", kind, ".rds"), paste0("plot_run_", kind, ".png"), "analysis_summary.txt"),
                setup = function(directory) {
                    stopifnot(file.copy(c(current, standard$response, standard$grna), file.path(directory, c("input.rds", "response.odm", "grna.odm"))))
                    for (i in seq_along(workers)) stopifnot(file.copy(
                        file.path(workers[[i]], c("result.rds", "precomputations.rds")),
                        file.path(directory, c(results[i], precomputations[i]))))
                }, validate = function(directory) {
                    result <- readRDS(file.path(directory, paste0("results_run_", kind, ".rds")))
                    saved <- readRDS(file.path(directory, "sceptre_object.rds"))
                    input <- readRDS(file.path(directory, "input.rds"))
                    checks <- list(pod_removed = !"pod" %in% names(result),
                        all_pairs_retained = nrow(result) == nrow(methods::slot(input, pair_slot)))
                    if (kind != "power_check") checks$global_significance <- identical(result$significant,
                        stats::p.adjust(result$p_value, method = input@multiple_testing_method) < input@multiple_testing_alpha)
                    if (kind == "discovery_analysis") checks$pruned_discovery_object <- is.null(saved)
                    if (kind == "calibration_check") checks$pruned_calibration_columns <- identical(
                        names(saved@calibration_result), c("p_value", "log_2_fold_change", "significant"))
                    checks
                })
            current <- if (!is.null(gather)) file.path(gather, "sceptre_object.rds") else NULL
        }
    }
    nuclear <- tryCatch(capture_prepare("nuclear", nuclear = TRUE), error = function(e) {
        blocked("nuclear-inputs", conditionMessage(e)); NULL
    })
    if (!is.null(nuclear)) {
        qc <- run_stage("nuclear-qc", "run_qc_nuclear.R",
            c(nuclear$object, nuclear$response, nuclear$grna, "0", "1", "0", "1", "default"),
            c("sceptre_object.rds", "plot_covariates.png", "plot_cellwise_qc.png", "analysis_summary.txt"))
        response_map <- data.frame(response_id = rownames(inputs$response_matrix), pod_id = rep(1:2, length.out = nrow(inputs$response_matrix)))
        map_path <- file.path(output_dir, "nuclear-response-map.rds")
        saveRDS(response_map, map_path, version = 3L)
        report$nuclear_response_map <<- list(path = map_path, sha256 = hash_file(map_path), rows = response_map, evidence = "supplied deterministic map")
        if (!requireNamespace("arrow", quietly = TRUE)) {
            blocked("nuclear-shards", "Arrow is unavailable; raw Parquet gate remains outstanding.")
        } else if (is.null(qc)) blocked("nuclear-shards", "Nuclear cellwise QC failed.")
        else {
            run_stage("nuclear-response-plan", "prepare_association_analysis_nuclear.R",
                c(file.path(qc, "sceptre_object.rds"), nuclear$response, nuclear$grna, "8"),
                c("discovery_analysis_pods", "response_to_pod_map.rds"), validate = function(directory) {
                    plan <- readRDS(file.path(directory, "response_to_pod_map.rds"))
                    list(all_responses = setequal(plan$response_id, rownames(inputs$response_matrix)),
                        one_pod_per_response = !anyDuplicated(plan$response_id),
                        matching_pods = identical(scan(file.path(directory, "discovery_analysis_pods"), quiet = TRUE), as.numeric(unique(plan$pod_id))))
                })
            for (trial in c(TRUE, FALSE)) for (pod in sort(unique(response_map$pod_id))) {
                label <- if (trial) paste0("nuclear-pod-", pod) else paste0("nuclear-nontrial-pod-", pod)
                run_stage(label, "run_association_analysis_nuclear.R",
                    c(file.path(qc, "sceptre_object.rds"), nuclear$response, nuclear$grna, map_path,
                        as.character(pod), "default", "default", as.character(trial), "TRUE"), paste0("result_", pod, ".parquet"))
            }
            zero_map <- file.path(output_dir, "nuclear-zero-response-map.rds")
            zero_ids <- rownames(inputs$response_matrix)[rowSums(inputs$response_matrix) == 0]
            stopifnot(length(zero_ids) > 0L)
            saveRDS(data.frame(response_id = zero_ids, pod_id = 3L), zero_map, version = 3L)
            run_stage("nuclear-zero-pass", "run_association_analysis_nuclear.R",
                c(file.path(qc, "sceptre_object.rds"), nuclear$response, nuclear$grna, zero_map,
                    "3", "default", "default", "FALSE", "TRUE"), "result_3.parquet", validate = function(directory) {
                        value <- as.data.frame(arrow::read_parquet(file.path(directory, "result_3.parquet")))
                        list(nonempty_failed_family = nrow(value) > 0L && !any(value$pass_qc),
                            only_empty_responses = all(as.character(value$response_id) %in% zero_ids),
                            missing_statistics = all(is.na(value$p_value)) && all(is.na(value$log_2_fold_change)),
                            baseline_na_schema = is.logical(value$p_value) && is.logical(value$log_2_fold_change),
                            no_significant = !"significant" %in% names(value))
                    })
        }
    }
}

fatal <- tryCatch({ main(); NULL }, error = function(e) conditionMessage(e))
if (!is.null(fatal)) report$fatal_error <- fatal
status <- vapply(report$stages, function(x) x$status, character(1L))
report$summary <- list(passed = sum(status == "passed"), failed = sum(status == "failed"),
    blocked = sum(status == "blocked"), fatal = fatal,
    comparison_payload = "comparison.rds",
    comparison_command = "all.equal(readRDS(BASELINE/comparison.rds), readRDS(CANDIDATE/comparison.rds), tolerance=0)")
write_reports()
cat(jsonlite::toJSON(report$summary, auto_unbox = TRUE, null = "null"), "\n")
quit(status = if (is.null(fatal) && length(status) && all(status == "passed")) 0L else 1L)
