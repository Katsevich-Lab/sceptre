#!/usr/bin/env Rscript
root <- commandArgs(TRUE)[1L]
if (is.na(root)) stop("Usage: compare-qualification.R STAGING_ROOT")
root <- normalizePath(root, mustWork = TRUE)
.libPaths(c("/private/tmp/sceptre-rf0-ondisc-library", .libPaths()))
stopifnot(requireNamespace("jsonlite", quietly = TRUE))
read_json <- function(path) jsonlite::read_json(file.path(root, path), simplifyVector = FALSE)
read_payload <- function(side, lane, file) readRDS(file.path(root, side, lane, file))
schema_and_exact_values <- function(x, path = "root") {
    if (isS4(x)) {
        slots <- methods::slotNames(x)
        attrs <- attributes(x)[setdiff(names(attributes(x)), slots)]
        return(list(type = typeof(x), class = class(x), attributes = attrs, slots = stats::setNames(
            lapply(slots, function(n) schema_and_exact_values(methods::slot(x, n), paste0(path, "/@", n))), slots)))
    }
    if (is.list(x)) return(list(type = typeof(x), class = class(x), attributes = attributes(x),
        items = stats::setNames(lapply(seq_along(x), function(i) schema_and_exact_values(x[[i]],
            paste0(path, "/", if (is.null(names(x))) i else names(x)[i]))), names(x))))
    if (is.double(x)) {
        exact <- grepl("idx|indices|indexes|pod_map|grna_assignment_matrix", path)
        return(list(type = typeof(x), class = class(x), attributes = attributes(x),
            missing = is.na(x), nan = is.nan(x), infinite = is.infinite(x),
            infinity_values = unname(x[is.infinite(x)]), exact_values = if (exact) x else NULL))
    }
    if (typeof(x) %in% c("environment", "externalptr")) stop("Nonportable value in comparison payload: ", path)
    list(type = typeof(x), class = class(x), attributes = attributes(x), exact_values = x)
}
compare_payload <- function(lane, file) {
    a <- read_payload("reference", lane, file)
    b <- read_payload("candidate", lane, file)
    numeric <- all.equal(a, b, tolerance = 1e-8, check.attributes = TRUE)
    exact <- all.equal(a, b, tolerance = 0, check.attributes = TRUE)
    structure <- identical(schema_and_exact_values(a), schema_and_exact_values(b))
    list(scientific_values_tolerance_1e8 = isTRUE(numeric), exact_payload = isTRUE(exact),
        exact_identifiers_order_decisions_classes_schemas_missingness = structure,
        tolerance_diagnostics = if (isTRUE(numeric)) NULL else numeric,
        exact_diagnostics = if (isTRUE(exact)) NULL else head(exact, 30L))
}
stages <- lapply(c("reference", "candidate"), function(side) read_json(paste0(side, "/stages/report.json")))
names(stages) <- c("reference", "candidate")
nextflow <- lapply(c("reference", "candidate"), function(side) read_json(paste0(side, "/nextflow/nextflow-report.json")))
names(nextflow) <- c("reference", "candidate")
books <- lapply(c("reference", "candidate"), function(side) read_json(paste0(side, "/book/book-report.json")))
names(books) <- c("reference", "candidate")
stage_ok <- vapply(stages, function(x) length(x$stages) == 30L && all(vapply(x$stages, function(s) identical(s$status, "passed"), logical(1L))), logical(1L))
nf_ok <- vapply(nextflow, function(x) isTRUE(x$passed) && identical(names(x$runs), c("standard", "nuclear")) &&
    all(vapply(x$runs, function(s) identical(s$status, "passed") && s$exit_status == 0L, logical(1L))), logical(1L))
expected_bad <- "import-data.qmd:28"
book_checks <- lapply(books, function(x) {
    failed <- names(x$chunks)[!vapply(x$chunks, function(s) identical(s$status, "passed"), logical(1L))]
    barcode_files <- x$input_files[endsWith(names(x$input_files), "barcodes.tsv.gz")]
    list(passed_cells = sum(vapply(x$chunks, function(s) identical(s$status, "passed"), logical(1L))),
        failed_cells = failed, known_failure_only = identical(failed, expected_bad),
        error_matches = grepl("^File is empty:", x$chunks[[expected_bad]]$error),
        two_empty_barcode_files = length(barcode_files) == 2L && all(vapply(barcode_files,
            function(s) identical(s$bytes, 20L) && identical(s$barcode_lines, 0L), logical(1L))),
        overall_book_not_claimed_pass = identical(x$passed, FALSE), input_files = x$input_files)
})
known_book <- all(vapply(book_checks, function(x) x$passed_cells == 14L &&
    x$known_failure_only && x$error_matches && x$two_empty_barcode_files && x$overall_book_not_claimed_pass, logical(1L))) &&
    identical(book_checks$reference$input_files, book_checks$candidate$input_files)
maps <- lapply(c("reference", "candidate"), function(side) readRDS(file.path(root, side, "nextflow/nuclear/published/response_to_pod_map.rds")))
map_equal <- identical(maps[[1L]], maps[[2L]])
comparisons <- list(actual_scripts = compare_payload("stages", "comparison.rds"),
    complete_nextflow = compare_payload("nextflow", "nextflow-comparison.rds"),
    completed_uncached_book_cells = compare_payload("book", "book-comparison.rds"))
data_versions <- list(
    actual_scripts = lapply(stages, function(x) list(package_path = x$package_path, dependencies = x$dependencies,
        fixture_hashes = x$fixture_hashes, inventory = x$consumer_inventory)),
    nextflow = lapply(nextflow, function(x) list(package = x$package, ondisc = x$ondisc, runtime = x$nextflow,
        source_revision = x$source$revision, archive_sha256 = x$source$archive_sha256, inputs = x$inputs,
        task_processes = lapply(x$runs, function(s) s$task_R_processes), nuclear_rows = x$runs$nuclear$arrow_dataset_rows)),
    book = lapply(books, function(x) list(package = x$package, ondisc = x$ondisc, sceptredata = x$sceptredata, revision = x$book_revision)))
result <- list(schema_version = 1L, observed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    scope = "Final clean-O2 original21f9ea0 versus production28df390; actual scripts, local seeded Nextflow, bounded uncached book cells.",
    execution = list(actual_script_stages_passed = as.list(stage_ok), nextflow_workflows_passed = as.list(nf_ok),
        book = book_checks, known_baseline_book_data_failure_preserved = known_book,
        exact_response_to_pod_plan = map_equal),
    comparisons = comparisons, provenance = data_versions,
    scoped_qualification_passed = all(stage_ok) && all(nf_ok) && known_book && map_equal &&
        all(vapply(comparisons, function(x) x$scientific_values_tolerance_1e8 &&
            x$exact_identifiers_order_decisions_classes_schemas_missingness, logical(1L))),
    limits = c("Whole-book success is not claimed: both retain the identical known CellRanger example-data failure.",
        "Historical containers, nf-test/hosted schedulers and unmodified default array configuration remain unqualified.",
        "Local Nextflow uses documented array=null, disabled containers/conda, offline cached plugin and fixed R task seed271828."))
jsonlite::write_json(result, file.path(root, "qualification-comparison.json"), pretty = TRUE,
    auto_unbox = TRUE, null = "null", na = "null", digits = NA)
cat(jsonlite::toJSON(list(scoped_qualification_passed = result$scoped_qualification_passed,
    stage_execution = as.list(stage_ok), nextflow_execution = as.list(nf_ok),
    known_book_failure = known_book, pod_maps_exact = map_equal, comparisons = comparisons), auto_unbox = TRUE), "\n")
quit(status = if (result$scoped_qualification_passed) 0L else 1L)
