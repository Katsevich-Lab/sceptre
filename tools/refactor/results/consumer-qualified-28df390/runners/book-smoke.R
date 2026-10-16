#!/usr/bin/env Rscript
# Execute a bounded whitelist of unchanged book R cells, with no knitr cache.
# Usage: Rscript tools/refactor/book-smoke.R LIB BOOK_CHECKOUT OUT
args <- commandArgs(TRUE)
child <- length(args) > 0L && identical(args[1L], "--child")
if (child) args <- args[-1L]
if (length(args) != 3L) stop("Usage: book-smoke.R LIB BOOK_CHECKOUT OUT")
library_dir <- normalizePath(args[1L], mustWork = TRUE)
book <- normalizePath(args[2L], mustWork = TRUE)
dir.create(args[3L], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[3L], mustWork = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
script <- normalizePath(script, mustWork = TRUE)
if (!child) {
    if (length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("OUT must be empty.")
    dir.create(file.path(out, "tmp"))
    status <- system2(file.path(R.home("bin"), "Rscript"),
        c("--vanilla", shQuote(script), "--child", shQuote(library_dir), shQuote(book), shQuote(out)),
        env = c(paste0("R_LIBS_USER=", shQuote(library_dir)),
            paste0("R_LIBS=", shQuote(paste(c(library_dir, .libPaths()), collapse = .Platform$path.sep))),
            paste0("TMPDIR=", shQuote(file.path(out, "tmp")))),
        stdout = file.path(out, "stdout.log"), stderr = file.path(out, "stderr.log"))
    quit(status = status)
}
setwd(out)
.libPaths(c(library_dir, .libPaths()))
suppressPackageStartupMessages(library("sceptre", lib.loc = library_dir, character.only = TRUE))
expected <- normalizePath(file.path(library_dir, "sceptre"), mustWork = TRUE)
stopifnot(identical(normalizePath(find.package("sceptre")), expected),
    requireNamespace("jsonlite", quietly = TRUE), requireNamespace("digest", quietly = TRUE))
hash <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
report <- list(schema_version = 1L, evidence = "executed", cache_used = FALSE,
    observed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package = list(path = expected, version = as.character(packageVersion("sceptre"))),
    ondisc = list(path = find.package("ondisc"), version = as.character(packageVersion("ondisc"))),
    sceptredata = list(path = find.package("sceptredata"), version = as.character(packageVersion("sceptredata"))),
    chunks = list(), outstanding = c("Whole-book rendering", "Full cis/trans/singleton analysis cells", "External rd7 dataset examples"))
payload <- list(schema_version = 1L, chunks = list())
save_reports <- function() {
    jsonlite::write_json(report, "book-report.json", pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
    saveRDS(payload, "book-comparison.rds", version = 3L)
}
chunks <- function(lines) {
    result <- list(); i <- 1L
    while (i <= length(lines)) {
        if (grepl("^```\\{r(?:[, }]|$)", lines[i], perl = TRUE)) {
            start <- i
            end <- i + which(grepl("^```[[:space:]]*$", lines[seq.int(i + 1L, length(lines))]))[1L]
            result[[length(result) + 1L]] <- list(start = start, end = end, header = lines[start], code = paste(lines[seq.int(start + 1L, end - 1L)], collapse = "\n"))
            i <- end
        }
        i <- i + 1L
    }
    result
}
observe <- function(env, chapter) {
    if (chapter == "import-data.qmd" && exists("sceptre_object", env, inherits = FALSE)) {
        object <- get("sceptre_object", env)
        response <- sceptre::get_response_matrix(object)
        guides <- sceptre::get_grna_matrix(object)
        covariates <- sceptre::get_cell_covariates(object)
        return(list(response_dimensions = dim(response), guide_dimensions = dim(guides),
            response_class = class(response), guide_class = class(guides),
            response_sha256 = digest::digest(response, algo = "sha256"),
            guide_sha256 = digest::digest(guides, algo = "sha256"),
            covariate_names = names(covariates), covariates = as.data.frame(covariates)))
    }
    if (chapter == "ondisc.qmd" && exists("gene_odm", env, inherits = FALSE)) {
        object <- get("gene_odm", env)
        answer <- list(dimensions = dim(object), rownames = rownames(object))
        if (exists("expression_vector", env, inherits = FALSE)) answer$expression_vector <- get("expression_vector", env)
        return(answer)
    }
    list()
}
main <- function() {
    inventory <- jsonlite::read_json(file.path(dirname(script), "consumer-inventory.json"))
    pin <- Filter(function(x) identical(x$name, "sceptre-book"), inventory$repositories)[[1L]]$revision
    current <- system2("git", c("-C", shQuote(book), "rev-parse HEAD"), stdout = TRUE)
    stopifnot(identical(current, pin))
    report$book_revision <<- pin
    data_root <- system.file("extdata", package = "sceptredata")
    data_files <- c(list.files(file.path(data_root, "highmoi_example"), recursive = TRUE, full.names = TRUE),
        list.files(file.path(data_root, "parse_example"), recursive = TRUE, full.names = TRUE))
    report$input_files <<- stats::setNames(lapply(data_files, function(path) {
        answer <- list(bytes = file.info(path)$size, sha256 = hash(path))
        if (basename(path) == "barcodes.tsv.gz") {
            connection <- gzfile(path)
            answer$barcode_lines <- length(readLines(connection))
            close(connection)
        }
        answer
    }), substring(data_files, nchar(data_root) + 2L))
    set.seed(271828)
    for (chapter in c("import-data.qmd", "ondisc.qmd")) {
        row <- Filter(function(x) identical(x$path, chapter), inventory$book_contracts)[[1L]]
        path <- file.path(book, chapter)
        stopifnot(identical(hash(path), row$sha256))
        candidates <- chunks(readLines(path, warn = FALSE))
        selected <- Filter(function(x) {
            if (chapter == "import-data.qmd") {
                grepl("library(sceptre)", x$code, fixed = TRUE) ||
                    (grepl("directories <- paste0", x$code, fixed = TRUE) && grepl("import_data_from_cellranger(", x$code, fixed = TRUE)) ||
                    grepl("sceptre_object <- import_data_from_parse(", x$code, fixed = TRUE)
            } else x$start >= 16L && x$end <= 109L
        }, candidates)
        stopifnot(length(selected) > 0L)
        env <- new.env(parent = .GlobalEnv)
        for (cell in selected) {
            label <- paste0(chapter, ":", cell$start)
            warnings <- character()
            problem <- tryCatch({
                withCallingHandlers(eval(parse(text = cell$code), envir = env), warning = function(w) {
                    warnings <<- c(warnings, conditionMessage(w))
                })
                stopifnot(identical(normalizePath(find.package("sceptre")), expected))
                observation <- observe(env, chapter)
                payload$chunks[[label]] <<- observation
                NULL
            }, error = function(e) conditionMessage(e))
            report$chunks[[label]] <<- list(status = if (is.null(problem)) "passed" else "failed", evidence = "executed uncached",
                path = chapter, source_sha256 = row$sha256, start_line = cell$start, end_line = cell$end,
                code_sha256 = digest::digest(cell$code, serialize = FALSE, algo = "sha256"), code = cell$code,
                warnings = warnings, error = problem)
            save_reports()
        }
    }
}
fatal <- tryCatch({ main(); NULL }, error = function(e) conditionMessage(e))
report$fatal_error <- fatal
report$passed <- is.null(fatal) && length(report$chunks) > 0L && all(vapply(report$chunks, function(x) identical(x$status, "passed"), logical(1L)))
save_reports()
cat(jsonlite::toJSON(list(passed = report$passed, chunks = length(report$chunks), fatal_error = fatal), auto_unbox = TRUE, null = "null"), "\n")
quit(status = if (report$passed) 0L else 1L)
