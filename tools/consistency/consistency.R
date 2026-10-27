#!/usr/bin/env Rscript
# Frozen behavior suite: developer command. Not part of the package build.
#
#   Rscript tools/consistency/consistency.R run [--lib DIR | --ref REF]
#       [--fixtures DIR] [--filter REGEX] [--tier strict|portable]
#       [--waive-internal]
#   Rscript tools/consistency/consistency.R regenerate --ref REF [--out DIR]
#   Rscript tools/consistency/consistency.R compare --reference REF [--ref REF]
#   Rscript tools/consistency/consistency.R record --lib DIR --out DIR
#       --commit SHA [--ref REF] [--saved-objects DIR]
#
# `run` checks a candidate against the stored fixtures. By default the
# candidate is the working tree: the tracked files, with their uncommitted
# changes (add new files with `git add` first). `--ref` selects a committed
# ref and `--lib` a library where sceptre is installed. `--tier` overrides the
# comparison tier chosen from the platform; `--waive-internal` skips the
# comparison of internal slots (README.md explains when that is allowed).
# `regenerate` records the fixtures from a reference commit. `compare` records
# a reference into a temporary directory and checks a candidate against it.
# `record` is the step that `regenerate` runs in a fresh R process.
# Everything except the requested fixture directory is written to tempdir().

script_file <- function() {
    file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
    normalizePath(file[[1L]], winslash = "/", mustWork = TRUE)
}
suite_dir <- dirname(script_file())
repo_dir <- normalizePath(
    file.path(suite_dir, "..", ".."),
    winslash = "/",
    mustWork = TRUE
)
default_fixtures <- file.path(suite_dir, "fixtures")

usage <- function() {
    lines <- readLines(script_file(), n = 12L)
    stop(
        "Usage:\n",
        paste(sub("^# ?", "", lines[3:9]), collapse = "\n"),
        call. = FALSE
    )
}

parse_args <- function(args) {
    if (!length(args)) {
        usage()
    }
    out <- list(command = args[[1L]])
    rest <- args[-1L]
    while (length(rest)) {
        key <- rest[[1L]]
        if (!startsWith(key, "--")) {
            usage()
        }
        name <- gsub("-", "_", substring(key, 3L))
        # A flag without a value, such as --waive-internal, is TRUE.
        if (length(rest) < 2L || startsWith(rest[[2L]], "--")) {
            out[[name]] <- TRUE
            rest <- rest[-1L]
        } else {
            out[[name]] <- rest[[2L]]
            rest <- rest[-(1:2)]
        }
    }
    out
}

single_thread_env <- function() {
    c(
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        VECLIB_MAXIMUM_THREADS = "1",
        R_DATATABLE_NUM_THREADS = "1",
        # Errors and warnings are compared as text: keep them untranslated.
        LANGUAGE = "en"
    )
}

run_command <- function(command, args, log = NULL, env = character()) {
    old <- Sys.getenv(names(env), unset = NA)
    if (length(env)) {
        do.call(Sys.setenv, as.list(env))
    }
    on.exit(
        {
            for (name in names(old)) {
                if (is.na(old[[name]])) {
                    Sys.unsetenv(name)
                } else {
                    do.call(
                        Sys.setenv,
                        stats::setNames(list(old[[name]]), name)
                    )
                }
            }
        },
        add = TRUE
    )
    output <- if (is.null(log)) "" else log
    status <- system2(command, args, stdout = output, stderr = output)
    if (!identical(as.integer(status), 0L)) {
        stop(
            basename(command),
            " failed (status ",
            status,
            ")",
            if (!is.null(log)) paste0("; see ", log),
            call. = FALSE
        )
    }
    invisible(status)
}

git <- function(...) {
    out <- system2(
        "git",
        c("-C", shQuote(repo_dir), ...),
        stdout = TRUE,
        stderr = TRUE
    )
    if (!is.null(attr(out, "status"))) {
        stop(paste(out, collapse = "\n"), call. = FALSE)
    }
    out
}

resolve_commit <- function(ref) {
    git("rev-parse", "--verify", shQuote(paste0(ref, "^{commit}")))
}

# Source tree of a committed ref or, when ref is NULL, of the working tree's
# tracked files (with uncommitted changes). Untracked files, such as scratch
# directories, never reach the installed package.
source_tree <- function(ref, work) {
    tree <- file.path(work, "source")
    dir.create(tree, recursive = TRUE)
    if (is.null(ref)) {
        files <- git("ls-files", "--cached")
        files <- files[file.exists(file.path(repo_dir, files))]
        for (directory in unique(dirname(files))) {
            dir.create(
                file.path(tree, directory),
                recursive = TRUE,
                showWarnings = FALSE
            )
        }
        copied <- file.copy(file.path(repo_dir, files), file.path(tree, files))
        if (!all(copied)) stop("Could not copy the working tree.")
    } else {
        archive <- file.path(work, "source.tar")
        git(
            "-c",
            "core.autocrlf=false",
            "archive",
            "--format=tar",
            paste0("--output=", shQuote(archive)),
            resolve_commit(ref)
        )
        utils::untar(archive, exdir = tree)
    }
    tree
}

install_sceptre <- function(ref, work) {
    tree <- source_tree(ref, work)
    library_dir <- file.path(work, "library")
    dir.create(library_dir)
    log <- file.path(work, "install.log")
    message(
        "Installing sceptre from ",
        if (is.null(ref)) "the working tree" else ref,
        " (log: ",
        log,
        ")"
    )
    run_command(
        file.path(R.home("bin"), "R"),
        c(
            "CMD",
            "INSTALL",
            "--no-docs",
            paste0("--library=", shQuote(library_dir)),
            shQuote(tree)
        ),
        log = log
    )
    library_dir
}

rscript <- function(args, env = character()) {
    run_command(
        file.path(R.home("bin"), "Rscript"),
        c("--vanilla", args),
        env = c(single_thread_env(), env)
    )
}

literal <- function(x) paste(deparse(x), collapse = "")

suite_files <- c("harness.R", "fixture-io.R", "export.R")

# --------------------------------------------------------------------------

command_record <- function(opts) {
    for (required in c("lib", "out", "commit")) {
        if (is.null(opts[[required]])) usage()
    }
    for (package in c("testthat", "withr", "jsonlite")) {
        if (!requireNamespace(package, quietly = TRUE)) {
            stop("Developer dependency missing: ", package)
        }
    }
    library_dir <- normalizePath(opts$lib, winslash = "/", mustWork = TRUE)
    .libPaths(c(library_dir, .libPaths()))
    suppressPackageStartupMessages(library(sceptre))
    if (
        !identical(
            normalizePath(find.package("sceptre"), winslash = "/"),
            normalizePath(file.path(library_dir, "sceptre"), winslash = "/")
        )
    ) {
        stop("sceptre was not loaded from ", library_dir)
    }
    out <- opts$out
    if (
        dir.exists(out) &&
            length(list.files(out, all.files = TRUE, no.. = TRUE))
    ) {
        stop("The output directory must be new or empty: ", out)
    }
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    out <- normalizePath(out, winslash = "/", mustWork = TRUE)
    for (file in suite_files) {
        source(file.path(suite_dir, file), local = globalenv())
    }

    saved_source <- if (is.null(opts$saved_objects)) {
        file.path(default_fixtures, "saved-objects")
    } else {
        opts$saved_objects
    }
    dir.create(file.path(out, "saved-objects"))
    saved <- file.path(saved_source, consistency_saved_object_file)
    if (
        !identical(consistency_sha256(saved), consistency_saved_object_sha256)
    ) {
        stop("The saved-object input does not have its recorded SHA256.")
    }
    file.copy(saved, file.path(out, "saved-objects"))

    consistency_write_inputs(file.path(out, "inputs"))
    cases <- consistency_cases()
    entries <- list()
    for (name in names(cases)) {
        started <- Sys.time()
        observation <- consistency_run_case(name, out)
        seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
        rds <- file.path("cases", paste0(name, ".rds"))
        dir.create(file.path(out, "cases"), showWarnings = FALSE)
        saveRDS(observation, file.path(out, rds), version = 3L, compress = "xz")
        written <- consistency_export_case(name, observation, out)
        entries[[name]] <- list(
            group = cases[[name]]$group,
            description = cases[[name]]$description,
            known_bug = consistency_known_bugs[[name]],
            error = observation$error,
            warnings = observation$warnings,
            rds = rds,
            tables = written
        )
        message(sprintf("%-38s %6.2fs", name, seconds))
    }

    imports <- c(
        "cowplot",
        "crayon",
        "data.table",
        "dplyr",
        "ggplot2",
        "Matrix",
        "Rcpp",
        "BH",
        "scales",
        "purrr",
        "parallelly",
        "withr",
        "testthat",
        "jsonlite"
    )
    versions <- lapply(c("sceptre", imports), function(package) {
        if (requireNamespace(package, quietly = TRUE)) {
            as.character(utils::packageVersion(package))
        } else {
            NULL
        }
    })
    names(versions) <- c("sceptre", imports)
    files <- list.files(out, recursive = TRUE)
    files <- sort(setdiff(files, "manifest.json"), method = "radix")
    manifest <- list(
        schema_version = consistency_schema_version,
        suite = "sceptre frozen behavior suite (tools/consistency)",
        reference = list(
            ref = if (is.null(opts$ref)) opts$commit else opts$ref,
            commit = opts$commit,
            package_version = versions$sceptre
        ),
        generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "UTC"),
        environment = list(
            R = R.version.string,
            platform = R.version$platform,
            os = utils::sessionInfo()$running,
            RNGkind = RNGkind(),
            lapack = La_version(),
            blas = extSoftVersion()[["BLAS"]],
            packages = versions
        ),
        harness = stats::setNames(
            as.list(consistency_sha256(file.path(suite_dir, suite_files))),
            suite_files
        ),
        seeds = list(
            inputs = consistency_input_seed,
            cases = consistency_case_seed
        ),
        comparison = list(
            doubles = paste(
                "elementwise relative: |actual - expected| <=",
                "tolerance * |expected|; where expected is exactly 0,",
                "|actual| <= 1e-12"
            ),
            relative_tolerance = as.list(consistency_tolerances),
            tier_policy = paste(
                "strict when R.version.string and platform equal",
                "environment$R and environment$platform; portable otherwise"
            ),
            exact = c(
                "integer, logical, character and factor values",
                "names, dimensions, classes, row names, NA/NaN/Inf positions",
                "errors, warnings and the RNG state after each case"
            )
        ),
        conventions = list(
            cell_and_feature_indices = "1-based, in the imported (input) order",
            tsv = "UTF-8, tab-separated, header row, NA for missing",
            doubles = "%.17g (round-trips exactly)"
        ),
        inputs = consistency_input_digests(),
        known_bug_cases = names(consistency_known_bugs),
        cases = entries,
        files = stats::setNames(
            as.list(consistency_sha256(file.path(out, files))),
            files
        )
    )
    jsonlite::write_json(
        manifest,
        file.path(out, "manifest.json"),
        pretty = TRUE,
        auto_unbox = TRUE,
        null = "null",
        digits = NA
    )
    message("Recorded ", length(entries), " cases in ", out)
}

command_regenerate <- function(opts) {
    if (is.null(opts$ref)) {
        usage()
    }
    commit <- resolve_commit(opts$ref)
    work <- tempfile("sceptre-consistency-reference-")
    dir.create(work)
    library_dir <- install_sceptre(commit, work)
    staging <- file.path(work, "fixtures")
    rscript(c(
        shQuote(script_file()),
        "record",
        "--lib",
        shQuote(library_dir),
        "--out",
        shQuote(staging),
        "--commit",
        commit,
        "--ref",
        shQuote(opts$ref),
        "--saved-objects",
        shQuote(file.path(default_fixtures, "saved-objects"))
    ))
    out <- if (is.null(opts$out)) default_fixtures else opts$out
    if (
        normalizePath(out, mustWork = FALSE) == normalizePath(default_fixtures)
    ) {
        # Replace the generated fixtures; inputs under saved-objects/ stay.
        unlink(
            file.path(out, c("manifest.json", "inputs", "cases", "tables")),
            recursive = TRUE
        )
    } else if (
        dir.exists(out) &&
            length(list.files(out, all.files = TRUE, no.. = TRUE))
    ) {
        stop("The output directory must be new or empty: ", out)
    }
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    for (item in list.files(staging)) {
        if (item == "saved-objects" && dir.exists(file.path(out, item))) {
            next
        }
        file.copy(file.path(staging, item), out, recursive = TRUE)
    }
    message("Fixtures for ", opts$ref, " (", commit, ") written to ", out)
    invisible(out)
}

command_run <- function(opts) {
    fixtures <- if (is.null(opts$fixtures)) default_fixtures else opts$fixtures
    fixtures <- normalizePath(fixtures, winslash = "/", mustWork = TRUE)
    library_dir <- if (!is.null(opts$lib)) {
        normalizePath(opts$lib, winslash = "/", mustWork = TRUE)
    } else {
        work <- tempfile("sceptre-consistency-candidate-")
        dir.create(work)
        install_sceptre(opts$ref, work)
    }
    filter <- if (is.null(opts$filter)) NULL else opts$filter
    env <- c(SCEPTRE_CONSISTENCY_FIXTURES = fixtures)
    if (!is.null(opts$tier)) {
        if (!opts$tier %in% c("strict", "portable")) {
            usage()
        }
        env[["SCEPTRE_CONSISTENCY_TIER"]] <- opts$tier
    }
    if (isTRUE(opts$waive_internal)) {
        env[["SCEPTRE_CONSISTENCY_WAIVE_INTERNAL"]] <- "true"
    }
    code <- c(
        paste0(".libPaths(c(", literal(library_dir), ", .libPaths()))"),
        "suppressPackageStartupMessages(library(sceptre))",
        "testthat::set_max_fails(Inf)",
        paste0(
            "stopifnot(identical(normalizePath(find.package('sceptre')), ",
            "normalizePath(file.path(",
            literal(library_dir),
            ", 'sceptre'))))"
        ),
        paste0(
            "results <- testthat::test_dir(",
            literal(file.path(suite_dir, "testthat")),
            ", package = 'sceptre', load_package = 'installed', ",
            "filter = ",
            literal(filter),
            ", stop_on_failure = FALSE)"
        ),
        "counts <- as.data.frame(results)",
        paste0(
            "if (sum(counts$failed) || any(counts$error)) ",
            "quit(status = 1L, save = 'no')"
        )
    )
    worker <- tempfile("sceptre-consistency-run-", fileext = ".R")
    writeLines(code, worker)
    rscript(shQuote(worker), env = env)
    message("Frozen behavior suite passed against ", fixtures)
}

command_compare <- function(opts) {
    if (is.null(opts$reference)) {
        usage()
    }
    out <- tempfile("sceptre-consistency-fixtures-")
    command_regenerate(list(ref = opts$reference, out = out))
    opts$fixtures <- out
    command_run(opts)
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
for (name in names(single_thread_env())) {
    if (!nzchar(Sys.getenv(name))) {
        do.call(Sys.setenv, as.list(single_thread_env()[name]))
    }
}
switch(
    opts$command,
    record = command_record(opts),
    regenerate = command_regenerate(opts),
    run = command_run(opts),
    compare = command_compare(opts),
    usage()
)
