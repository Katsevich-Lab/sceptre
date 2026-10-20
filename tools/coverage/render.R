# Summaries count unique physical executable lines, separately for R/native.
# No exclusions beyond covr's standard generated Rcpp wrappers.
render_coverage <- function(coverage, root, out) {
    lines <- covr::tally_coverage(coverage, by = "line")
    lines <- stats::aggregate(value ~ filename + line, lines, min)
    lines <- lines[order(lines$filename, lines$line), ]
    lines$language <- ifelse(grepl("^R/", lines$filename), "R", "native")
    utils::write.csv(lines, file.path(out, "lines.csv"), row.names = FALSE)
    summarise <- function(x) data.frame(executable = nrow(x), covered = sum(x$value > 0),
        percent = if (nrow(x)) 100 * mean(x$value > 0) else NA_real_)
    files <- do.call(rbind, lapply(split(lines, lines$filename), function(x)
        cbind(file = x$filename[[1L]], language = x$language[[1L]], summarise(x))))
    rownames(files) <- NULL
    utils::write.csv(files, file.path(out, "files.csv"), row.names = FALSE)
    totals <- lapply(c("R", "native"), function(lang) summarise(lines[lines$language == lang, ]))
    names(totals) <- c("R", "native")
    inventory <- list.files(root, pattern = "\\.(R|cpp|h|hpp)$", recursive = TRUE)
    inventory <- inventory[grepl("^(R|src)/", inventory)]
    generated <- c("R/RcppExports.R", "src/RcppExports.cpp")
    missing <- setdiff(inventory, c(files$file, generated))
    sha <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)
    status <- system2("git", c("-C", shQuote(root), "status", "--porcelain"), stdout = TRUE)
    inputs <- c(list.files(file.path(root, "R"), full.names = TRUE, recursive = TRUE),
        list.files(file.path(root, "src"), pattern = "\\.(cpp|h|hpp)$", full.names = TRUE, recursive = TRUE),
        list.files(file.path(root, "tests"), full.names = TRUE, recursive = TRUE))
    hashes <- as.list(unname(tools::md5sum(inputs)))
    names(hashes) <- substring(inputs, nchar(root) + 2L)
    state_file <- file.path(out, "source-state.rds")
    if (file.exists(state_file)) {
        state <- readRDS(state_file)
        sha <- state$commit
        status <- state$dirty
    }
    metadata <- list(commit = sha, dirty = status, measured_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
        test_scope = "Installed package testthat suite; tests only, not examples/vignettes or standalone native bridges",
        metric = "Unique executable physical lines hit at least once; not branch coverage or statistical validity",
        exclusions = generated, totals = totals, unmeasured_source_files = missing,
        R = R.version.string, platform = R.version$platform, covr = as.character(utils::packageVersion("covr")),
        compiler_flags = as.list(getOption("covr.flags")), gcov = getOption("covr.gcov"), gcov_args = getOption("covr.gcov_args"),
        ondisc = if (requireNamespace("ondisc", quietly = TRUE)) as.character(utils::packageVersion("ondisc")) else NULL,
        ondisc_remote_sha = if (requireNamespace("ondisc", quietly = TRUE)) utils::packageDescription("ondisc")$RemoteSha else NULL,
        ondisc_description_md5 = if (requireNamespace("ondisc", quietly = TRUE)) unname(tools::md5sum(system.file("DESCRIPTION", package = "ondisc"))) else NULL,
        source_test_md5 = hashes)
    jsonlite::write_json(metadata, file.path(out, "summary.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
    capture.output(sessionInfo(), file = file.path(out, "session-info.txt"))
    covr::to_cobertura(coverage, filename = file.path(out, "cobertura.xml"))
    escape <- function(x) {
        x <- gsub("&", "&amp;", x, fixed = TRUE)
        x <- gsub("<", "&lt;", x, fixed = TRUE)
        gsub(">", "&gt;", x, fixed = TRUE)
    }
    body <- c("<!doctype html><html><head><meta charset='utf-8'><title>SCEPTRE coverage</title>",
        "<style>body{font:16px system-ui;max-width:1200px;margin:2em auto}pre{font:13px monospace;overflow:auto}.miss{background:#ffe0df}.hit{background:#e8f2ff}summary{cursor:pointer;padding:.5em}a{color:#164b94}</style></head><body>",
        "<h1>SCEPTRE installed-package test coverage</h1>", paste0("<p>Source commit: <code>", sha, "</code>. See summary.json for dirty state and exact input hashes.</p>"),
        "<p>Executable lines, not branches or a scientific-validity score. Generated Rcpp wrappers excluded by covr. Unmeasured files are listed below; absence is not zero coverage.</p>")
    for (lang in names(totals)) {
        z <- totals[[lang]]
        label <- if (is.na(z$percent)) "unmeasured" else sprintf("%.2f%% (%s/%s)", z$percent, z$covered, z$executable)
        body <- c(body, paste0("<p><b>", lang, ": ", label, "</b></p>"))
        badge <- sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="240" height="24" role="img" aria-label="%s line coverage %s"><rect width="240" height="24" fill="#435467"/><text x="120" y="17" fill="white" text-anchor="middle" font-family="sans-serif" font-size="12">%s lines: %s</text></svg>', lang, escape(label), lang, escape(label))
        writeLines(badge, file.path(out, paste0(tolower(lang), "-coverage.svg")))
    }
    body <- c(body, paste0("<p>Unmeasured source files: ", escape(paste(missing, collapse = ", ")), "</p>"))
    for (f in files$file[order(files$percent)]) {
        z <- files[files$file == f, ]
        source <- readLines(file.path(root, f), warn = FALSE)
        hits <- lines[lines$filename == f, ]
        rows <- vapply(seq_along(source), function(i) {
            value <- hits$value[match(i, hits$line)]
            cls <- if (is.na(value)) "" else if (value > 0) "hit" else "miss"
            sprintf('<span class="%s">%5d %s</span>', cls, i, escape(source[[i]]))
        }, character(1))
        body <- c(body, sprintf("<details><summary>%s: %.2f%% (%s/%s)</summary><pre>", escape(f), z$percent, z$covered, z$executable), rows, "</pre></details>")
    }
    writeLines(c(body, "</body></html>"), file.path(out, "index.html"))
    print(files[order(files$percent), ], row.names = FALSE)
    print(totals)
    invisible(metadata)
}
