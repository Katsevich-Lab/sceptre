#!/usr/bin/env Rscript
# Explicit temporary handoff command; never called by package tests or CI.
# Render a committed source snapshot into two self-contained HTML pages and a ZIP.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) stop("Usage: render.R /fresh/output [git-ref; default HEAD]")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite first.")
quarto <- Sys.which("quarto")
if (!nzchar(quarto)) stop("Install Quarto first; see the handoff README.")
if (!nzchar(Sys.which("zip"))) stop("Install the zip command first.")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
repo <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
out <- path.expand(args[[1L]])
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a fresh empty output directory.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, winslash = "/", mustWork = TRUE)
if (startsWith(paste0(out, "/"), paste0(repo, "/"))) stop("Write outside the checkout.")
ref <- if (length(args) == 2L) args[[2L]] else "HEAD"
pin <- system2("git", c("-C", shQuote(repo), "rev-parse", "--verify", shQuote(paste0(ref, "^{commit}"))), stdout = TRUE)
if (!is.null(attr(pin, "status"))) stop("Git ref does not resolve.")
archive <- file.path(out, "source.tar")
stopifnot(system2("git", c("-C", shQuote(repo), "archive", paste0("--output=", shQuote(archive)), pin)) == 0L)
source <- file.path(out, "source"); dir.create(source)
utils::untar(archive, exdir = source)
setwd(file.path(source, "tools/refactor"))
for (name in c("walkthrough", "bug-fixes")) {
    log <- file.path(out, paste0(name, "-render.log"))
    status <- system2(quarto, c("render", paste0(name, ".qmd")), stdout = log, stderr = log)
    if (status != 0L) stop("Render failed; read ", log)
    path <- paste0(name, ".html")
    html <- readLines(path, warn = FALSE)
    links <- unique(unlist(regmatches(html, gregexpr('href="[^"#]+"', html))))
    for (link in links) {
        target <- substring(link, 7L, nchar(link) - 1L)
        if (grepl("^[a-z]+:|^/", target) || grepl("\\.html(#.*)?$", target)) next
        candidate <- normalizePath(file.path(getwd(), sub("#.*$", "", target)), winslash = "/", mustWork = FALSE)
        if (!startsWith(candidate, paste0(source, "/")) || !file.exists(candidate)) next
        relative <- substring(candidate, nchar(source) + 2L)
        anchor <- if (grepl("#", target)) sub("^[^#]*", "", target) else ""
        replacement <- paste0('href="https://github.com/Katsevich-Lab/sceptre/blob/', pin, "/", relative, anchor, '"')
        html <- gsub(link, replacement, html, fixed = TRUE)
    }
    writeLines(html, file.path(out, path))
}
writeLines(c("SCEPTRE refactor reading guide", paste("Source commit:", pin), "",
    "Open walkthrough.html. The bug appendix is bug-fixes.html.",
    "Both pages are self-contained; keep them together for their cross-links.",
    "Source/test links point to the exact committed source on GitHub.",
    "These pages do not run SCEPTRE or constitute new numerical qualification."), file.path(out, "READ-ME.txt"))
jsonlite::write_json(list(source_commit = pin, quarto = system2(quarto, "--version", stdout = TRUE),
    html_md5 = as.list(tools::md5sum(file.path(out, c("walkthrough.html", "bug-fixes.html"))))),
    file.path(out, "render.json"), pretty = TRUE, auto_unbox = TRUE)
setwd(out)
status <- utils::zip("sceptre-refactor-tour.zip", c("walkthrough.html", "bug-fixes.html", "READ-ME.txt", "render.json"))
if (status != 0L) stop("Could not create the handout ZIP.")
cat("Open ", file.path(out, "walkthrough.html"), " or share sceptre-refactor-tour.zip\n", sep = "")
