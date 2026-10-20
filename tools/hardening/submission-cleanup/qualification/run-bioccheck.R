# Run in a fresh process with the requested checker library first in R_LIBS.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
archive <- normalizePath(args[[1L]])
output <- args[[2L]]
checker <- args[[3L]]
bioc <- args[[4L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output)
stopifnot(as.character(packageVersion("BiocCheck")) == checker)
source_archive <- archive
archive <- file.path(output, basename(source_archive))
if (!file.exists(archive)) {
    stopifnot(file.copy(source_archive, archive))
}
stopifnot(identical(
    digest::digest(file = source_archive, algo = "sha256"),
    digest::digest(file = archive, algo = "sha256")))
Sys.setenv(R_BIOC_VERSION = bioc)
cache <- file.path(output, "cache")
dir.create(cache, showWarnings = FALSE)
BiocFileCache::setBFCOption("CACHE", cache)
writeLines(c(R.version.string, capture.output(sessionInfo())),
    file.path(output, "session.txt"))
jsonlite::write_json(list(
    source_archive = source_archive,
    archive = archive,
    archive_sha256 = digest::digest(file = archive, algo = "sha256"),
    checker = checker,
    cache = cache,
    Bioconductor = bioc,
    flags = c("no-check-R-ver", "no-check-bioc-help"),
    platform = R.version$platform
), file.path(output, "invocation.json"), pretty = TRUE, auto_unbox = TRUE)
setwd(output)
result <- BiocCheck::BiocCheck(archive,
    debug = TRUE, callr = FALSE,
    `no-check-R-ver` = TRUE, `no-check-bioc-help` = TRUE,
    `quit-with-status` = FALSE)
saveRDS(result, "result.rds")
result$toJSON("result.json")
