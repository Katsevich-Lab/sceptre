# Metadata qualification only: this does not execute R 4.1 or its dependencies.
# Usage: Rscript .../audit-r-floor-snapshot.R [downloaded PACKAGES.gz]
repo <- "https://packagemanager.posit.co/cran/2024-03-01"
args <- commandArgs(trailingOnly = TRUE)
index <- if (length(args)) args[[1L]] else tempfile(fileext = ".gz")
if (!length(args)) download.file(paste0(repo, "/src/contrib/PACKAGES.gz"), index)
db <- read.dcf(gzfile(index))
compatible <- function(spec) {
    if (is.na(spec)) return(TRUE)
    required <- tools:::.split_dependencies(spec)$R
    if (is.null(required) || is.null(required$version)) return(TRUE)
    do.call(required$op, list(package_version("4.1.0"), package_version(required$version)))
}
# Snapshot indices include multiple R-specific Recommended rows. Match the
# R-version filtering performed by available.packages(), rather than choosing
# a new-R recommended row just because it appears first in the index.
db <- db[vapply(db[, "Depends"], compatible, logical(1)), , drop = FALSE]
db <- db[!duplicated(db[, "Package"]), , drop = FALSE]
d <- read.dcf("DESCRIPTION")
roots <- unique(unlist(lapply(d[1L, c("Imports", "LinkingTo", "Suggests")],
    function(x) names(tools:::.split_dependencies(x)))))
roots <- setdiff(roots, "ondisc")
closure <- unique(c(roots, unlist(tools::package_dependencies(roots, db = db,
    which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE))))
closure <- setdiff(closure, c("base", "compiler", "datasets", "graphics", "grDevices",
    "grid", "methods", "parallel", "splines", "stats", "stats4", "tcltk", "tools", "utils"))
rows <- db[match(sort(closure), db[, "Package"]),
    c("Package", "Version", "Depends", "Imports", "LinkingTo"), drop = FALSE]
stopifnot(!anyNA(rows[, "Package"]))
write.table(rows, "tools/hardening/dependencies/evidence/r41-snapshot-dependency-closure.tsv",
    sep = "\t", quote = TRUE, row.names = FALSE)
cat(nrow(rows), "R-compatible repository entries in recursive closure; testthat separately pinned 3.1.7.\n")
