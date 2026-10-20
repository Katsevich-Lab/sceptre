#!/usr/bin/env Rscript
# Run from the repository root with an explicitly isolated installed library.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: Rscript tools/validation/run-scientific.R OUTPUT REPLICATES [N_CELLS] [CASE_IDS] [SETUP_R_OR_DASH] [FIRST_REPLICATE]")
out <- normalizePath(args[1L], mustWork = FALSE)
reps <- as.integer(args[2L])
n <- if (length(args) >= 3L) as.integer(args[3L]) else 600L
stopifnot(length(reps) == 1L, !is.na(reps), reps > 0L, reps <= 100000L,
          length(n) == 1L, !is.na(n), n >= 100L)
if (dir.exists(out)) stop("Output directory already exists; choose a new run directory.")
source("tools/validation/scientific-harness.R")
cases <- validation_cases()
if (length(args) >= 4L && args[4L] != "all") {
    requested <- strsplit(args[4L], ",", fixed = TRUE)[[1L]]
    stopifnot(all(requested %in% cases$id))
    cases <- cases[match(requested, cases$id), , drop = FALSE]
}
first <- if (length(args) >= 6L) as.integer(args[6L]) else 1L
stopifnot(length(first) == 1L, !is.na(first), first >= 1L, first + reps - 1L <= 100000L)
replicate_ids <- seq.int(first, length.out = reps)
setup <- if (length(args) >= 5L && args[5L] != "-") normalizePath(args[5L], mustWork = TRUE) else NULL
if (!is.null(setup)) source(setup)
library(sceptre)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
dir.create(out, recursive = TRUE)
revision <- function() paste(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), collapse = "")
installed <- find.package("sceptre")
installed_files <- list.files(installed, recursive = TRUE, full.names = TRUE)
source_files <- c("DESCRIPTION", "NAMESPACE",
    list.files("R", pattern = "[.]R$", full.names = TRUE),
    list.files("src", pattern = "[.](cpp|h)$", full.names = TRUE))
manifest <- list(started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    source_revision = revision(), source_status = system2("git", c("status", "--porcelain"), stdout = TRUE),
    source_root = getwd(), installed_package = installed,
    source_sha256 = stats::setNames(vapply(source_files, digest::digest, character(1),
        algo = "sha256", file = TRUE), source_files),
    installed_sha256 = stats::setNames(vapply(installed_files, digest::digest, character(1),
        algo = "sha256", file = TRUE), substring(installed_files, nchar(installed) + 2L)),
    harness_sha256 = vapply(c("tools/validation/scientific-harness.R", "tools/validation/run-scientific.R"),
                            digest::digest, character(1), algo = "sha256", file = TRUE),
    setup = setup, setup_sha256 = if (!is.null(setup)) digest::digest(setup, file = TRUE, algo = "sha256") else NULL,
    cases = cases, replicates = reps, first_replicate = first,
    last_replicate = first + reps - 1L, n_cells = n, rng = RNGkind(),
    seed_version = "v2-disjoint-case-ranges", seed_rule = validation_seed_rule(),
    replication_unit = "One independently regenerated cell dataset per case and replicate; endpoints never pooled for binomial precision",
    interval = "Pointwise exact Clopper-Pearson 95% conditional on valid p-values; failures/QC and per-attempt rejection separately reported",
    session = utils::sessionInfo())
if (!is.null(setup)) file.copy(setup, file.path(out, "candidate-setup.R"))
saveRDS(manifest, file.path(out, "manifest.rds"))
capture.output(str(manifest), file = file.path(out, "manifest.txt"))
all_rows <- list()
for (i in seq_len(nrow(cases))) {
    case <- cases[i, , drop = FALSE]
    for (r in replicate_ids) {
        seed <- validation_seed(case$id, r)
        result <- validation_replicate(case, r, seed, n)
        all_rows[[length(all_rows) + 1L]] <- result
        utils::write.table(result, file.path(out, "replicates.csv"), sep = ",", row.names = FALSE,
                           col.names = length(all_rows) == 1L, append = length(all_rows) != 1L, na = "NA")
        if (r %% 25L == 0L || r == tail(replicate_ids, 1L)) cat(case$id, r - first + 1L, "/", reps, "replicate", r, "\n")
    }
    utils::write.csv(validation_summary(do.call(rbind, all_rows)), file.path(out, "summary.csv"), row.names = FALSE)
}
manifest$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
manifest$replicates_sha256 <- digest::digest(file.path(out, "replicates.csv"), file = TRUE, algo = "sha256")
saveRDS(manifest, file.path(out, "manifest.rds"))
capture.output(str(manifest), file = file.path(out, "manifest.txt"))
cat("Completed", nrow(cases) * reps, "independent case-datasets; summary:", file.path(out, "summary.csv"), "\n")
