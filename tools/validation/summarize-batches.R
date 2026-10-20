#!/usr/bin/env Rscript
# Aggregate independent-index batches only after validating terminal identities.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: summarize-batches.R NEW_OUTPUT_DIR BATCH_DIR [BATCH_DIR ...]")
source("tools/validation/scientific-harness.R")
out <- normalizePath(args[1L], mustWork = FALSE)
if (dir.exists(out)) stop("Output directory already exists")
batches <- normalizePath(args[-1L], mustWork = TRUE)
records <- lapply(batches, function(path) {
    m <- readRDS(file.path(path, "manifest.rds"))
    if (is.null(m$finished_utc)) stop("Incomplete batch: ", path)
    hash <- digest::digest(file.path(path, "replicates.csv"), algo = "sha256", file = TRUE)
    if (!identical(hash, m$replicates_sha256)) stop("Changed batch CSV: ", path)
    rows <- utils::read.csv(file.path(path, "replicates.csv"), stringsAsFactors = FALSE)
    if (nrow(rows) != m$replicates * nrow(m$cases) * 3L) stop("Incomplete batch rows: ", path)
    list(manifest = m, rows = rows)
})
# Scientific source, installed binary, DGP and hook identity must match exactly.
identity <- function(m) m[c("source_sha256", "installed_sha256", "harness_sha256",
    "setup_sha256", "cases", "n_cells", "rng", "seed_rule")]
reference <- identity(records[[1L]]$manifest)
if (!all(vapply(records, function(x) identical(identity(x$manifest), reference), logical(1))))
    stop("Batch scientific source/binary/settings differ")
rows <- do.call(rbind, lapply(records, `[[`, "rows"))
key <- paste(rows$case, rows$replicate, rows$response_id, sep = "/")
if (anyDuplicated(key)) stop("Duplicated case/replicate/endpoint: do not pool repeated pilot seeds")
# Every global replicate must contain the same three requested endpoints.
if (any(table(paste(rows$case, rows$replicate, sep = "/")) != 3L)) stop("Incomplete replicate")
dir.create(out, recursive = TRUE)
utils::write.csv(rows, file.path(out, "replicates.csv"), row.names = FALSE)
utils::write.csv(validation_summary(rows), file.path(out, "summary.csv"), row.names = FALSE)
saveRDS(list(batches = batches, input_manifests = lapply(records, `[[`, "manifest"),
    input_csv_sha256 = vapply(records, function(x) x$manifest$replicates_sha256, character(1)),
    replicated_datasets_per_case = table(unique(rows[c("case", "replicate")])$case),
    combined_csv_sha256 = digest::digest(file.path(out, "replicates.csv"), file = TRUE, algo = "sha256"),
    finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(out, "combined-manifest.rds"))
cat("Combined", length(unique(paste(rows$case, rows$replicate))), "independent case-datasets\n")
