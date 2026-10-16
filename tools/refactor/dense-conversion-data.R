#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply a fresh dense-cost fixture directory.")
output <- args[[1L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(output))) stop("Dense-cost fixtures require a fresh directory.")
hash <- function(path) digest::digest(file = path, algo = "sha256")
manifest <- list(recipe = "DENSE_CONVERSION_COST.md", R = R.version.string, fixtures = list())
for (shape in list(c(1000L, 1500L), c(2000L, 3000L))) {
    elements <- prod(shape)
    input <- matrix(0, shape[[1L]], shape[[2L]],
        dimnames = list(sprintf("r%05d", seq_len(shape[[1L]])),
                        sprintf("c%05d", seq_len(shape[[2L]]))))
    positions <- seq.int(1L, elements, by = 10L)
    input[positions] <- ((seq_along(positions) %% 97L) + 1L) / 7
    name <- paste0(elements, "-dense")
    path <- file.path(output, paste0(name, ".rds"))
    saveRDS(input, path, compress = FALSE, version = 3)
    manifest$fixtures[[name]] <- list(file = basename(path), sha256 = hash(path),
        class = "matrix", dimensions = shape, elements = elements,
        input_payload_bytes = elements * 8, nnz = length(positions), density = 0.1,
        values_range = range(input))
    rm(input, positions)
    gc()
}
jsonlite::write_json(manifest, file.path(output, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
