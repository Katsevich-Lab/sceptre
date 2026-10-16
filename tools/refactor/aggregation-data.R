#!/usr/bin/env Rscript
# Deterministic sparse inputs only; no expected output or full dense assay.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: aggregation-data.R output-directory")
directory <- args[1L]
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(Matrix))
set.seed(20261016L)
metadata <- function(ids, repeats = FALSE) {
    index <- rev(seq_along(ids))
    result <- data.frame(grna_id = ids[index],
        grna_target = paste0("target", (index - 1L) %/% 4L + 1L),
        vector_id = paste0("v", (index - 1L) %/% 4L + 1L))
    if (repeats) result <- rbind(result, result[seq.int(1L, nrow(result), by = 50L), ])
    result
}
synthetic <- function(guides, cells, triplet = FALSE, logical = FALSE) {
    per_guide <- as.integer(cells / 200L)
    i <- rep(seq_len(guides), each = per_guide)
    j <- as.integer(((rep(seq_len(per_guide) - 1L, guides) * 197 + i * 311) %% cells) + 1L)
    x <- as.numeric(1L + (seq_along(i) * 13L) %% 17L)
    if (triplet) {
        x <- x / 8
        x[seq.int(1L, length(x), by = 31L)] <- 0
        duplicate <- seq.int(1L, length(i), by = 100L)
        i <- c(i, i[duplicate]); j <- c(j, j[duplicate]); x <- c(x, rep(.25, length(duplicate)))
        order <- rev(seq_along(i)); i <- i[order]; j <- j[order]; x <- x[order]
    }
    matrix <- sparseMatrix(i = i, j = j, x = x, dims = c(guides, cells),
        dimnames = list(paste0("g", seq_len(guides)), paste0("c", seq_len(cells))),
        repr = if (triplet) "T" else "C")
    if (logical) matrix <- methods::new("lgCMatrix", i = matrix@i, p = matrix@p,
        x = seq_along(matrix@x) %% 8L != 0L, Dim = matrix@Dim, Dimnames = matrix@Dimnames)
    stopifnot(methods::validObject(matrix))
    list(m = matrix, d = metadata(rownames(matrix), repeats = triplet),
        recipe = list(kind = "deterministic synthetic aggregation workload", guides = guides,
            cells = cells, stored_entries = length(matrix@x), stored_density = length(matrix@x) / (guides * cells),
            nonzero_values = sum(matrix@x != 0), four_guides_per_vector = TRUE,
            metadata_repetitions = if (triplet) "one additional selection every 50 metadata rows" else "none",
            triplet_duplicates = triplet, logical_stored_FALSE = logical,
            seed = 20261016L, source_cell_order = "c1...cN", vector_order = "first-seen reversed metadata",
            dense_input_counterfactual_bytes = as.double(guides) * cells * (if (logical) 4 else 8),
            dense_output_bytes = (guides / 4) * cells * 8))
}
write_case <- function(name, input) {
    path <- file.path(directory, paste0(name, ".rds"))
    saveRDS(input, path, compress = FALSE)
    list(name = name, file = basename(path), sha256 = digest::digest(file = path, algo = "sha256"), recipe = input$recipe)
}
cases <- list(write_case("large_numeric_C", synthetic(3200L, 24000L)),
              write_case("medium_numeric_T_repeats", synthetic(800L, 6000L, triplet = TRUE)),
              write_case("medium_logical_C", synthetic(800L, 6000L, logical = TRUE)))
real_status <- "sceptredata unavailable; real-data setting not executed"
if (requireNamespace("sceptredata", quietly = TRUE)) {
    scope <- new.env(parent = emptyenv())
    utils::data("lowmoi_example_data", package = "sceptredata", envir = scope)
    real <- scope$lowmoi_example_data$grna_matrix
    d <- metadata(rownames(real))
    input <- list(m = real, d = d, recipe = list(kind = "actual Papalexi example guide counts; synthetic groups of four",
        source = "sceptredata::lowmoi_example_data$grna_matrix", sceptredata_version = as.character(packageVersion("sceptredata")),
        guides = nrow(real), cells = ncol(real),
        stored_entries = if (methods::is(real, "sparseMatrix")) length(real@x) else sum(real != 0),
        stored_density = if (methods::is(real, "sparseMatrix")) length(real@x) / prod(dim(real)) else mean(real != 0), vector_count = length(unique(d$vector_id)),
        real_vector_metadata = FALSE, dense_input_counterfactual_bytes = prod(as.double(dim(real))) * 8,
        dense_output_bytes = length(unique(d$vector_id)) * ncol(real) * 8))
    cases[[length(cases) + 1L]] <- write_case("actual_papalexi_guides", input)
    real_status <- "Actual local guide counts included; vector grouping is synthetic, not biological evidence."
}
jsonlite::write_json(list(schema_version = 1L, cases = cases, real_status = real_status,
    generator = "aggregation-data.R", generator_sha256 = digest::digest(file = sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)), algo = "sha256"),
    R = R.version.string, Matrix = as.character(packageVersion("Matrix"))), file.path(directory, "manifest.json"),
    pretty = TRUE, auto_unbox = TRUE, null = "null")
