#!/usr/bin/env Rscript
# Installed-namespace replay; outputs are temporary experiment results, not goldens.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: aggregation-characterization.R library fixtures.rds output.rds")
.libPaths(c(normalizePath(args[1L]), .libPaths()))
suppressPackageStartupMessages(library(sceptre))
data.table::setDTthreads(1L)
stopifnot(identical(normalizePath(find.package("sceptre")),
                    file.path(normalizePath(args[1L]), "sceptre")))
stopifnot(identical(digest::digest(file = args[2L], algo = "sha256"),
    "f53c93eabbb30a53d5fed01788ee8ef0ef4447b94ed308540410e21baa624e2d"))
fixtures <- readRDS(args[2L])
stopifnot(length(fixtures) == 28L)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(94201L)
for (k in seq_len(200L)) {
    m <- methods::new("dgTMatrix", i = sample.int(4L, 30L, TRUE) - 1L,
        j = sample.int(3L, 30L, TRUE) - 1L,
        x = sample(c(0, .1, .2, .3, 1, 2, 2^53, 1e-20), 30L, TRUE),
        Dim = c(4L, 3L), Dimnames = list(paste0("g", 1:4), c("z", "a", "m")))
    d <- data.frame(grna_id = paste0("g", sample.int(4L, 8L, TRUE)),
        grna_target = "t", vector_id = sample(c("vZ", "vA"), 8L, TRUE))
    fixtures[[paste0("rounding-", k)]] <- list(m = m, d = d)
}
stopifnot(length(fixtures) == 228L)
generated_input_sha256 <- digest::digest(fixtures, algo = "sha256")
fn <- getFromNamespace("collapse_grna_matrix", "sceptre")
outputs <- lapply(fixtures, function(input) {
    input <- unserialize(serialize(input, NULL))
    before <- serialize(input, NULL)
    result <- fn(input$m, input$d)
    stopifnot(identical(before, serialize(input, NULL)))
    result
})
saveRDS(list(outputs = outputs, generated_input_sha256 = generated_input_sha256,
    fixture_sha256 = digest::digest(file = args[2L], algo = "sha256"),
    package_path = find.package("sceptre"), R = R.version.string,
    DLL_sha256 = digest::digest(file = getLoadedDLLs()[["sceptre"]][["path"]], algo = "sha256"),
    count = length(outputs), input_nonmutation = TRUE), args[3L])
cat(length(outputs), "installed-namespace cases; inputs unchanged.\n")
