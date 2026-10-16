#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply a fresh fixture directory.")
output <- args[[1L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(output))) stop("Triplet fixtures must be generated in a fresh directory.")
data.table::setDTthreads(1L)
invisible(loadNamespace("Matrix"))
hash <- function(path) digest::digest(file = path, algo = "sha256")
slot_hash <- function(object) {
    path <- tempfile(tmpdir = output)
    connection <- file(path, "wb")
    serialize(object, connection, version = 3)
    close(connection)
    on.exit(unlink(path))
    hash(path)
}
manifest <- list(recipe = "approved triplet-resource-plan.md", R = R.version.string,
    dimensions = c(20000L, 50000L), fixtures = list())
for (nnz in c(1000000L, 8000000L)) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(271828L)
    k <- as.double(seq_len(nnz))
    q <- (104729 * k + 37) %% 1000000000
    stopifnot(!anyDuplicated(q))
    permutation <- sample.int(nnz, nnz, replace = FALSE)
    permuted_k <- k[permutation]
    coordinates <- q[permutation]
    i <- as.integer(coordinates %% 20000)
    j <- as.integer(coordinates %/% 20000)
    stopifnot(all(i >= 0L & i < 20000L), all(j >= 0L & j < 50000L))
    dimnames <- list(sprintf("r%05d", seq_len(20000L)), sprintf("c%05d", seq_len(50000L)))
    for (representation in c("numeric", "logical")) {
        values <- if (representation == "numeric") {
            x <- ((permuted_k %% 97) + 1) / 7
            x[permuted_k %% 19 == 0] <- 0
            x
        } else permuted_k %% 7 != 0
        object <- methods::new(if (representation == "numeric") "dgTMatrix" else "lgTMatrix",
            i = i, j = j, x = values, Dim = c(20000L, 50000L), Dimnames = dimnames)
        methods::validObject(object)
        name <- paste(nnz, representation, sep = "-")
        path <- file.path(output, paste0(name, ".rds"))
        saveRDS(object, path, compress = FALSE, version = 3)
        manifest$fixtures[[name]] <- list(file = basename(path), sha256 = hash(path),
            class = class(object)[[1L]], nnz = nnz, density = nnz / 1000000000,
            zero_or_false = sum(!values),
            i_range = range(i), j_range = range(j), values_range = range(values),
            slot_sha256 = lapply(list(i = i, j = j, x = values, Dim = object@Dim,
                Dimnames = dimnames), slot_hash), RNGkind = RNGkind(), seed = 271828L)
        rm(object, values)
        if (exists("x", inherits = FALSE)) rm(x)
        gc()
    }
    rm(k, q, permutation, permuted_k, coordinates, i, j, dimnames)
    gc()
}
jsonlite::write_json(manifest, file.path(output, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
