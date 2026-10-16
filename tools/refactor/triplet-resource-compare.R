#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Supply original output candidate output comparison JSON.")
invisible(loadNamespace("Matrix"))
original <- readRDS(args[[1L]])
candidate <- readRDS(args[[2L]])
stopifnot(identical(class(candidate), class(original)),
          identical(dim(candidate), dim(original)),
          identical(dimnames(candidate), dimnames(original)))
slots <- c("p", if (methods::is(candidate, "dgRMatrix")) "j" else "i", "x")
same <- vapply(slots, function(name) identical(methods::slot(original, name),
                                             methods::slot(candidate, name)), logical(1))
stopifnot(all(same))
jsonlite::write_json(list(exact_class_dimensions_dimnames = TRUE, exact_slots = same,
    stored = length(candidate@x), zeros = sum(candidate@x == 0),
    fractions = sum(candidate@x != floor(candidate@x))), args[[3L]],
    pretty = TRUE, auto_unbox = TRUE)
