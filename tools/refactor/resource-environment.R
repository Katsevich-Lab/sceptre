#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: resource-environment.R <library> <output.json>")
.libPaths(c(normalizePath(args[[1L]]), .libPaths()))
packages <- c("sceptre", "ondisc", "Matrix", "Rcpp", "BH", "data.table", "dplyr", "parallelly", "withr", "jsonlite", "digest", "sceptredata", "MASS", "stats")
information <- stats::setNames(lapply(packages, function(name) {
    path <- find.package(name)
    dlls <- list.files(file.path(path, "libs"), pattern = "\\.(so|dylib|dll)$", recursive = TRUE, full.names = TRUE)
    list(version = as.character(packageVersion(name)), path = path,
         description_md5 = unname(tools::md5sum(file.path(path, "DESCRIPTION"))),
         binary_md5 = stats::setNames(as.list(unname(tools::md5sum(dlls))), basename(dlls)))
}), packages)
jsonlite::write_json(list(R = R.version.string, external = extSoftVersion(),
    platform = R.version$platform, packages = information,
    note = "Dependency resolution for the frozen package library; package source commits are recorded in the run plan. Sceptre code/path differs deliberately."),
    args[[2L]], pretty = TRUE, auto_unbox = TRUE)
