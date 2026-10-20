# Use in a fresh process with the intended checker first in R_LIBS.
# Args: extracted package root, output CSV. Retain duplicate function names.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[[1L]])
inventory <- list()
for (path in sort(list.files(file.path(root, 'R'), pattern='\\.[Rr]$', full.names=TRUE))) {
    parsed <- getParseData(parse(path, keep.source=TRUE))
    functions <- BiocCheck:::getFunctionLengths(parsed)
    if (!length(names(functions))) next
    for (i in seq_along(functions)) {
        entry <- functions[[i]]
        inventory[[length(inventory) + 1L]] <- data.frame(
            file=paste0('R/', basename(path)), name=names(functions)[i],
            as.list(entry), stringsAsFactors=FALSE
        )
    }
}
write.csv(do.call(rbind, inventory), args[[2L]], row.names=FALSE)
