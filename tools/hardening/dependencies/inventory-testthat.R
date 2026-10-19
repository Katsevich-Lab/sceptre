# Run from the repository root. Record every syntactic testthat call, including
# named arguments, regardless of whether a conditional branch executes locally.
exports <- getNamespaceExports("testthat")
files <- sort(c(list.files("tests", "\\.[Rr]$", recursive = TRUE, full.names = TRUE),
    list.files("tools", "\\.[Rr]$", recursive = TRUE, full.names = TRUE)))
# Detect future expectation/skip APIs even when the running floor does not
# export them. Exclude project-defined scientific assertion helpers.
parsed <- lapply(files, parse)
project_helpers <- character()
find_helpers <- function(x) {
    if (!is.call(x) && !is.expression(x) && !is.pairlist(x)) return(invisible(NULL))
    if (is.call(x) && length(x) >= 3L && is.symbol(x[[1L]]) &&
        as.character(x[[1L]]) %in% c("<-", "=") && is.symbol(x[[2L]]) &&
        is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("function"))) {
        project_helpers <<- union(project_helpers, as.character(x[[2L]]))
    }
    for (part in as.list(x)) {
        if (missing(part)) next
        if (!is.symbol(part)) find_helpers(part)
    }
}
for (expressions in parsed) find_helpers(expressions)
rows <- list()
walk <- function(x, file) {
    if (!is.call(x) && !is.expression(x) && !is.pairlist(x)) return(invisible(NULL))
    if (is.call(x)) {
        head <- x[[1L]]
        api <- NULL
        if (is.symbol(head)) {
            name <- as.character(head)
            if (name %in% exports || (grepl("^(expect_|skip_|local_|with_|capture_)", name) &&
                !name %in% project_helpers)) api <- name
        }
        if (is.call(head) && identical(head[[1L]], as.name("::")) &&
            identical(as.character(head[[2L]]), "testthat")) api <- as.character(head[[3L]])
        if (length(api)) {
            args <- names(as.list(x)[-1L]); if (is.null(args)) args <- character()
            rows[[length(rows) + 1L]] <<- data.frame(file = file, api = api,
                arguments = paste(args[nzchar(args)], collapse = ","),
                call = substr(paste(deparse(x, width.cutoff = 500L), collapse = " "), 1L, 300L))
        }
    }
    for (part in as.list(x)) {
        if (missing(part)) next
        if (!is.symbol(part)) walk(part, file)
    }
}
for (i in seq_along(files)) walk(parsed[[i]], files[[i]])
result <- unique(do.call(rbind, rows))
write.table(result, "tools/hardening/dependencies/testthat-calls.tsv", sep = "\t", quote = TRUE, qmethod = "double", row.names = FALSE)
writeLines(sort(unique(result$api)), "tools/hardening/dependencies/testthat-apis.txt")
cat(nrow(result), "unique call sites;", length(unique(result$api)), "APIs\n")
