args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
entries <- jsonlite::fromJSON(args[1L], simplifyVector = FALSE)
for (entry in entries) {
    decode <- function(x) {
        value <- parse(text = paste0('"', x, '"'))
        stopifnot(length(value) == 1L, is.character(value[[1L]]))
        value[[1L]]
    }
    original <- decode(entry$raw)
    pieces <- lapply(entry$chunks, decode)
    stopifnot(identical(original, do.call(paste0, pieces)))
    if (entry$head %in% c('stop', 'warning', 'message')) {
        capture <- function(args) {
            tryCatch(do.call(entry$head, args),
                error = function(x) conditionMessage(x),
                warning = function(x) conditionMessage(x),
                message = function(x) conditionMessage(x))
        }
        stopifnot(identical(capture(list(original)), capture(pieces)))
    }
}
cat('PASS ', length(entries),
    ' split literals preserve evaluated text and signal conditionMessage output\n', sep = '')
