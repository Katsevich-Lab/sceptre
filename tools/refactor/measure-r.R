#!/usr/bin/env Rscript
# Static parsing only: never source or evaluate the measured package.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript measure-r.R SOURCE_ROOT OUTPUT_DIR")
root <- normalizePath(args[[1L]], mustWork = TRUE)
dir.create(args[[2L]], recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(args[[2L]], mustWork = TRUE)
stopifnot(as.character(packageVersion("cyclocomp")) == "1.1.2")
sha <- function(x) digest::digest(x, algo = "sha256", serialize = FALSE)
head_name <- function(x) if (is.call(x) && is.symbol(x[[1L]])) as.character(x[[1L]]) else ""
# Missing arguments can occur in calls such as x[, i]. Do not evaluate them.
is_missing_part <- function(parts, i) identical(parts[i], alist(.missing = ) |> unname())
prune_nested <- function(x) {
    if (!is.call(x)) return(x)
    if (head_name(x) == "function") return(quote(.nested_function_literal))
    parts <- as.list(x)
    for (i in seq_along(parts)) {
        if (!is_missing_part(parts, i)) parts[i] <- list(prune_nested(parts[[i]]))
    }
    as.call(parts)
}
control_depth <- function(x, depth = 0L) {
    if (!is.call(x) || head_name(x) == "function") return(depth)
    next_depth <- depth + as.integer(head_name(x) %in% c("if", "for", "while", "repeat", "switch"))
    parts <- as.list(x)
    result <- next_depth
    for (i in seq_along(parts)) {
        if (!is_missing_part(parts, i)) result <- max(result, control_depth(parts[[i]], next_depth))
    }
    result
}
wrapped_body <- function(body) as.call(list(as.name("function"), pairlist(), body))
# Small analyzer-contract checks, not tests of package implementation.
probe <- quote(function() { if (x) y; f <- function() { if (z) y } })
stopifnot(cyclocomp::cyclocomp(wrapped_body(prune_nested(probe[[3L]]))) == 2L,
          cyclocomp::cyclocomp(wrapped_body(probe[[3L]])) == 3L,
          control_depth(quote(if (x) while (y) z)) == 2L)
records <- list()
file_records <- list()
files <- sort(list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE))
files <- files[basename(files) != "RcppExports.R"]
for (file in files) {
    rel <- substring(file, nchar(root) + 2L)
    parsed <- parse(file, keep.source = TRUE)
    pd <- utils::getParseData(parsed)
    text_lines <- readLines(file, warn = FALSE)
    tokens <- pd[pd$terminal & pd$token != "COMMENT", ]
    code_lines <- if (nrow(tokens)) unique(unlist(Map(seq.int, tokens$line1, tokens$line2))) else integer()
    file_records[[length(file_records) + 1L]] <- list(path = rel,
        sha256 = digest::digest(file = file, algo = "sha256", serialize = FALSE),
        physical_lines = length(text_lines), noncomment_token_lines = length(code_lines))
    declarations <- pd[pd$token == "FUNCTION", ]
    declarations <- declarations[order(declarations$line1, declarations$col1), ]
    local_records <- list()
    for (j in seq_len(nrow(declarations))) {
        token <- declarations[j, ]
        fnrow <- pd[pd$id == token$parent, ]
        source <- utils::getParseText(pd, fnrow$id)
        fn <- parse(text = source, keep.source = FALSE)[[1L]]
        stopifnot(head_name(fn) == "function")
        siblings <- pd[pd$parent == fnrow$parent, ]
        siblings <- siblings[order(siblings$line1, siblings$col1), ]
        label <- NULL
        kind <- "anonymous"
        if (any(siblings$token %in% c("LEFT_ASSIGN", "EQ_ASSIGN"))) {
            lhs <- siblings[siblings$id != fnrow$id & siblings$token == "expr", ]
            if (nrow(lhs) == 1L) {
                label <- utils::getParseText(pd, lhs$id)
                kind <- "assigned"
            }
        }
        fn_position <- which(siblings$id == fnrow$id)
        if (is.null(label) && length(fn_position) == 1L && fn_position >= 3L &&
            siblings$token[[fn_position - 1L]] == "EQ_SUB") {
            label <- siblings$text[[fn_position - 2L]]
            kind <- "named_argument"
        }
        parent_row <- pd[pd$id == fnrow$parent, ]
        if (is.null(label) && nrow(parent_row)) {
            parent_expr <- tryCatch(parse(text = utils::getParseText(pd, parent_row$id), keep.source = FALSE)[[1L]], error = function(e) NULL)
            if (head_name(parent_expr) == "setMethod" && is.character(parent_expr[[2L]])) {
                label <- paste0("setMethod[", parent_expr[[2L]], "]")
                kind <- "S4_method"
            }
        }
        if (is.null(label)) label <- paste0("<anonymous@", fnrow$line1, ":", fnrow$col1, ">")
        start <- fnrow$line1 * 1000000 + fnrow$col1
        end <- fnrow$line2 * 1000000 + fnrow$col2
        parents <- which(vapply(local_records, function(r) r$position_start < start && r$position_end >= end, logical(1L)))
        parent <- if (length(parents)) tail(parents, 1L) else NULL
        id <- if (is.null(parent)) paste(rel, label, sep = "::") else paste(local_records[[parent]]$id, label, sep = "::")
        own <- prune_nested(fn[[3L]])
        local_records[[length(local_records) + 1L]] <- list(id = id, name = label, file = rel,
            kind = kind, top_level = is.null(parent), parent_id = if (is.null(parent)) NULL else local_records[[parent]]$id,
            start_line = fnrow$line1, start_column = fnrow$col1, end_line = fnrow$line2, end_column = fnrow$col2,
            position_start = start, position_end = end,
            source_sha256 = sha(source), ast_sha256 = sha(paste(deparse(fn, width.cutoff = 500L), collapse = "\n")),
            cyclomatic_exclusive_body = as.integer(cyclocomp::cyclocomp(wrapped_body(own))),
            cyclomatic_inclusive_body = as.integer(cyclocomp::cyclocomp(wrapped_body(fn[[3L]]))),
            control_nesting_exclusive_body = control_depth(own),
            physical_span_lines_including_nested = fnrow$line2 - fnrow$line1 + 1L)
    }
    records <- c(records, local_records)
}
result <- list(schema_version = 1L, language = "R", analyzer = list(name = "cyclocomp", version = as.character(packageVersion("cyclocomp")),
    R = R.version.string, parser = "base::parse(keep.source=TRUE) plus utils::getParseData"),
    rules = list(scope = "All function literals in handwritten R/*.R, including S4 methods and nested/anonymous functions; RcppExports.R excluded. Package source is parsed, never evaluated.",
        primary_complexity = "cyclocomp on function body wrapped in function() with nested function literals replaced by a symbol. Formal default expressions excluded. Each function literal is a separate exclusive unit.",
        inclusive_complexity = "Additional cyclocomp on full function body including nested bodies, for top-level-owner inspection only. Never sum inclusive and child scores.",
        cyclocomp_rule = "Unmodified cyclocomp 1.1.2 control-flow graph E-N+2; if/for/while/repeat/short-circuit && and || plus its return/break/next handling. switch is not specially modeled by this analyzer.",
        nesting = "Custom static AST maximum simultaneous depth of if/for/while/repeat/switch, excluding nested functions and defaults. Braces and boolean operators add zero; else-if is another nested if. A syntactic engineering diagnostic, not cognitive complexity.",
        LOC = "Physical function spans include comments and nested functions and overlap; do not sum spans. File noncomment_token_lines counts unique lines touched by noncomment parser tokens.",
        limitations = "Quoted/unexecuted syntactic function literals and controls are counted. Runtime dispatch, dynamic evaluation, semantic duplication, direct-call coupling, cognitive complexity and human comprehension are unmeasured."),
    files = file_records, functions = records)
jsonlite::write_json(result, file.path(output, "r-metrics.json"), pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
cat(length(records), "R function literals measured\n")
