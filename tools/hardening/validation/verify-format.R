# Compare mechanical formatting with the committed implementation, allowing only
# single-expression braces. Air inserts these around function/control bodies.
# All operators, literals, arguments, order and other expression nodes must match.
reference <- '74a421a'
formatted <- 'f9f8142bd404c6da1dc39eca2c993f78cb6e5eb1'
normalize_braces <- function(x) {
    if (rlang::is_missing(x)) return(rlang::missing_arg())
    if (is.call(x)) {
        if (identical(x[[1L]], as.name('{')) && length(x) == 2L) {
            return(normalize_braces(x[[2L]]))
        }
        return(as.call(lapply(as.list(x), normalize_braces)))
    }
    if (is.expression(x)) return(as.expression(lapply(x, normalize_braces)))
    if (is.pairlist(x)) return(as.pairlist(lapply(x, normalize_braces)))
    x
}
files <- list.files('R', pattern = '[.]R$', full.names = TRUE)
before <- lapply(files, function(file) parse(text = system2('git',
    c('show', paste0(reference, ':', file)), stdout = TRUE), keep.source = FALSE))
after <- lapply(files, function(file) parse(text = system2('git',
    c('show', paste0(formatted, ':', file)), stdout = TRUE), keep.source = FALSE))
raw_equal <- mapply(identical, before, after)
normalized_equal <- mapply(function(a, b) identical(normalize_braces(a), normalize_braces(b)),
    before, after)
stopifnot(all(normalized_equal))
jsonlite::write_json(list(formatter = 'air 0.9.0', reference = reference, formatted_revision = formatted,
    files = length(files), raw_ast_changed_files = files[!raw_equal],
    normalization = 'single-expression brace blocks only',
    normalized_ast_identical = all(normalized_equal)),
    'tools/hardening/validation/format-ast-proof.json', pretty = TRUE, auto_unbox = TRUE)
cat('Every parsed expression matches after single-expression brace normalization\n')
