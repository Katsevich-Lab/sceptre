#' Internal backed count store
#'
#' An anchored backing descriptor and an optional process-local handle.
#' Dimensions and labels are available from the descriptor while closed.
#' @name BackedCountStore-class
#' @aliases dim,BackedCountStore-method dimnames,BackedCountStore-method
#' @keywords internal
NULL

# The package stores a small typed descriptor and a process-local handle; ondisc
# remains optional and is loaded only to reconstruct a live numerical view.
#' @noRd
as_count_store <- function(x, descriptor = NULL) {
    if (methods::is(x, "dgRMatrix")) return(x)
    if (methods::is(x, "BackedCountStore")) {
        if (!is.null(descriptor) && !identical(x@descriptor, descriptor)) {
            stop("Count store descriptor differs from its authoritative descriptor.")
        }
        return(x)
    }
    if (is.null(x)) {
        if (is.null(descriptor)) return(NULL)
        return(methods::new("BackedCountStore", descriptor = descriptor, handle = NULL))
    }
    if (inherits(x, "odm")) {
        if (is.null(descriptor)) descriptor <- backed_store_descriptor(x) else {
            validate_backed_descriptor(x, descriptor)
            descriptor$path <- normalizePath(attr(x, "h5_file"), mustWork = TRUE)
        }
        pointer <- attr(x, "ptr")
        if (identical(pointer, methods::new("externalptr"))) pointer <- NULL
        return(methods::new("BackedCountStore", descriptor = descriptor, handle = pointer))
    }
    stop("Count stores must be normalized row-accessible Matrix or ondisc objects.")
}

#' @noRd
count_store_value <- function(x) {
    if (is.null(x)) stop("Count store is closed or unavailable; reopen its backing files.")
    if (methods::is(x, "dgRMatrix")) return(x)
    if (!methods::is(x, "BackedCountStore")) stop("Unrecognized count store.")
    if (is.null(x@handle) || identical(x@handle, methods::new("externalptr"))) {
        stop("Count store is closed; reopen its backing files.")
    }
    require_ondisc_capabilities()
    d <- x@descriptor
    methods::new("odm", h5_file = d$path, dimension = as.integer(d$dimensions),
        feature_ids = as.character(d$dimnames[[1L]]), integer_id = as.integer(d$integer_id),
        ptr = x@handle)
}

#' @noRd
count_store_descriptor <- function(x) {
    if (methods::is(x, "BackedCountStore")) x@descriptor else NULL
}

#' @noRd
close_count_store <- function(x) {
    if (!methods::is(x, "BackedCountStore")) stop("Expected a backed count store.")
    x@handle <- NULL
    x
}

setMethod("dim", "BackedCountStore", function(x) x@descriptor$dimensions)
setMethod("dimnames", "BackedCountStore", function(x) x@descriptor$dimnames)
setValidity("BackedCountStore", function(object) {
    d <- object@descriptor
    if (!all(c("path", "content_md5", "dimensions", "dimnames", "integer_id") %in% names(d)) ||
        !is.character(d$path) || length(d$path) != 1L || is.na(d$path) || !nzchar(d$path) ||
        !is.character(d$content_md5) || length(d$content_md5) != 1L ||
        is.na(d$content_md5) || !grepl("^[[:xdigit:]]{32}$", d$content_md5) ||
        !is.integer(d$dimensions) || length(d$dimensions) != 2L || anyNA(d$dimensions) || any(d$dimensions < 0L) ||
        !is.list(d$dimnames) || length(d$dimnames) != 2L ||
        !is.integer(d$integer_id) || length(d$integer_id) != 1L || is.na(d$integer_id)) {
        return("Invalid backed count-store descriptor")
    }
    for (axis in seq_len(2L)) {
        labels <- d$dimnames[[axis]]
        if (!is.null(labels) && (!is.character(labels) || length(labels) != d$dimensions[axis])) {
            return("Backed count-store labels must align with descriptor dimensions")
        }
    }
    if (!is.null(object@handle) && !methods::is(object@handle, "externalptr")) {
        return("Backed count-store handle must be an external pointer or NULL")
    }
    TRUE
})
