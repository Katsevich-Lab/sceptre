# ondisc is optional. Inspect its actual API only on a backed path, before file
# creation/reopening or native execution. Version text is diagnostic, not proof
# of private-API or numerical compatibility.
ondisc_capability_cache <- new.env(parent = emptyenv())

ondisc_required_api <- function() {
    prefix <- c(
        "file_name_in",
        "f_row_ptr",
        "n_genes",
        "n_cells_orig",
        "n_cells_sub",
        "grna_group_idxs"
    )
    suffix <- c(
        "all_nt_idxs",
        "to_analyze_response_idxs",
        "to_analyze_grna_idxs",
        "control_group_complement",
        "cells_in_use"
    )
    list(
        create_odm_from_r_matrix_internal = c(
            "mat",
            "file_to_write",
            "chunk_size",
            "compression_level",
            "integer_id"
        ),
        initialize_odm_from_backing_file = "odm_file",
        threshold_count_matrix_cpp = c(
            "file_name_in",
            "f_row_ptr",
            "row_idx",
            "threshold"
        ),
        compute_n_ok_pairs_ondisc = c(prefix, suffix, "unique_response_idxs"),
        compute_nt_nonzero_matrix_and_n_ok_pairs_ondisc = c(
            prefix,
            "indiv_nt_grna_idxs",
            suffix
        ),
        compute_n_trt_cells_matrix_ondisc = c(
            "file_name_in",
            "f_row_ptr",
            "n_cells_orig",
            "n_cells_sub",
            "n_genes",
            "nt_grna_group_idxs",
            "cells_in_use"
        )
    )
}

ondisc_required_native <- function() {
    c(
        `_ondisc_create_odm_r_matrix_cpp` = 10L,
        `_ondisc_load_row_cpp` = 4L,
        `_ondisc_read_integer_vector` = 3L,
        `_ondisc_read_feature_ids` = 2L,
        `_ondisc_read_row_ptr` = 2L,
        `_ondisc_threshold_count_matrix_cpp` = 4L,
        `_ondisc_compute_n_ok_pairs_ondisc` = 12L,
        `_ondisc_compute_nt_nonzero_matrix_and_n_ok_pairs_ondisc` = 12L,
        `_ondisc_compute_n_trt_cells_matrix_ondisc` = 7L
    )
}

ondisc_capability_snapshot <- function(namespace) {
    functions <- lapply(names(ondisc_required_api()), function(name) {
        value <- get0(name, envir = namespace, inherits = FALSE)
        if (is.function(value)) names(formals(value)) else NULL
    })
    names(functions) <- names(ondisc_required_api())
    definition <- methods::getClassDef("odm", where = namespace)
    slots <- if (is.null(definition)) {
        character()
    } else {
        vapply(
            definition@slots,
            function(type) as.character(type)[[1L]],
            character(1)
        )
    }
    dll <- getLoadedDLLs()[["ondisc"]]
    native <- if (is.null(dll)) list() else getDLLRegisteredRoutines(dll)$.Call
    has_method <- function(generic, signature) {
        method <- methods::selectMethod(
            generic,
            signature = signature,
            optional = TRUE
        )
        methods::is(method, "MethodDefinition") &&
            identical(as.character(method@defined[[1L]]), "odm")
    }
    methods <- c(
        row_read = has_method(
            "[",
            c(x = "odm", i = "ANY", j = "missing", drop = "missing")
        ),
        dim = has_method("dim", "odm"),
        dimnames = has_method("dimnames", "odm")
    )
    list(
        formals = functions,
        exports = getNamespaceExports(namespace),
        slots = slots,
        methods = methods,
        native = vapply(
            native,
            function(routine) as.integer(routine$numParameters),
            integer(1)
        ),
        version = as.character(getNamespaceVersion(namespace))
    )
}

# Pure validation makes incompatible-install tests independent of the user's
# installed library; it does not call suspect native code to discover an ABI.
ondisc_capability_issues <- function(snapshot) {
    issues <- ondisc_api_argument_issues(snapshot)
    if (!"initialize_odm_from_backing_file" %in% snapshot$exports) {
        issues <- c(
            issues,
            "`initialize_odm_from_backing_file()` is not exported"
        )
    }
    expected_slots <- c(
        h5_file = "character",
        dimension = "integer",
        feature_ids = "character",
        ptr = "externalptr",
        integer_id = "integer"
    )
    if (!identical(snapshot$slots[names(expected_slots)], expected_slots)) {
        issues <- c(
            issues,
            paste0(
                "incompatible `odm` slots (path, dimensions, feature IDs, ",
                "handle or shared import ID)"
            )
        )
    }
    for (method in c("row_read", "dim", "dimnames")) {
        if (!isTRUE(unname(snapshot$methods[method]))) {
            issues <- c(issues, paste0("missing `odm` method for ", method))
        }
    }
    expected_native <- ondisc_required_native()
    for (name in names(expected_native)) {
        actual <- snapshot$native[name]
        if (
            length(actual) != 1L ||
                is.na(actual) ||
                actual != expected_native[[name]]
        ) {
            issues <- c(
                issues,
                paste0(
                    "missing or incompatible registered native routine `",
                    name,
                    "`"
                )
            )
        }
    }
    issues
}

ondisc_capability_error <- function(issues, version = "unavailable") {
    stop(
        "The optional `ondisc` installation (version ",
        version,
        ") cannot support SCEPTRE backed counts: ",
        paste(issues, collapse = "; "),
        ". Install the qualified source ",
        "`timothy-barry/ondisc@589803edea2a66ad15e7f07186bee16c7739ab8c`",
        " and restart R. See the optional disk-backed setup in SCEPTRE's ",
        "README.",
        " In-memory analyses do not require ondisc.",
        call. = FALSE
    )
}

ondisc_namespace <- function() {
    if (!requireNamespace("ondisc", quietly = TRUE)) {
        return(NULL)
    }
    asNamespace("ondisc")
}

require_ondisc_capabilities <- function() {
    namespace <- ondisc_namespace()
    if (is.null(namespace)) {
        ondisc_capability_error(
            "package is not installed or its namespace could not be loaded"
        )
    }
    if (identical(ondisc_capability_cache$namespace, namespace)) {
        return(invisible(NULL))
    }
    snapshot <- tryCatch(
        ondisc_capability_snapshot(namespace),
        error = identity
    )
    if (inherits(snapshot, "error")) {
        ondisc_capability_error(
            paste("API inspection failed:", conditionMessage(snapshot)),
            as.character(getNamespaceVersion(namespace))
        )
    }
    issues <- ondisc_capability_issues(snapshot)
    if (length(issues)) {
        ondisc_capability_error(issues, snapshot$version)
    }
    ondisc_capability_cache$namespace <- namespace
    invisible(NULL)
}

ondisc_api_argument_issues <- function(snapshot) {
    issues <- character()
    for (name in names(ondisc_required_api())) {
        actual <- snapshot$formals[[name]]
        if (is.null(actual)) {
            issues <- c(issues, paste0("missing function `", name, "()`"))
        } else if (!identical(actual, ondisc_required_api()[[name]])) {
            issues <- c(
                issues,
                paste0("incompatible arguments for `", name, "()`")
            )
        }
    }
    issues
}
