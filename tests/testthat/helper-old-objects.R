# An object of 0.99.1 and earlier: the matrices in one-element lists and the
# covariates in their own slot, in place of the `data` slot.
old_layout_object <- function(sceptre_object) {
    old <- sceptre_object
    attr(old, "response_matrix") <- list(get_response_matrix(sceptre_object))
    attr(old, "grna_matrix") <- list(get_grna_matrix(sceptre_object))
    attr(old, "covariate_data_frame") <- get_cell_covariates(sceptre_object)
    attr(old, "data") <- NULL
    attr(old, "seed") <- NULL
    attr(old, "version") <- NULL
    old
}


# The library holding the installed sceptre under test, or NULL when sceptre
# is loaded from source (for example, by devtools::load_all()), in which case
# a separate R process could not load the code under test.
installed_sceptre_library <- function() {
    lib <- dirname(system.file(package = "sceptre"))
    if (file.exists(file.path(lib, "sceptre", "Meta", "package.rds"))) {
        lib
    } else {
        NULL
    }
}


# Runs `f(...)` in a fresh R process that has loaded only the installed
# sceptre, and returns its value together with the namespaces loaded before
# `f` ran.
run_in_fresh_session <- function(f, ...) {
    testthat::skip_if_not_installed("callr")
    lib <- installed_sceptre_library()
    testthat::skip_if(is.null(lib), "sceptre is not installed")
    # `f` must not carry the test environment into the other process
    environment(f) <- globalenv()
    callr::r(
        function(lib, f, args) {
            suppressPackageStartupMessages(library(sceptre, lib.loc = lib))
            loaded <- loadedNamespaces()
            list(loaded_before = loaded, value = do.call(f, args))
        },
        args = list(lib = lib, f = f, args = list(...))
    )
}
