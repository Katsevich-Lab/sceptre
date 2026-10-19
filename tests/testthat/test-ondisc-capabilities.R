compatible_ondisc_snapshot <- function() list(
    formals = ondisc_required_api(),
    exports = "initialize_odm_from_backing_file",
    slots = c(h5_file = "character", dimension = "integer", feature_ids = "character",
        ptr = "externalptr", integer_id = "integer"),
    methods = c(row_read = TRUE, dim = TRUE, dimnames = TRUE),
    native = ondisc_required_native(), version = "999.0")

test_that("capability acceptance examines API rather than a version whitelist", {
    snapshot <- compatible_ondisc_snapshot()
    expect_identical(ondisc_capability_issues(snapshot), character())
    snapshot$version <- "0.0.0"
    expect_identical(ondisc_capability_issues(snapshot), character())
    snapshot$formals$create_odm_from_r_matrix_internal <- NULL
    expect_match(ondisc_capability_issues(snapshot)[[1L]], "missing function.*create_odm_from_r_matrix_internal")
})

test_that("changed class, exported reader or native ABI fails explicitly", {
    snapshot <- compatible_ondisc_snapshot()
    snapshot$slots[["ptr"]] <- "integer"
    expect_true(any(grepl("incompatible `odm` slots", ondisc_capability_issues(snapshot))))
    snapshot <- compatible_ondisc_snapshot()
    snapshot$exports <- character()
    expect_true(any(grepl("not exported", ondisc_capability_issues(snapshot))))
    snapshot <- compatible_ondisc_snapshot()
    snapshot$formals$initialize_odm_from_backing_file <- c("filename", "...")
    expect_true(any(grepl("incompatible arguments.*initialize_odm", ondisc_capability_issues(snapshot))))
    snapshot <- compatible_ondisc_snapshot()
    snapshot$methods[["row_read"]] <- FALSE
    expect_true(any(grepl("missing `odm` method for row_read", ondisc_capability_issues(snapshot))))
    snapshot <- compatible_ondisc_snapshot()
    snapshot$native[["_ondisc_compute_n_ok_pairs_ondisc"]] <- 11L
    expect_true(any(grepl("registered native routine.*compute_n_ok_pairs", ondisc_capability_issues(snapshot))))
})

test_that("missing optional installation has an actionable backed-only error", {
    local_mocked_bindings(ondisc_namespace = function() NULL)
    expect_error(require_ondisc_capabilities(), "package is not installed.*qualified source.*restart R.*In-memory")
    sparse <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 2, dims = c(1L, 2L), repr = "R")
    expect_identical(count_store_value(as_count_store(sparse)), sparse)
    expect_identical(storage_read_response(sparse, 1L), c(2, 0))
})

test_that("an incompatible backend is refused before count files are created", {
    snapshot <- compatible_ondisc_snapshot()
    snapshot$formals$create_odm_from_r_matrix_internal <- c("mat", "file_to_write")
    local_mocked_bindings(ondisc_namespace = function() new.env(parent = emptyenv()),
        ondisc_capability_snapshot = function(namespace) snapshot)
    directory <- tempfile("incompatible-ondisc-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    expect_error(storage_create_backed_pair(NULL, NULL, directory, 12L),
        "version 999.0.*incompatible arguments.*create_odm.*qualified source")
    expect_identical(list.files(directory), character())
})

test_that("the installed optional backend exposes the required operations", {
    skip_if_not_installed("ondisc")
    snapshot <- ondisc_capability_snapshot(asNamespace("ondisc"))
    expect_identical(ondisc_capability_issues(snapshot), character())
    expect_no_error(require_ondisc_capabilities())
})
