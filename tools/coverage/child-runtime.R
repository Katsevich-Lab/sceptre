# Called by covr's command child with the same isolated-library configuration.
# This is a direct receipt for that child, after the separate installed test run.
write_coverage_child_runtime <- function(expected_sceptre_path, receipt_file) {
    sceptre_ns <- asNamespace("sceptre")
    sceptre_path <- normalizePath(getNamespaceInfo(sceptre_ns, "path"), mustWork = TRUE)
    expected_sceptre_path <- normalizePath(expected_sceptre_path, mustWork = TRUE)
    backend <- NULL
    if (requireNamespace("ondisc", quietly = TRUE)) {
        path <- normalizePath(getNamespaceInfo(asNamespace("ondisc"), "path"), mustWork = TRUE)
        desc <- utils::packageDescription("ondisc", lib.loc = dirname(path))
        dll <- getLoadedDLLs()[["ondisc"]]
        dll_path <- if (is.null(dll)) NULL else dll[["path"]]
        backend <- list(path = path, version = desc$Version, RemoteSha = desc$RemoteSha,
            description_md5 = unname(tools::md5sum(file.path(path, "DESCRIPTION"))),
            dll_path = dll_path,
            dll_md5 = if (is.null(dll_path)) NULL else unname(tools::md5sum(dll_path)))
    }
    expected_sha <- Sys.getenv("SCEPTRE_EXPECT_ONDISC_SHA")
    sha_matches <- if (nzchar(expected_sha)) identical(backend$RemoteSha, expected_sha) else NULL
    receipt <- list(scope = "Actual covr command child after installed test run; same isolated library configuration, not a historical test-process namespace capture",
        pid = Sys.getpid(), R = R.version.string, platform = R.version$platform,
        library_paths = .libPaths(), sceptre_path = sceptre_path,
        expected_sceptre_path = expected_sceptre_path,
        sceptre_path_matches = identical(sceptre_path, expected_sceptre_path),
        ondisc = backend, expected_ondisc_sha = expected_sha,
        expected_ondisc_sha_matches = sha_matches)
    # Preserve the actual identity even when the subsequent assertion rejects it.
    jsonlite::write_json(receipt, receipt_file, pretty = TRUE, auto_unbox = TRUE, null = "null")
    if (!receipt$sceptre_path_matches) stop("Coverage child loaded sceptre from the wrong installation")
    if (nzchar(expected_sha) && !isTRUE(sha_matches))
        stop("Coverage child ondisc RemoteSha does not match SCEPTRE_EXPECT_ONDISC_SHA")
    invisible(receipt)
}
