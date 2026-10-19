# Run in a fresh R process with R_LIBS_USER pointing at the isolated floor lib.
# Never load latest testthat then change .libPaths: namespaces are process-wide.
stopifnot(identical(as.character(packageVersion("testthat")), "3.1.7"))
api <- readLines("tools/hardening/dependencies/testthat-apis.txt")
stopifnot(all(api %in% getNamespaceExports("testthat")))
forms <- lapply(api, function(name) names(formals(getExportedValue("testthat", name))))
names(forms) <- api
out <- Sys.getenv("SCEPTRE_DEPENDENCY_EVIDENCE", unset = tempdir())
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(forms, file.path(out, "testthat-floor-formals.rds"))
# Exercise newer unconditional APIs and argument forms before package tests.
# A real installed package binding is required for the mocking contract.
f <- function() {
    original <- stats::runif
    testthat::local_mocked_bindings(runif = function(n, ...) rep(0.25, n),
        .package = "stats")
    testthat::expect_identical(stats::runif(2L), c(0.25, 0.25))
    testthat::expect_no_error(1 + 1)
    testthat::expect_no_warning(1 + 1)
    testthat::expect_error(stop("controlled floor error"), regexp = "controlled", fixed = TRUE)
    testthat::expect_warning(warning("controlled floor warning"), regexp = "controlled", fixed = TRUE)
    testthat::expect_error(stop("UPPER"), regex = "upper", ignore.case = TRUE, info = "regexp alias")
    testthat::expect_match("literal.+", "literal.+", fixed = TRUE)
    testthat::expect_equal(1, 1 + 1e-12, tolerance = 1e-8, info = "floor tolerance")
    original
}
original <- f()
testthat::expect_identical(stats::runif, original)
writeLines(c(paste("testthat", packageVersion("testthat")),
    paste("library", find.package("testthat")), capture.output(sessionInfo())),
    file.path(out, "floor-probe-session.txt"))
cat("All", length(api), "inventory APIs exported; focused API behavior passed.\n")
