pkgload::load_all(".", recompile = FALSE)
testthat::test_local(
    ".",
    filter = paste0(
        "^(distributed-execution|perturbation-model|schema-migration|",
        "run-plan-migration|fit-cache-migration|refactor-context-execution)$"
    ),
    reporter = "summary", stop_on_failure = TRUE
)
for (name in c("prepare_execution", "execute_work_unit", "collect_execution")) {
    cat("\nEXAMPLE", name, "\n")
    path <- tempfile(fileext = ".R")
    tools::Rd2ex(tools::parse_Rd(file.path("man", paste0(name, ".Rd"))), out = path)
    sys.source(path, envir = new.env(parent = globalenv()))
    unlink(path)
}
for (file in c("R/execution.R", "R/distributed_execution.R",
               "R/schema_migration.R", "R/perturbation_fits.R")) {
    lengths <- BiocCheck:::getFunctionLengths(
        utils::getParseData(parse(file, keep.source = TRUE))
    )
    stopifnot(all(vapply(lengths, function(x) x[["length"]] <= 50, logical(1))))
    cat(file, max(vapply(lengths, function(x) x[["length"]], numeric(1))), "lines\n")
}
print(sessionInfo())
