# Run in a fresh R process from the repository root.
parity_library <- Sys.getenv("SCEPTRE_PARITY_LIBRARY")
if (nzchar(parity_library)) {
    library(sceptre, lib.loc = parity_library, character.only = FALSE)
} else {
    pkgload::load_all(".", recompile = FALSE)
}
source("tests/testthat/helper-refactor.R")
validation_dir <- Sys.getenv("SCEPTRE_PARITY_OUTPUT",
    file.path("tools", "hardening", "validation"))
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
inputs <- readRDS("tests/testthat/fixtures/refactor/inputs.rds")
cases <- readRDS("tests/testthat/fixtures/refactor/cases.rds")
outputs <- lapply(cases, function(case) refactor_run_case(inputs, case))
saveRDS(outputs, file.path(validation_dir, "parity-current.rds"))
jsonlite::write_json(list(
    mode = if (nzchar(parity_library)) "installed" else "source",
    package_path = getNamespaceInfo(asNamespace("sceptre"), "path"),
    package_version = as.character(packageVersion("sceptre")),
    source_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    archive_sha256 = Sys.getenv("SCEPTRE_PARITY_ARCHIVE_SHA256", unset = NA_character_),
    session = capture.output(sessionInfo())
), file.path(validation_dir, "parity-environment.json"), pretty = TRUE, auto_unbox = TRUE)
source(file.path("tools", "hardening", "validation", "compare-parity.R"))
