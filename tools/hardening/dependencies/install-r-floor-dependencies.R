# Intended for a disposable R 4.1 CI library, not a user's ordinary library.
# This source snapshot predates Matrix's R>=4.4 requirement. Resolve against
# that repository coherently instead of mixing current packages with old R.
stopifnot(getRversion() >= "4.1", getRversion() < "4.2")
repo <- "https://packagemanager.posit.co/cran/2024-03-01"
options(repos = c(CRAN = repo), pkgType = "source")
description <- read.dcf("DESCRIPTION")
fields <- intersect(c("Imports", "LinkingTo", "Suggests"), colnames(description))
deps <- unique(unlist(lapply(description[1L, fields], function(x)
    names(tools:::.split_dependencies(x)))))
# Optional storage/interoperability has its own current-ecosystem lanes.
# The advertised R-floor lane deliberately tests installation without it.
optional_backends <- c("ondisc", "rhdf5", "SingleCellExperiment",
    "SummarizedExperiment", "S4Vectors", "SeuratObject")
deps <- setdiff(deps, c(optional_backends,
    rownames(installed.packages(priority = "base"))))
install.packages(deps, dependencies = NA)
# Install the exact advertised test floor after dependencies, without allowing
# any later dependency installation to replace it.
archive <- "https://cran.r-project.org/src/contrib/Archive/testthat/testthat_3.1.8.tar.gz"
install.packages(archive, repos = NULL, type = "source")
stopifnot(all(vapply(deps, requireNamespace, logical(1), quietly = TRUE)),
    identical(as.character(packageVersion("testthat")), "3.1.8"),
    !requireNamespace("ondisc", quietly = TRUE))
writeLines(c(paste("repository", repo), capture.output(sessionInfo())),
    "r-floor-dependencies.txt")
write.table(installed.packages()[, c("Package", "Version", "Built")],
    "r-floor-installed.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
