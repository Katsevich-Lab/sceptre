# Run only after legitimate host capacity admission. No dependency installs.
# Rebuild actual pkgdown contributor sections in the existing docs destination.
# Existing user vignette pages are retained without rerunning their analyses.
# Reference help is rendered without running examples or loading source code.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || args[[1L]] != "--capacity-admitted") {
    stop("Host capacity admission is required; this flag is not itself a grant.")
}
root <- normalizePath(getwd(), mustWork = TRUE)
stopifnot(file.exists("_pkgdown.yml"), file.exists("tools/architecture/pkgdown-pages.json"))
pkg <- pkgdown::as_pkgdown(root)
pkgdown::init_site(pkg)
pkgdown::build_home(pkg, preview = FALSE, quiet = FALSE)
pkgdown::build_reference(pkg, lazy = FALSE, examples = FALSE, devel = FALSE, preview = FALSE)
pages <- jsonlite::fromJSON("tools/architecture/pkgdown-pages.json")$pages
for (article in pages$article) {
    name <- sub("[.]Rmd$", "", sub("^vignettes/", "", article))
    pkgdown::build_article(name, pkg, lazy = FALSE, new_process = FALSE, quiet = FALSE)
}
pkgdown::build_articles_index(pkg)
pkgdown::build_search(pkg)
cat("Actual pkgdown contributor sections rebuilt.\n")
cat("Reference help refreshed without examples; existing user-vignette output retained.\n")
cat("No full-site or package-check claim.\n")
