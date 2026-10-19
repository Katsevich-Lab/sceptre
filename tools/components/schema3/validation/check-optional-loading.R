# Fresh R process: neither memory migration nor closed descriptor access loads ondisc.
pkgload::load_all(".", quiet = TRUE, recompile = FALSE)
stopifnot(!"ondisc" %in% loadedNamespaces())
legacy <- readRDS("tests/testthat/fixtures/saved-objects/sceptre-0.99.0-21f9ea0.rds")
current <- sceptre:::migrate_sceptre_object(legacy)
stopifnot(methods::is(sceptre::get_response_matrix(current), "dgRMatrix"),
    !"ondisc" %in% loadedNamespaces())
closed <- methods::new("BackedCountStore", descriptor = list(path = "/unused/closed.h5",
    content_md5 = paste(rep("0", 32), collapse = ""), dimensions = c(1L, 2L),
    dimnames = list("r1", c("c1", "c2")), integer_id = 1L), handle = NULL)
stopifnot(identical(dim(closed), c(1L, 2L)), identical(rownames(closed), "r1"),
    identical(colnames(closed), c("c1", "c2")), !"ondisc" %in% loadedNamespaces())
cat("PASS: package loading, memory migration/access, and closed descriptor inspection do not load ondisc.\n")
