# Structural fixtures only: these fixed toy p-values are NOT calibration data.
source("tools/validation/scientific-harness.R")
root <- tempfile("batch-integrity-")
dir.create(root)
make_fixture <- function(name, replicate, binary = "same") {
    path <- file.path(root, name)
    dir.create(path)
    rows <- data.frame(case = "toy", replicate = replicate,
        response_id = c("null", "down", "up"), true_log_effect = c(0, -1, 1),
        error = "", pass_qc = TRUE, p_value = c(.5, .01, .01))
    utils::write.csv(rows, file.path(path, "replicates.csv"), row.names = FALSE)
    saveRDS(list(finished_utc = "fixture", replicates = 1L,
        cases = data.frame(id = "toy"), installed_sha256 = binary,
        replicates_sha256 = digest::digest(file.path(path, "replicates.csv"), file = TRUE, algo = "sha256")),
        file.path(path, "manifest.rds"))
    path
}
a <- make_fixture("a", 1L)
b <- make_fixture("b", 2L)
c <- make_fixture("c", 3L, "different")
run <- function(name, inputs) {
    log <- file.path(root, paste0(name, ".log"))
    system2(file.path(R.home("bin"), "Rscript"),
        c("tools/validation/summarize-batches.R", shQuote(file.path(root, name)),
          vapply(inputs, shQuote, character(1))), stdout = log, stderr = log)
}
stopifnot(run("combined", c(a, b)) == 0L,
          nrow(utils::read.csv(file.path(root, "combined", "replicates.csv"))) == 6L,
          run("duplicate", c(a, a)) != 0L,
          run("wrong-binary", c(a, c)) != 0L)
m <- readRDS(file.path(b, "manifest.rds"))
m$finished_utc <- NULL
saveRDS(m, file.path(b, "manifest.rds"))
stopifnot(run("incomplete", c(a, b)) != 0L)
cat("Batch integrity checks passed: distinct indices accepted; duplicates, binary mismatch and nonterminal inputs rejected.\n")
unlink(root, recursive = TRUE)
