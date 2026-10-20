# Read the actual W02 diagnostic; no fitting, draws or synthetic p-values.
args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1L] else "tools/validation/evidence/degenerate-design-review-20261019/w02-baseline-receipt.rds"
x <- readRDS(path)
cat("Input SHA256:", digest::digest(path, file = TRUE, algo = "sha256"), "\n")
cat("DGP:", x$distribution, "\nDesign:", x$design, "\n")
Z <- stats::model.matrix(~ treatment, x$covariates)
X <- x$covariates$treatment
cat("Rank Z:", qr(Z)$rank, "Rank [Z,X]:", qr(cbind(Z, X))$rank,
    "Max absolute X projection residual:", max(abs(qr.resid(qr(Z), X))), "\n")
for (mechanism in names(x$results)) {
    result <- x$results[[mechanism]]
    cat("\nMechanism:", mechanism, "Status:", result$status, "\n")
    print(result$raw_observed)
    print(result$warnings)
    if (!is.null(result$table)) {
        tab <- as.data.frame(result$table)
        print(tab[, intersect(c("response_id", "pass_qc", "p_value", "stage", "z_orig", "significant"),
                              names(tab)), drop = FALSE])
    }
    if (!is.null(result$message)) cat(result$message, "\n")
}
