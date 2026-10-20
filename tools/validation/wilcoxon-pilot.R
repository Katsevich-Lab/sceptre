# Reproducible small scientific pilot: independent simulated screens, not null
# p-values drawn from Uniform. Example-source hash and full settings are saved.
args <- commandArgs(trailingOnly = TRUE)
replicates <- if (length(args)) as.integer(args[1L]) else 200L
out <- if (length(args) >= 2L) args[2L] else tempfile("wilcoxon-pilot-")
stopifnot(length(replicates) == 1L, !is.na(replicates), replicates >= 10L)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
library(sceptre)
source("tools/examples/wilcoxon.R")
wilson <- function(x, n) {
    z <- stats::qnorm(0.975)
    center <- (x / n + z^2 / (2 * n)) / (1 + z^2 / n)
    half <- z * sqrt((x / n) * (1 - x / n) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    c(lower = center - half, upper = center + half)
}
rows <- list()
start <- proc.time()
for (scenario in c("null", "effect", "confounded")) {
    for (i in seq_len(replicates)) {
        object <- suppressMessages(make_wilcoxon_example(100000L + i, scenario = scenario))
        result <- get_result(run_wilcoxon(object), "run_discovery_analysis")
        rows[[length(rows) + 1L]] <- data.frame(scenario = scenario,
            replicate = i, response_id = result$response_id,
            p_value = result$p_value, significant = result$significant)
    }
}
results <- do.call(rbind, rows)
summary <- do.call(rbind, lapply(split(results, paste(results$scenario, results$response_id)), function(x) {
    k <- sum(x$p_value < 0.05)
    ci <- wilson(k, nrow(x))
    data.frame(scenario = x$scenario[1L], response_id = x$response_id[1L],
        screens = nrow(x), rejected_raw_005 = k, rate = k / nrow(x),
        lower_95 = ci[1L], upper_95 = ci[2L])
}))
utils::write.csv(results, file.path(out, "pvalues.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
saveRDS(list(seed = "100000 + replicate (paired across scenarios)", n_cells = 400L,
    library_umis = 1000L, replicates = replicates, alpha = 0.05,
    method = "wilcoxon_unadjusted", normalized = "log1p_CP10K_full_library",
    example_sha256 = digest::digest(file = "tools/examples/wilcoxon.R", algo = "sha256"),
    harness_sha256 = digest::digest(file = "tools/validation/wilcoxon-pilot.R", algo = "sha256"),
    elapsed = proc.time() - start, session = sessionInfo()), file.path(out, "settings.rds"))
print(summary)
cat("Confounded scenario has no causal treatment effect on the null gene, but
",
    "treatment and batch differ: marginal distributions violate the rank null.
",
    "Its rejection rate illustrates the unadjusted method's causal limitation.
",
    "Intervals describe independent screens, not genes. Raw alpha .05 is not
",
    "the family's FDR; two-pair BH calls are retained separately in pvalues.csv.
",
    "This pilot is not real-data, donor-level, small-sample or general FDR validation.
", sep = "")
