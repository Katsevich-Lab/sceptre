#!/usr/bin/env Rscript
# Actual null outcomes with treatment exactly in nuisance design. This is an
# unsupported-identifiability boundary diagnostic, not a modelled calibration DGP.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: positivity-boundary.R NEW_OUTPUT_DIR [DATASETS=40]")
out <- normalizePath(args[1L], mustWork = FALSE)
if (dir.exists(out)) stop("Choose a new output directory")
reps <- if (length(args) >= 2L) as.integer(args[2L]) else 40L
stopifnot(length(reps) == 1L, !is.na(reps), reps > 0L, reps <= 1000L)
library(sceptre)
data.table::setDTthreads(1L)
dir.create(out, recursive = TRUE)
installed <- find.package("sceptre")
files <- list.files(installed, recursive = TRUE, full.names = TRUE)
manifest <- list(source_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    source_status = system2("git", c("status", "--porcelain"), stdout = TRUE),
    installed = installed, installed_sha256 = stats::setNames(vapply(files, digest::digest,
        character(1), file = TRUE, algo = "sha256"), substring(files, nchar(installed) + 2L)),
    script_sha256 = digest::digest("tools/validation/positivity-boundary.R", file = TRUE, algo = "sha256"),
    seed_rule = "23703 + dataset - 1; testing seed = seed + 100000",
    datasets = reps, n_cells = 160L, response_genes = 4L,
    DGP = "Independent Poisson(5); all genes null; fixed guides with exact treatment indicator in nuisance design",
    request = "Each gene tested in its OWN public one-pair discovery request, to expose finite accepted scores otherwise hidden by another gene error",
    alpha = .05, interval_scope = "No pooled binomial precision: genes and mechanisms share each dataset",
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), session = utils::sessionInfo())
saveRDS(manifest, file.path(out, "manifest.rds"))
rows <- list()
for (i in seq_len(reps)) {
    seed <- 23703L + i - 1L
    n <- 160L
    cells <- paste0("cell", seq_len(n))
    assigned <- rep(1:4, length.out = n)
    X <- as.integer(assigned <= 2L)
    guide <- matrix(0, 4L, n, dimnames = list(c("a1", "a2", "nt1", "nt2"), cells))
    guide[cbind(assigned, seq_len(n))] <- 50
    response <- withr::with_seed(seed, matrix(stats::rpois(4L * n, 5), 4L,
        dimnames = list(paste0("r", 1:4), cells)))
    covariates <- data.frame(treatment = X, row.names = cells)
    metadata <- data.frame(grna_id = rownames(guide),
        grna_target = c("A", "A", "non-targeting", "non-targeting"))
    Z <- stats::model.matrix(~ treatment, covariates)
    residual <- max(abs(qr.resid(qr(Z), X)))
    input_hash <- digest::digest(list(response, guide, covariates), algo = "sha256")
    for (mechanism in c("permutations", "crt")) for (gene in rownames(response)) {
        warnings <- character()
        value <- tryCatch(withCallingHandlers(suppressMessages({
            object <- sceptre::import_data(response, guide, metadata, moi = "low", extra_covariates = covariates)
            object <- sceptre::set_analysis_parameters(object,
                discovery_pairs = data.frame(grna_target = "A", response_id = gene),
                formula_object = ~ treatment, control_group = "complement", side = "both",
                resampling_mechanism = mechanism, resampling_approximation = "no_approximation")
            object <- sceptre::assign_grnas(object, method = "thresholding", threshold = 40L, print_progress = FALSE)
            object <- sceptre::run_qc(object, response_n_umis_range = c(0, 1),
                response_n_nonzero_range = c(0, 1), n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L)
            fitted <- withr::with_seed(seed + 100000L, sceptre::run_discovery_analysis(object,
                print_progress = FALSE, output_amount = 2L))
            sceptre::get_result(fitted, "run_discovery_analysis")
        }), warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
        }), error = function(e) e)
        row <- data.frame(dataset = i, seed = seed, input_sha256 = input_hash,
            mechanism = mechanism, response_id = gene, treatment_projection_max_residual = residual,
            error = "", warning = paste(unique(warnings), collapse = " | "),
            pass_qc = NA, observed_score = NA_real_, p_value = NA_real_, significant = NA)
        if (inherits(value, "error")) row$error <- conditionMessage(value) else {
            row$pass_qc <- value$pass_qc[1L]
            row$observed_score <- value$z_orig[1L]
            row$p_value <- value$p_value[1L]
            row$significant <- value$significant[1L]
        }
        rows[[length(rows) + 1L]] <- row
        utils::write.table(row, file.path(out, "requests.csv"), sep = ",", row.names = FALSE,
            col.names = length(rows) == 1L, append = length(rows) != 1L, na = "NA")
    }
    if (i %% 10L == 0L || i == reps) cat("Boundary dataset", i, "/", reps, "\n")
}
all_rows <- do.call(rbind, rows)
summary <- do.call(rbind, lapply(split(all_rows, all_rows$mechanism), function(x) {
    completed <- x$error == ""
    finite <- completed & is.finite(x$observed_score)
    data.frame(mechanism = x$mechanism[1L], datasets = reps, one_pair_requests = nrow(x),
        errors = sum(!completed), completed = sum(completed), finite_completed = sum(finite),
        nonfinite_completed = sum(completed & !is.finite(x$observed_score)),
        finite_p_below_005 = sum(finite & is.finite(x$p_value) & x$p_value < .05),
        min_finite_completed_p = if (any(finite)) min(x$p_value[finite], na.rm = TRUE) else NA_real_)
}))
utils::write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
manifest$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
manifest$requests_sha256 <- digest::digest(file.path(out, "requests.csv"), file = TRUE, algo = "sha256")
saveRDS(manifest, file.path(out, "manifest.rds"))
print(summary)
