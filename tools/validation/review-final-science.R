#!/usr/bin/env Rscript
# Independent OFFLINE readback of terminal final-source science. No package
# inference or simulated-data generation; do not duplicate a live study.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(8L, 9L)) stop("Usage: review-final-science.R SCIENCE_ROOT FINAL_SOURCE_ROOT ARCHIVE PACKAGE_LIB BACKEND_LIB NEW_OUTPUT EXPECTED_COMMIT EXPECTED_ARCHIVE_SHA [allow-shared-inputs]")
science <- normalizePath(args[1L]); source_root <- normalizePath(args[2L])
archive <- normalizePath(args[3L]); lib <- normalizePath(args[4L]); backend <- normalizePath(args[5L])
out <- normalizePath(args[6L], mustWork = FALSE)
if (dir.exists(out)) stop("Review output already exists")
revision <- args[7L]
expected_archive <- args[8L]
allow_shared <- length(args) == 9L && identical(args[9L], "allow-shared-inputs")
stopifnot(grepl("^[a-f0-9]{40}$", revision), grepl("^[a-f0-9]{64}$", expected_archive))
hash <- function(path) digest::digest(path, file = TRUE, algo = "sha256")
stopifnot(file.exists(file.path(science, "COMPLETE")), identical(hash(archive), expected_archive))
head <- system2("git", c("-C", shQuote(source_root), "rev-parse", "HEAD"), stdout = TRUE)
stopifnot(identical(head, revision), !length(system2("git", c("-C", shQuote(source_root), "status", "--porcelain"), stdout = TRUE)))
loaded_dll_hash <- hash(file.path(lib, "sceptre", "libs", "sceptre.so"))
case_ids <- c("low_nt_empirical", "high_crt_empirical", "high_crt_skew", "confounded_crt")
case_indices <- c(low_nt_empirical = 1L, high_crt_empirical = 3L, high_crt_skew = 4L, confounded_crt = 6L)
batch_dirs <- file.path(science, paste0("batch-", seq.int(1L, 901L, 100L)))
read_rows <- function(path) utils::read.csv(path, stringsAsFactors = FALSE, colClasses = c(error = "character", warning = "character"))
check_manifest <- function(path) {
    m <- readRDS(file.path(path, "manifest.rds"))
    stopifnot(!is.null(m$finished_utc), identical(m$source_revision, revision), !length(m$source_status))
    package_path <- if (!is.null(m$installed_package)) m$installed_package else m$installed
    stopifnot(identical(normalizePath(package_path), normalizePath(file.path(lib, "sceptre"))))
    stopifnot(identical(unname(m$installed_sha256[["libs/sceptre.so"]]), loaded_dll_hash))
    # Actual installed library file identities must still match measurement.
    actual <- vapply(file.path(lib, "sceptre", names(m$installed_sha256)), hash, character(1))
    stopifnot(identical(unname(actual), unname(m$installed_sha256)))
    if (!is.null(m$source_sha256)) {
        actual_source <- vapply(file.path(source_root, names(m$source_sha256)), hash, character(1))
        stopifnot(identical(unname(actual_source), unname(m$source_sha256)))
    }
    if (!is.null(m$harness_sha256)) {
        actual_harness <- vapply(file.path(source_root, names(m$harness_sha256)), hash, character(1))
        stopifnot(identical(unname(actual_harness), unname(m$harness_sha256)))
    }
    m
}
manifest_rows <- lapply(batch_dirs, check_manifest)
stopifnot(all(vapply(manifest_rows, function(m) is.null(m$setup), logical(1))))
for (i in seq_along(batch_dirs)) {
    m <- manifest_rows[[i]]
    stopifnot(m$replicates == 100L, m$first_replicate == 1L + 100L * (i - 1L),
              m$last_replicate == 100L * i, m$n_cells == 600L,
              identical(as.character(m$cases$id), case_ids),
              identical(hash(file.path(batch_dirs[i], "replicates.csv")), m$replicates_sha256))
}
rows <- read_rows(file.path(science, "combined", "replicates.csv"))
combined <- readRDS(file.path(science, "combined", "combined-manifest.rds"))
stopifnot(!is.null(combined$finished_utc), identical(hash(file.path(science, "combined", "replicates.csv")), combined$combined_csv_sha256))
batch_rows <- do.call(rbind, lapply(batch_dirs, function(p) read_rows(file.path(p, "replicates.csv"))))
stopifnot(identical(rows, batch_rows), nrow(rows) == 12000L, !anyDuplicated(rows[c("case", "replicate", "response_id")]),
          setequal(rows$case, case_ids), all(rows$error == ""), all(rows$pass_qc),
          all(is.finite(rows$p_value)), all(rows$p_value >= 0 & rows$p_value <= 1), all(rows$n_cells == 600L))
truth <- c(null = 0, down = log(.5), up = log(2))
seed_versions <- vapply(manifest_rows, function(m) if (is.null(m$seed_version)) "v1-legacy" else m$seed_version, character(1))
stopifnot(length(unique(seed_versions)) == 1L, seed_versions[1L] %in% c("v1-legacy", "v2-disjoint-case-ranges"))
stride <- if (seed_versions[1L] == "v2-disjoint-case-ranges") 1000000L else 1000L
stopifnot(all(abs(rows$true_log_effect - truth[rows$response_id]) < 1e-14),
          all(rows$seed == 20261000L + stride * case_indices[rows$case] + 2L * rows$replicate))
for (id in case_ids) for (endpoint in names(truth)) {
    x <- rows[rows$case == id & rows$response_id == endpoint, ]
    stopifnot(nrow(x) == 1000L, identical(sort(x$replicate), 1:1000))
}
# Each dataset has one input hash across its three tested outcomes; no pooling.
stopifnot(all(vapply(split(rows$input_sha256, paste(rows$case, rows$replicate)),
                     function(x) length(unique(x)) == 1L, logical(1))))
input_ids <- unique(rows[c("case", "replicate", "input_sha256")])
stopifnot(nrow(input_ids) == 4000L, !anyDuplicated(input_ids[c("case", "replicate")]))
unique_inputs <- length(unique(input_ids$input_sha256))
if (unique_inputs != 4000L && !allow_shared)
    stop("Study has ", unique_inputs, " unique inputs across4000case evaluations; final independent target unmet")
if (!is.null(combined$unique_input_datasets)) stopifnot(combined$unique_input_datasets == unique_inputs)
summary <- do.call(rbind, lapply(case_ids, function(id) do.call(rbind, lapply(names(truth), function(endpoint) {
    x <- rows[rows$case == id & rows$response_id == endpoint, ]
    do.call(rbind, lapply(c(.01, .05, .1), function(alpha) {
        k <- sum(x$p_value <= alpha); ci <- stats::binom.test(k, 1000L)$conf.int
        data.frame(case = id, response_id = endpoint, true_log_effect = truth[endpoint],
            alpha = alpha, attempted = 1000L, valid = 1000L, rejected = k, estimate = k / 1000L,
            lower_95 = ci[1L], upper_95 = ci[2L], mc_se = sqrt(k / 1000 * (1 - k / 1000) / 1000), row.names = NULL)
    }))
}))))
original <- utils::read.csv(file.path(science, "combined", "summary.csv"))
key <- function(x) paste(x$case, x$response_id, x$alpha, sep = "/")
idx <- match(key(summary), key(original))
stopifnot(!anyNA(idx), nrow(summary) == nrow(original))
for (field in c("attempted", "valid", "rejected", "estimate", "lower_95", "upper_95", "mc_se"))
    stopifnot(all(abs(summary[[field]] - original[[field]][idx]) < 1e-12))
stopifnot(all(original$errors == 0), all(original$invalid_p == 0), all(original$qc_failed == 0))
# Separate actual-count tied regression and unsupported-design boundary.
tied_manifest <- check_manifest(file.path(science, "all-tied"))
tied <- read_rows(file.path(science, "all-tied", "replicates.csv"))
stopifnot(identical(hash(file.path(science, "all-tied", "replicates.csv")), tied_manifest$replicates_sha256),
          nrow(tied) == 120L, all(tied$error == ""), all(tied$pass_qc), all(tied$p_value == 1))
boundary_manifest <- check_manifest(file.path(science, "positivity"))
boundary <- read_rows(file.path(science, "positivity", "requests.csv"))
stopifnot(identical(hash(file.path(science, "positivity", "requests.csv")), boundary_manifest$requests_sha256), nrow(boundary) == 320L)
completed <- boundary$error == ""
stopifnot(all(grepl("Observed association score is not finite", boundary$error[!completed], fixed = TRUE)),
          all(boundary$pass_qc[completed]), all(is.finite(boundary$observed_score[completed])),
          all(is.finite(boundary$p_value[completed])), all(boundary$p_value[completed] >= .05 & boundary$p_value[completed] <= 1))
backend_desc <- read.dcf(file.path(backend, "ondisc", "DESCRIPTION"))
stopifnot(identical(unname(backend_desc[1L, "RemoteSha"]), "589803edea2a66ad15e7f07186bee16c7739ab8c"))
dir.create(out, recursive = TRUE)
utils::write.csv(summary, file.path(out, "independent-summary.csv"), row.names = FALSE)
review <- list(status = if (unique_inputs == 4000L) "independently_verified_terminal_unique4000" else "verified_case_specific_historical_SHARED_INPUTS_TARGET_UNMET", reviewed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    source_revision = revision, archive_sha256 = expected_archive, loaded_sceptre_DLL_sha256 = loaded_dll_hash,
    case_evaluations = 4000L, unique_input_datasets = unique_inputs, per_case = 1000L, endpoint_rows = 12000L,
    final_unique_input_target_met = unique_inputs == 4000L, seed_version = seed_versions[1L],
    across_case_shared_input_excess = 4000L - unique_inputs,
    interval = "Pointwise exact Clopper-Pearson95%; datasets per endpoint/case; no pooled pilot duplicates",
    batch_indices = "1:1000 exactly once; pilotfirst40 deliberately not added", no_setup = TRUE,
    tied = c(datasets = 40L, rows = nrow(tied), errors = sum(tied$error != ""), p1 = sum(tied$p_value == 1)),
    boundary = list(datasets = 40L, one_pair_requests = 320L, errors = sum(!completed), finite_completed = sum(completed),
        finite_significant = sum(boundary$p_value[completed] < .05),
        scope = "Unsupported-identifiability diagnostic; finite score not proof of estimability; genes/mechanisms not independent replicas"),
    backend_available = list(pin = backend_desc[1L, "RemoteSha"], version = backend_desc[1L, "Version"],
        DLL_sha256 = hash(file.path(backend, "ondisc", "libs", "ondisc.so")),
        loaded_in_scientific_snapshot = "ondisc" %in% names(manifest_rows[[1L]]$session$loadedOnly),
        interpretation = "In-memory scientific runs, no backed/ondisc statistical validity claim"),
    rows_sha256 = hash(file.path(science, "combined", "replicates.csv")),
    warnings_by_case = table(rows$case[nzchar(rows$warning)]),
    original_summary_recomputed = TRUE, actual_batch_rows_equal_combined = TRUE,
    manifest_hashes = vapply(c(file.path(batch_dirs, "manifest.rds"),file.path(science,"combined","combined-manifest.rds")), hash, character(1)))
saveRDS(review, file.path(out, "independent-review.rds"))
jsonlite::write_json(review, file.path(out, "independent-review.json"), auto_unbox = TRUE, pretty = TRUE)
print(summary[summary$alpha == .05, ])
cat("Independent offline readback complete; uniqueinputs:", unique_inputs, "/4000; no inference rerun.\n")
