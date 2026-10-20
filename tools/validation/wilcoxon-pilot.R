# Run from the repository root against the exact installed candidate.
# Independent simulated screens, not synthetic Uniform null p-values.
wilson <- function(x, n) {
    z <- stats::qnorm(0.975)
    center <- (x / n + z^2 / (2 * n)) / (1 + z^2 / n)
    half <- z * sqrt((x / n) * (1 - x / n) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    c(lower = center - half, upper = center + half)
}

validate_pilot_result <- function(result) {
    fields <- c("response_id", "pass_qc", "p_value", "significant", "inference_method")
    if (!all(fields %in% names(result)) || nrow(result) != 2L ||
        !setequal(result$response_id, c("null", "effect")) ||
        !is.logical(result$pass_qc) || anyNA(result$pass_qc) || !all(result$pass_qc) ||
        !is.numeric(result$p_value) || anyNA(result$p_value) || any(!is.finite(result$p_value)) ||
        any(result$p_value < 0 | result$p_value > 1) ||
        !is.logical(result$significant) || anyNA(result$significant) ||
        anyNA(result$inference_method) ||
        !all(result$inference_method == "wilcoxon_unadjusted")) {
        stop("Unexpected pilot QC, family, method or probability result.")
    }
    invisible(NULL)
}

pilot_main <- function(args = commandArgs(trailingOnly = TRUE)) {
    replicates <- if (length(args)) as.integer(args[1L]) else 200L
    out <- if (length(args) >= 2L) args[2L] else tempfile("wilcoxon-pilot-")
    stopifnot(length(replicates) == 1L, !is.na(replicates), replicates >= 10L)
    if (file.exists(out)) stop("Pilot output already exists; choose a fresh directory.")
    if (!dir.create(out, recursive = TRUE)) stop("Cannot create pilot output.")
    out <- normalizePath(out)
    utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    started <- utc()
    state <- list(status = "running", started_utc = started, completed_screens = 0L)
    checkpoint <- function(status, condition = NULL) {
        state$status <- status
        state$updated_utc <- utc()
        state$condition <- condition
        temporary <- tempfile("status-", tmpdir = out)
        saveRDS(state, temporary)
        if (!file.rename(temporary, file.path(out, "status.rds"))) {
            stop("Cannot commit pilot checkpoint.")
        }
    }
    checkpoint("running")
    tryCatch({
        library(sceptre)
        source("tools/examples/wilcoxon.R")
        git <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
        dirty <- system2("git", c("status", "--porcelain"), stdout = TRUE)
        stopifnot(length(git) == 1L, grepl("^[a-f0-9]{40}$", git), !length(dirty))
        archive <- Sys.getenv("SCEPTRE_PILOT_ARCHIVE", "")
        settings <- list(
            started_utc = started, git_commit = git, git_status = dirty,
            source_directory = normalizePath("."),
            installed_package = find.package("sceptre"),
            runtime = sceptre:::execution_build_certificate(),
            archive_sha256 = if (nzchar(archive)) {
                digest::digest(file = archive, algo = "sha256")
            } else NULL,
            seed = "100000 + replicate; paired across scenarios",
            n_cells = 400L, library_umis = 1000L, replicates = replicates,
            alpha_raw = 0.05, alpha_BH_family = 0.1,
            method = "wilcoxon_unadjusted",
            normalized = "log1p_CP10K_full_library",
            example_sha256 = digest::digest(file = "tools/examples/wilcoxon.R", algo = "sha256"),
            harness_sha256 = digest::digest(file = "tools/validation/wilcoxon-pilot.R", algo = "sha256"),
            session = sessionInfo()
        )
        saveRDS(settings, file.path(out, "settings.rds"))
        rows <- list()
        start <- proc.time()
        for (scenario in c("null", "effect", "confounded")) {
            for (i in seq_len(replicates)) {
                state$scenario <- scenario
                state$replicate <- i
                object <- suppressMessages(make_wilcoxon_example(100000L + i, scenario = scenario))
                result <- get_result(run_wilcoxon(object), "run_discovery_analysis")
                validate_pilot_result(result)
                rows[[length(rows) + 1L]] <- data.frame(scenario = scenario,
                    replicate = i, response_id = result$response_id,
                    p_value = result$p_value, significant = result$significant)
                state$completed_screens <- state$completed_screens + 1L
                if (i %% 25L == 0L) checkpoint("running")
            }
            saveRDS(rows, file.path(out, "checkpoint-rows.rds"))
            checkpoint("running")
        }
        results <- do.call(rbind, rows)
        stopifnot(nrow(results) == 6L * replicates)
        summary <- do.call(rbind, lapply(split(results, paste(results$scenario, results$response_id)), function(x) {
            stopifnot(nrow(x) == replicates, !anyDuplicated(x$replicate))
            k <- sum(x$p_value < 0.05)
            ci <- wilson(k, nrow(x))
            data.frame(scenario = x$scenario[1L], response_id = x$response_id[1L],
                screens = nrow(x), rejected_raw_005 = k, rate = k / nrow(x),
                lower_95 = ci[1L], upper_95 = ci[2L])
        }))
        utils::write.csv(results, file.path(out, "pvalues.csv"), row.names = FALSE)
        utils::write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
        state$elapsed <- proc.time() - start
        state$artifact_sha256 <- lapply(c("pvalues.csv", "summary.csv", "settings.rds"), function(name) {
            digest::digest(file = file.path(out, name), algo = "sha256")
        })
        names(state$artifact_sha256) <- c("pvalues.csv", "summary.csv", "settings.rds")
        checkpoint("complete")
        print(summary)
        message("TERMINAL_SUCCESS: ", state$completed_screens, " independent screens; ", out)
        message("Confounded conditional null violates the marginal rank null; it is not a type-I error measurement.")
        message("Intervals use independent screens per endpoint. Raw alpha .05 rates are not family FDR.")
        message("This pilot does not establish real-data, donor-level, small-sample or general FDR validity.")
    }, error = function(error) {
        checkpoint("failed", conditionMessage(error))
        stop(error)
    }, interrupt = function(interrupt) {
        checkpoint("interrupted", conditionMessage(interrupt))
        stop(interrupt)
    })
}

if (sys.nframe() == 0L) pilot_main()
