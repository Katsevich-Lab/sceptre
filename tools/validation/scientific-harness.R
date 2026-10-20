# Independent datasets, not reference-output equality, are the replication unit.
# This file is sourceable without starting a simulation or loading sceptre.
validation_cases <- function() {
    data.frame(
        id = c("low_nt_empirical", "low_complement_empirical", "high_crt_empirical",
               "high_crt_skew", "high_permutation_skew", "confounded_crt",
               "omitted_confounder", "sparse_qc", "all_tied"),
        moi = c("low", "low", rep("high", 5), "low", "low"),
        control = c("nt_cells", rep("complement", 6), "nt_cells", "nt_cells"),
        mechanism = c("permutations", "permutations", "crt", "crt",
                      "permutations", "crt", "crt", "permutations", "permutations"),
        approximation = c(rep("no_approximation", 3), rep("skew_normal", 2),
                          rep("no_approximation", 4)),
        confounding = c(rep(FALSE, 5), TRUE, TRUE, FALSE, FALSE),
        omit_z = c(rep(FALSE, 6), TRUE, FALSE, FALSE),
        kind = c(rep("nb", 7), "sparse", "tied"),
        interpretation = c(rep("modelled", 6), "intentional_misspecification",
                           "eligibility_stress", "degenerate_stress"),
        stringsAsFactors = FALSE
    )
}

# Version2 leaves room for every admitted global replicate index within a case.
# Version1 used stride1000, overlapping adjacent case ranges beyond replicate500.
validation_seed <- function(case_id, replicate) {
    index <- match(case_id, validation_cases()$id)
    stopifnot(length(index) == 1L, !is.na(index), length(replicate) > 0L,
              all(is.finite(replicate)), all(replicate == floor(replicate)),
              all(replicate >= 1 & replicate <= 100000L))
    20261000L + 1000000L * index + 2L * as.integer(replicate)
}

validation_seed_rule <- function() {
    "v2: 20261000 + 1000000 * index in validation_cases() + replicate * 2; testing seed = DGP seed + 1; max global replicate100000"
}

validation_input <- function(case, seed, n = 600L) {
    withr::with_seed(seed, {
        z <- stats::rnorm(n)
        batch <- stats::rbinom(n, 1L, .5)
        probability <- stats::plogis(-1.4 + if (case$confounding) 1.3 * z + .6 * batch else 0)
        x <- stats::rbinom(n, 1L, probability)
        guide_ids <- c("gA", "gB", paste0("nt", 1:4))
        guides <- matrix(0, 6L, n, dimnames = list(guide_ids, paste0("c", seq_len(n))))
        guides[1L, ] <- 50 * x
        if (case$moi == "low") {
            other <- stats::rbinom(n, 1L, .3)
            group <- ifelse(x == 1L, 1L, ifelse(other == 1L, 2L, sample(3:6, n, TRUE)))
            guides[cbind(group, seq_len(n))] <- 50
        } else {
            # Independent additional guides create overlapping high-MOI groups.
            guides[2:6, ] <- 50 * matrix(stats::rbinom(5L * n, 1L, .3), 5L, n)
        }
        effects <- c(null = 0, down = log(.5), up = log(2))
        baseline <- exp(1 + .8 * z + .4 * batch)
        response <- t(vapply(effects, function(effect) {
            stats::rnbinom(n, size = 3, mu = baseline * exp(effect * x))
        }, numeric(n)))
        if (case$kind == "sparse") response[,] <- stats::rpois(length(response), .02)
        if (case$kind == "tied") response[,] <- 2L
        # Background outcomes keep complete-library QC/depth well-defined;
        # they are not tested and do not enter the explicitly specified design.
        background <- matrix(stats::rpois(8L * n, 5), 8L, n,
                             dimnames = list(paste0("background", 1:8), NULL))
        response <- rbind(response, background)
        colnames(response) <- colnames(guides)
        list(response = response, guides = guides,
             metadata = data.frame(grna_id = guide_ids,
                 grna_target = c("A", "B", rep("non-targeting", 4))),
             covariates = data.frame(z = z, batch = factor(batch), row.names = colnames(guides)),
             pairs = data.frame(grna_target = "A", response_id = names(effects)),
             truth = if (case$kind == "nb") effects else stats::setNames(rep(0, 3), names(effects)),
             assignment = x, assignment_probability = probability)
    })
}

validation_prepare <- function(input, case) {
    object <- sceptre::import_data(response_matrix = input$response,
        grna_matrix = input$guides, grna_target_data_frame = input$metadata,
        extra_covariates = input$covariates, moi = case$moi, use_ondisc = FALSE)
    object <- sceptre::set_analysis_parameters(object, discovery_pairs = input$pairs,
        formula_object = if (case$kind == "tied") ~ 1 else if (case$omit_z) ~ batch else ~ z + batch,
        side = "both", control_group = case$control,
        resampling_mechanism = case$mechanism,
        resampling_approximation = case$approximation)
    object <- sceptre::assign_grnas(object, method = "thresholding", threshold = 5L,
                                    print_progress = FALSE)
    object <- sceptre::run_qc(object, response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1), n_nonzero_trt_thresh = 7L,
        n_nonzero_cntrl_thresh = 7L)
    # Configure only the fully prepared object, as candidate fit APIs require.
    configure <- getOption("sceptre.validation.configure")
    if (!is.null(configure)) object <- configure(object, case)
    object
}

validation_replicate <- function(case, replicate, seed, n = 600L) {
    input <- validation_input(case, seed, n)
    warnings <- character()
    started <- proc.time()[[3L]]
    value <- tryCatch(withCallingHandlers(suppressMessages({
        object <- validation_prepare(input, case)
        # Preserve the public resampling policy, including adaptive stage exits.
        runner <- getOption("sceptre.validation.run")
        withr::with_seed(seed + 1L, {
            if (is.null(runner)) {
                fitted <- sceptre::run_discovery_analysis(object,
                    print_progress = FALSE, output_amount = 2L, parallel = FALSE)
                sceptre::get_result(fitted, "run_discovery_analysis")
            } else runner(object, case)
        })
    }), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
    elapsed <- proc.time()[[3L]] - started
    rows <- data.frame(case = case$id, replicate = replicate, seed = seed,
        response_id = names(input$truth), true_log_effect = unname(input$truth),
        input_sha256 = digest::digest(input, algo = "sha256"), n_cells = n,
        n_treated = sum(input$assignment), elapsed_seconds = elapsed,
        warning = paste(unique(warnings), collapse = " | "), error = "",
        pass_qc = FALSE, p_value = NA_real_, stage = NA_integer_)
    if (inherits(value, "error")) {
        if (identical(conditionMessage(value), "Zero pairs pass pairwise QC. Cannot run analysis.")) {
            rows$warning <- paste(rows$warning, "All requested pairs failed QC", sep = " | ")
        } else rows$error <- conditionMessage(value)
    } else {
        idx <- match(rows$response_id, value$response_id)
        rows$pass_qc <- value$pass_qc[idx]
        rows$p_value <- value$p_value[idx]
        if ("stage" %in% names(value)) rows$stage <- value$stage[idx]
    }
    rows
}

validation_interval <- function(k, n, level = .95) {
    if (!n) return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
    ci <- stats::binom.test(k, n, conf.level = level)$conf.int
    c(estimate = k / n, lower = ci[1L], upper = ci[2L])
}

validation_summary <- function(rows, alphas = c(.01, .05, .1)) {
    groups <- split(rows, interaction(rows$case, rows$response_id, drop = TRUE))
    do.call(rbind, lapply(groups, function(group) {
        valid <- group$error == "" & !is.na(group$pass_qc) & group$pass_qc &
            is.finite(group$p_value) & group$p_value >= 0 & group$p_value <= 1
        do.call(rbind, lapply(alphas, function(alpha) {
            k <- sum(group$p_value[valid] <= alpha)
            ci <- validation_interval(k, sum(valid))
            data.frame(case = group$case[1L], response_id = group$response_id[1L],
                true_log_effect = group$true_log_effect[1L], alpha = alpha,
                attempted = nrow(group), valid = sum(valid), errors = sum(group$error != ""),
                qc_failed = sum(!group$pass_qc & group$error == "", na.rm = TRUE),
                invalid_p = sum(group$pass_qc & group$error == "" &
                    (!is.finite(group$p_value) | group$p_value < 0 | group$p_value > 1), na.rm = TRUE),
                rejected = k, estimate = ci[1L], lower_95 = ci[2L], upper_95 = ci[3L],
                mc_se = if (sum(valid)) sqrt(ci[1L] * (1 - ci[1L]) / sum(valid)) else NA_real_,
                # The unconditional denominator keeps QC/errors visible.
                rejection_per_attempt = k / nrow(group), row.names = NULL)
        }))
    }))
}
