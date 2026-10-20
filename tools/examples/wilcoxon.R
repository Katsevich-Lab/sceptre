# Run from the repository root after installing the current package, or use
# pkgload::load_all("."). No external data or ordinary library edits are needed.
# A fixed multinomial library makes the null response exchangeable between arms.
# Increasing the effect gene draws mass only from the reservoir response.
make_wilcoxon_example <- function(seed = 821L, n = 400L,
                                 scenario = "effect") {
    set.seed(seed)
    assigned <- rep(1:4, length.out = n)
    treated <- assigned <= 2L
    batch <- if (scenario == "confounded") {
        stats::rbinom(n, 1, ifelse(treated, 0.9, 0.1))
    } else stats::rbinom(n, 1, 0.5)
    response <- vapply(seq_len(n), function(i) {
        effect_p <- if (scenario == "effect" && treated[i]) 0.04 else 0.02
        null_p <- if (scenario == "confounded" && batch[i] == 1L) 0.04 else 0.01
        stats::rmultinom(1L, 1000L,
            c(null_p, effect_p, 0.02, 1 - null_p - effect_p - 0.02))[, 1L]
    }, integer(4))
    dimnames(response) <- list(c("null", "effect", "other", "reservoir"), paste0("cell", seq_len(n)))
    guide <- matrix(0, 4L, n, dimnames = list(c("a1", "a2", "nt1", "nt2"), colnames(response)))
    guide[cbind(assigned, seq_len(n))] <- 50
    object <- sceptre::import_data(response_matrix = response, grna_matrix = guide,
        grna_target_data_frame = data.frame(grna_id = rownames(guide),
            grna_target = c("A", "A", "non-targeting", "non-targeting")),
        extra_covariates = data.frame(batch = batch, row.names = colnames(response)), moi = "low")
    pairs <- data.frame(grna_target = "A", response_id = c("null", "effect"))
    object <- sceptre::set_analysis_parameters(object, discovery_pairs = pairs,
        positive_control_pairs = pairs[2L, , drop = FALSE],
        formula_object = ~ batch, side = "both", control_group = "nt_cells")
    object <- sceptre::assign_grnas(object, method = "thresholding", threshold = 40,
        print_progress = FALSE)
    sceptre::run_qc(object, response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1), n_nonzero_trt_thresh = 7L,
        n_nonzero_cntrl_thresh = 7L)
}

# Keep the original for covariate-adjusted SCEPTRE comparisons.
if (sys.nframe() == 0L) {
    prepared <- suppressMessages(make_wilcoxon_example())
    alternative <- sceptre::run_wilcoxon(prepared)
    print(sceptre::get_result(alternative, "run_discovery_analysis"))
    # Known-null check and power check use this same explicitly selected method.
    alternative <- suppressMessages(sceptre::run_wilcoxon(alternative,
        "calibration_check", n_calibration_pairs = 4L, calibration_group_size = 1L))
    alternative <- sceptre::run_wilcoxon(alternative, "power_check")
    print(alternative)
    # For outputs/plots: write_outputs_to_directory(alternative, new_directory).
    # This example is not a calibration study; use the pilot harness below.
}
