# Shared deterministic inputs for package regression tests, adapted from
# tests/testthat/helper-refactor.R on the refactor branch (#234).

refactor_make_inputs <- function() {
    withr::with_seed(20261016, {
        n <- 160L
        guide_ids <- c("a1", "a2", "b1", "b2", "nt1", "nt2", "nt3", "nt4")
        targets <- c("A", "A", "B", "B", rep("non-targeting", 4L))
        guide <- matrix(
            0,
            length(guide_ids),
            n,
            dimnames = list(guide_ids, paste0("cell_", seq_len(n)))
        )
        assigned <- rep(seq_along(guide_ids), length.out = n)
        guide[cbind(assigned, seq_len(n))] <- 50
        guide[1L, which(assigned == 1L)[seq_len(4L)]] <- 43
        response <- matrix(
            stats::rpois(8L * n, 3),
            8L,
            n,
            dimnames = list(paste0("r", seq_len(8L)), colnames(guide))
        )
        response[1L, assigned %in% 1:2] <- response[1L, assigned %in% 1:2] + 20
        response[2L, assigned %in% 3:4] <- 0
        response[3L, ] <- 0
        response[3L, assigned == 1L] <- 5
        response[3L, assigned >= 5L] <- 3
        response[4L, ] <- 0
        high <- guide
        high[1L, seq(10L, n, 16L)] <- 50
        high[3L, seq(5L, n, 16L)] <- 50
        high[6L, seq(3L, n, 16L)] <- 50
        list(
            response_matrix = response,
            low_grna_matrix = guide,
            high_grna_matrix = high,
            grna_target_data_frame = data.frame(
                grna_id = guide_ids,
                grna_target = targets
            ),
            extra_covariates = data.frame(
                batch = factor(rep(c("b", "a", "b", "b", "a"), length.out = n)),
                z = sin(seq_len(n)),
                row.names = colnames(guide)
            ),
            discovery_pairs = data.frame(
                grna_target = c("B", "A", "A", "B", "A"),
                response_id = c("r5", "r3", "r1", "r2", "r4")
            ),
            positive_control_pairs = data.frame(
                grna_target = c("A", "B"),
                response_id = c("r1", "r2")
            ),
            remove_cells = c(2L, 17L, 49L, 103L, 158L)
        )
    })
}

refactor_case_definitions <- function() {
    list(
        low_nt_perm = list(
            moi = "low",
            control = "nt_cells",
            mechanism = "permutations",
            integration = "union",
            side = "both",
            approximation = "skew_normal",
            output = 3L
        ),
        low_complement_perm = list(
            moi = "low",
            control = "complement",
            mechanism = "permutations",
            integration = "union",
            side = "left",
            approximation = "no_approximation",
            output = 2L
        ),
        high_complement_crt = list(
            moi = "high",
            control = "complement",
            mechanism = "crt",
            integration = "union",
            side = "right",
            approximation = "skew_normal",
            output = 3L
        ),
        low_complement_crt = list(
            moi = "low",
            control = "complement",
            mechanism = "crt",
            integration = "union",
            side = "both",
            approximation = "no_approximation",
            output = 1L
        ),
        high_perm = list(
            moi = "high",
            control = "complement",
            mechanism = "permutations",
            integration = "union",
            side = "both",
            approximation = "no_approximation",
            output = 3L
        ),
        low_bonferroni = list(
            moi = "low",
            control = "nt_cells",
            mechanism = "permutations",
            integration = "bonferroni",
            side = "left",
            approximation = "skew_normal",
            output = 2L
        ),
        low_singleton = list(
            moi = "low",
            control = "complement",
            mechanism = "permutations",
            integration = "singleton",
            side = "right",
            approximation = "no_approximation",
            output = 1L
        )
    )
}

refactor_import <- function(inputs, case) {
    sceptre::import_data(
        response_matrix = inputs$response_matrix,
        grna_matrix = inputs[[paste0(case$moi, "_grna_matrix")]],
        grna_target_data_frame = inputs$grna_target_data_frame,
        extra_covariates = inputs$extra_covariates,
        moi = case$moi
    )
}

refactor_prepare <- function(inputs, case, assignment_threshold = 40) {
    object <- refactor_import(inputs, case)
    object <- sceptre::set_analysis_parameters(
        object,
        discovery_pairs = inputs$discovery_pairs,
        positive_control_pairs = inputs$positive_control_pairs,
        formula_object = stats::as.formula("~ batch + z", env = baseenv()),
        side = case$side,
        grna_integration_strategy = case$integration,
        control_group = case$control,
        resampling_mechanism = case$mechanism,
        resampling_approximation = case$approximation
    )
    object <- sceptre::assign_grnas(
        object,
        method = "thresholding",
        threshold = assignment_threshold,
        print_progress = FALSE
    )
    sceptre::run_qc(
        object,
        response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1),
        n_nonzero_trt_thresh = 7L,
        n_nonzero_cntrl_thresh = 7L,
        additional_cells_to_remove = inputs$remove_cells
    )
}
