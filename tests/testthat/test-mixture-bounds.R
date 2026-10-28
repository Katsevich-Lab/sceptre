mixture_parameters <- function() {
    list(
        n_em_rep = 1,
        pi_guess_range = c(.01, .1),
        g_pert_guess_range = c(0, 20),
        n_nonzero_cells_cutoff = 1,
        backup_threshold = 1,
        probability_threshold = .8,
        formula_object = ~1
    )
}

test_that("mixture perturbation bounds use the documented perturbation range", {
    object <- methods::new("sceptre_object", low_moi = TRUE)
    for (upper in c(19.9, 20)) {
        parameters <- mixture_parameters()
        parameters$g_pert_guess_range[2L] <- upper
        expect_null(check_assign_grna_inputs(
            object,
            "mixture",
            parameters,
            "auto"
        ))
    }
    parameters$g_pert_guess_range[2L] <- 20.1
    expect_error(
        check_assign_grna_inputs(object, "mixture", parameters, "auto"),
        "`g_pert_guess_range` should be a numeric vector of length 2",
        fixed = TRUE
    )
    parameters$pi_guess_range[2L] <- 2
    expect_error(
        check_assign_grna_inputs(object, "mixture", parameters, "auto"),
        "`pi_guess_range` should be a numeric vector of length 2",
        fixed = TRUE
    )
})
