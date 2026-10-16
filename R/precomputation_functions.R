# Private compatibility entry points used by older component scripts. Current
# package workers use the fitting/score modules directly.
perform_response_precomputation <- function(expressions, covariate_matrix) {
    fit_response_current(expressions, covariate_matrix)
}

perform_grna_precomputation <- function(trt_idxs, covariate_matrix, return_fitted_values) {
    fit_assignment_current(trt_idxs, covariate_matrix, return_fitted_values)
}

compute_precomputation_pieces <- function(expression_vector, covariate_matrix,
                                        fitted_coefs, theta, full_test_stat) {
    prepare_score_current(expression_vector, covariate_matrix,
        list(fitted_coefs = fitted_coefs, theta = theta), full_test_stat)
}
