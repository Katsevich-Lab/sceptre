# response matrix class union
#' @import Matrix
setClassUnion("response_matrix_class", c("matrix", "dgCMatrix", "dgRMatrix", "dgTMatrix"))
setClassUnion("grna_matrix_class", c("matrix", "dgCMatrix", "dgRMatrix", "dgTMatrix", "lgCMatrix", "lgRMatrix", "lgTMatrix"))


# sceptre object class
setClass("sceptre_object",
         slots = list(
           # raw data
           response_matrix = "response_matrix_class",
           grna_matrix = "grna_matrix_class",
           covariate_data_frame = "data.frame",
           covariate_matrix = "matrix",
           grna_target_data_frame = "data.frame",
           low_moi = "logical",
           user_specified_covariates = "character",
           response_names = "character",

           # analysis parameters
           discovery_pairs = "data.frame",
           positive_control_pairs = "data.frame",
           formula_object = "formula",
           side_code = "integer",
           fit_parametric_curve = "logical",
           control_group_complement = "logical",
           run_permutations = "logical",
           n_nonzero_trt_thresh = "integer",
           n_nonzero_cntrl_thresh = "integer",
           B1 = "integer", B2 = "integer", B3 = "integer",
           grna_grouping_strategy = "character",
           grna_assignment_method = "character",
           grna_assignment_hyperparameters = "list",
           multiple_testing_alpha = "numeric",
           multiple_testing_method = "character",
           cell_removal_metrics = "integer",

           # computed objects
           mitochondrial_gene = "logical",
           M_matrix = "matrix",
           n_nonzero_tot_vector = "integer",
           discovery_pairs_with_info = "data.frame",
           positive_control_pairs_with_info = "data.frame",
           negative_control_pairs = "data.frame",
           initial_grna_assignment_list = "list",
           grna_assignments_raw = "list",
           grna_assignments = "list",
           grnas_per_cell = "integer",
           cells_w_multiple_grnas = "integer",
           cells_in_use = "integer",
           n_ok_discovery_pairs = "integer",
           n_ok_positive_control_pairs = "integer",
           calibration_group_size = "integer",
           n_calibration_pairs = "integer",

           # cached objects
           response_precomputations = "list",

           # flags
           last_function_called = "character",
           functs_called = "logical",

           # results
           calibration_result = "data.frame",
           power_result = "data.frame",
           discovery_result = "data.frame"))
