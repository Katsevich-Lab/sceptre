utils::globalVariables(c("n_nonzero_trt", "n_nonzero_cntrl", "pair_str", "assignment", "g",
                         "multiple_grnas", "x", "grna_id", "grna_expressions_bin", "bin_counts",
                         "significant", "lab", "p_values", "pass_qc", "fraction_cells_removed",
                         "grna_id", "pass_qc", "grna_group", "p_value", "gene_id", "response_id",
                         "log_2_fold_change", "reject", "y", "grna_target", "gene_position_data_frame_grch38",
                         "any_pass_qc"))
#' sceptre
#'
#' `sceptre` is an R package for single-cell CRISPR screen data analysis that emphasizes statistical rigor, computational efficiency, and ease of use.
#'
#' @useDynLib sceptre, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import BH
#' @docType package
#' @name sceptre
#' @examples
#' \dontrun{
#' ##########################
#' # Low-MOI CRISPRko example
#' ##########################
#' # 1. create the sceptre object
#' data("lowmoi_example_data")
#' sceptre_object <- import_data(
#' response_matrix = lowmoi_example_data$response_matrix,
#' grna_matrix = lowmoi_example_data$grna_matrix,
#' extra_covariates = lowmoi_example_data$extra_covariates,
#' grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
#' moi = "low")
#' print(sceptre_object)
#'
#' # 2. set the analysis parameters
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' discovery_pairs <- construct_trans_pairs(sceptre_object = sceptre_object,
#' positive_control_pairs = positive_control_pairs,
#' pairs_to_exclude = "pc_pairs")
#'
#' sceptre_object <- set_analysis_parameters(
#' sceptre_object = sceptre_object,
#' discovery_pairs = discovery_pairs,
#' positive_control_pairs = positive_control_pairs)
#' print(sceptre_object)
#'
#' # 3. assign grnas
#' plot_grna_count_distributions(sceptre_object)
#' sceptre_object <- sceptre_object |> assign_grnas()
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 4. run qc
#' plot_covariates(sceptre_object, p_mito_threshold = 0.075)
#' sceptre_object <- sceptre_object |> run_qc(p_mito_threshold = 0.075)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 5. run the calibration check
#' sceptre_object <- run_calibration_check(sceptre_object, parallel = TRUE)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 6. run power check
#' sceptre_object <- run_power_check(sceptre_object, parallel = TRUE)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 7. run discovery analysis
#' sceptre_object <- run_discovery_analysis(sceptre_object, parallel = TRUE)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 8. write results
#' write_outputs_to_directory(sceptre_object = sceptre_object, "~/sceptre_outputs_lowmoi/")
#'
#' ##########################
#' # High-MOI CRISPRi example
#' ##########################
#' # 1. create the sceptre object from cellranger output
#' directories <- paste0(system.file("extdata", package = "sceptre"),
#' "/highmoi_example/gem_group_", 1:2)
#' data(grna_target_data_frame_highmoi)
#' sceptre_object <- import_data_from_cellranger(directories = directories,
#' moi = "high",
#' grna_target_data_frame = grna_target_data_frame_highmoi)
#' print(sceptre_object)
#'
#' # 2. set the analysis parameters
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' discovery_pairs <- construct_cis_pairs(sceptre_object,
#' positive_control_pairs = positive_control_pairs,
#' distance_threshold = 5e6)
#'
#' sceptre_object <- set_analysis_parameters(
#' sceptre_object = sceptre_object,
#' discovery_pairs = discovery_pairs,
#' positive_control_pairs = positive_control_pairs,
#' side = "left")
#' print(sceptre_object)
#'
#' # 3. assign grnas
#' plot_grna_count_distributions(sceptre_object)
#' sceptre_object <- sceptre_object |> assign_grnas(parallel = TRUE)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 4. run qc
#' plot_covariates(sceptre_object, p_mito_threshold = 0.075)
#' sceptre_object <- sceptre_object |> run_qc(p_mito_threshold = 0.075)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 5. run the calibration check
#' sceptre_object <- run_calibration_check(sceptre_object, parallel = TRUE)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 6. run the power check
#' sceptre_object <- run_power_check(sceptre_object, parallel = TRUE)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 7. run discovery analysis
#' sceptre_object <- run_discovery_analysis(sceptre_object, parallel = TRUE)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 8. write results
#' write_outputs_to_directory(sceptre_object = sceptre_object, "~/sceptre_outputs_highmoi/")
#'
#' #####################
#' # out-of-core example
#' #####################
#' # 0. set file paths
#' base_dir <- "/Users/timbarry/research_offsite/external/replogle-2022/raw/rd7/rpe1_other"
#' directories <- list.files(base_dir, full.names = TRUE)[1:3]
#' directory_to_write <- "/Users/timbarry/research_offsite/external/replogle-2022/processed/rd7/small"
#' moi <- "low"
#' grna_target_data_frame <- readRDS("/Users/timbarry/research_offsite/external/replogle-2022/raw/rd7/grna_target_data_frame.rds")
#'
#' # 1. import data
#' sceptre_object <- import_data_from_cellranger(
#'    directories = directories,
#'    directory_to_write = directory_to_write,
#'    moi = moi,
#'    grna_target_data_frame = grna_target_data_frame,
#'    use_ondisc = TRUE)
#'
#' # 2. set analysis parameters
#' discovery_pairs <- construct_trans_pairs(
#'   sceptre_object = sceptre_object) |> dplyr::sample_n(1000000)
#' sceptre_object <- set_analysis_parameters(
#'   sceptre_object = sceptre_object,
#'   discovery_pairs = discovery_pairs)
#'
#' # 3. set analysis parameters
#' sceptre_object <- sceptre_object |> set_analysis_parameters(
#'   discovery_pairs = discovery_pairs
#' )
#'
#' # 4. write sceptre_object to disk
#' write_ondisc_backed_sceptre_object(
#'   sceptre_object,
#'   paste0(directory_to_write, "/sceptre_object.rds")
#' )
#'
#' # 5. read sceptre_object from disk
#' sceptre_object <- read_ondisc_backed_sceptre_object(
#'   paste0(directory_to_write, "/sceptre_object.rds"),
#'   paste0(directory_to_write, "/gene.odm"),
#'   paste0(directory_to_write, "/grna.odm")
#' )
#'
#' # 6. assign gRNAs
#' sceptre_object <- sceptre_object |>
#'   assign_grnas(
#'     method = "mixture",
#'     parallel = TRUE
#'   )
#'
#' # 7. run qc
#' sceptre_object <- run_qc(sceptre_object)
#'
#' }
NULL
