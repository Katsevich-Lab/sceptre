# import from cellranger disk
import_data_from_cellranger_disk <- function(directories, moi, grna_target_data_frame, extra_covariates, directory_to_write) {
  # 1. call the corresponding ondisc function
  output <- ondisc::create_odm_from_cellranger(directories_to_load = directories,
                                               directory_to_write = directory_to_write,
                                               write_cellwise_covariates = FALSE)
  # 2. update fields on the sceptre_object
  sceptre_object <- methods::new("sceptre_object")
  sceptre_object@response_matrix <- output$gene
  sceptre_object@grna_matrix <- output$grna
  sceptre_object@grna_target_data_frame <- grna_target_data_frame |> dplyr::mutate(grna_id = as.character(grna_id), grna_target = as.character(grna_target))
  sceptre_object@low_moi <- (moi == "low")
  sceptre_object@integer_id <- output$gene@integer_id
  # 3. devise the initial gRNA assignment list and process cellwise covariates
  cellwise_covariates <- output$cellwise_covariates
  sceptre_object@ondisc_grna_assignment_info <- list(max_grna = cellwise_covariates$grna_feature_w_max_expression,
                                                     max_grna_frac_umis = cellwise_covariates$grna_frac_umis_max_feature)
  cellwise_covariates$grna_feature_w_max_expression <- cellwise_covariates$grna_frac_umis_max_feature <- NULL
  colnames(cellwise_covariates) <- gsub(pattern = "gene", replacement = "response", fixed = TRUE, x = colnames(cellwise_covariates))
  sceptre_object@covariate_data_frame <- cellwise_covariates
  # 4. initialize flags
  sceptre_object@last_function_called <- "import_data"
  sceptre_object@functs_called <- c(import_data = TRUE, set_analysis_parameters = FALSE,
                                    assign_grnas = FALSE, run_qc = FALSE, run_calibration_check = FALSE,
                                    run_power_check = FALSE, run_discovery_analysis = FALSE)
  return(sceptre_object)
}


read_ondisc_backed_sceptre_object <- function(sceptre_object_fp, response_odm_file_fp, grna_odm_file_fp) {
  # read in objects
  sceptre_object <- readRDS(sceptre_object_fp)
  response_odm <- ondisc::initialize_odm_from_backing_file(response_odm_file_fp)
  grna_odm <- ondisc::initialize_odm_from_backing_file(grna_odm_file_fp)
  # check for concordance of integer ids
  if (sceptre_object@integer_id != response_odm@integer_id) {
    stop("The `sceptre_object` and `response_odm` have distinct IDs. The `sceptre_object` likely is associated with a different backing .odm file.")
  }
  if (sceptre_object@integer_id != grna_odm@integer_id) {
    stop("The `sceptre_object` and `grna_odm` have distinct IDs. The `sceptre_object` likely is associated with a different backing .odm file.")
  }
  # update response_odm and grna_odm and return
  sceptre_object@response_matrix <- response_odm
  sceptre_object@grna_matrix <- grna_odm
  return(sceptre_object)
}


write_ondisc_backed_sceptre_object <- function(sceptre_object, sceptre_object_fp) {
  sceptre_object@grna_matrix <- matrix()
  sceptre_object@response_matrix <- matrix()
  saveRDS(sceptre_object, sceptre_object_fp)
}
