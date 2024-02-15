# import from cellranger disk
import_data_from_cellranger_disk <- function(directories, moi, grna_target_data_frame, extra_covariates, directory_to_write) {
  # 1. call the corresponding ondisc function
  output <- ondisc::create_odm_from_cellranger(directories_to_load = directories,
                                               directory_to_write = directory_to_write,
                                               write_cellwise_covariates = FALSE)
  # 2 check the inputs
  check_import_data_inputs(output$gene, output$grna, grna_target_data_frame, moi, extra_covariates) |> invisible()

  # 3. update fields on the sceptre_object
  sceptre_object <- methods::new("sceptre_object")
  sceptre_object <- set_response_matrix(sceptre_object, output$gene)
  sceptre_object <- set_grna_matrix(sceptre_object, output$grna)
  if (!("vector_id" %in% colnames(grna_target_data_frame))) {
    sceptre_object@grna_target_data_frame <- grna_target_data_frame |> dplyr::mutate(grna_id = as.character(grna_id),
                                                                                     grna_target = as.character(grna_target))
  } else {
    sceptre_object@grna_target_data_frame_with_vector <- grna_target_data_frame |> dplyr::mutate(grna_id = as.character(grna_id),
                                                                                                 grna_target = as.character(grna_target),
                                                                                                 vector_id = as.character(vector_id))
  }
  sceptre_object@low_moi <- (moi == "low")
  sceptre_object@integer_id <- output$gene@integer_id
  sceptre_object@nf_pipeline <- FALSE
  # 4. devise the initial gRNA assignment list and process cellwise covariates
  cellwise_covariates <- output$cellwise_covariates
  sceptre_object@ondisc_grna_assignment_info <- list(max_grna = cellwise_covariates$grna_feature_w_max_expression,
                                                     max_grna_frac_umis = cellwise_covariates$grna_frac_umis_max_feature)
  cellwise_covariates$grna_feature_w_max_expression <- cellwise_covariates$grna_frac_umis_max_feature <- NULL
  colnames(cellwise_covariates) <- gsub(pattern = "gene", replacement = "response", fixed = TRUE, x = colnames(cellwise_covariates))
  sceptre_object@covariate_data_frame <- cellwise_covariates
  sceptre_object@covariate_names <- sort(colnames(sceptre_object@covariate_data_frame))
  # 5. initialize flags
  sceptre_object@last_function_called <- "import_data"
  sceptre_object@functs_called <- c(import_data = TRUE, set_analysis_parameters = FALSE,
                                    assign_grnas = FALSE, run_qc = FALSE, run_calibration_check = FALSE,
                                    run_power_check = FALSE, run_discovery_analysis = FALSE)
  return(sceptre_object)
}


#' Read ondisc-backed sceptre object
#'
#' Reads and initializes a `sceptre_object` from backing .odm files.
#'
#' @param sceptre_object_fp file path to a `sceptre_object.rds` file
#' @param response_odm_file_fp file path to a backing `.odm` file for the response modality
#' @param grna_odm_file_fp file path to a backing `.odm` file for the gRNA modality
#'
#' @return an ondisc-backed `sceptre_object`
#' @export
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
  sceptre_object <- set_response_matrix(sceptre_object, response_odm)
  sceptre_object <- set_grna_matrix(sceptre_object, grna_odm)
  return(sceptre_object)
}


#' Write ondisc-backed sceptre object
#'
#' Write an ondisc-backed sceptre object to disk.
#'
#' @param sceptre_object a `sceptre_object`
#' @param directory_to_write directory in which to write `sceptre_object.rds`
#'
#' @return NULL
#' @rdname read_ondisc_backed_sceptre_object
#' @export
write_ondisc_backed_sceptre_object <- function(sceptre_object, directory_to_write) {
  sceptre_object <- set_grna_matrix(sceptre_object, list())
  sceptre_object <- set_response_matrix(sceptre_object, list())
  # check that dir exists
  if (!dir.exists(directory_to_write)) dir.create(directory_to_write, recursive = TRUE)
  saveRDS(sceptre_object, paste0(directory_to_write, "/sceptre_object.rds"))
}


###########################
# PIPELINE HELPER FUNCTIONS
###########################
get_id_vect <- function(v, pod_size) {
  n_elements <- length(v)
  breaks <- round(n_elements/pod_size)
  if (breaks >= 2) {
    as.integer(cut(seq(1, n_elements), breaks))
  } else {
    rep(1L, n_elements)
  }
}

write_vector <- function(vector, file_name) {
  file_con <- file(file_name)
  writeLines(as.character(vector), file_con)
  close(file_con)
}
