##########################
# READ AND WRITE FUNCTIONS
##########################
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
