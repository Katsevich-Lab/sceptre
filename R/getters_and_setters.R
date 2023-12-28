#' Getter functions
#'
#' The getter functions (i.e., `get_*`) return of a specified field from a `sceptre_object`.
#'
#' @param sceptre_object a `sceptre_object`
#'
#' @return
#' `get_response_matrix()` returns the response matrix contained within a `sceptre_object`.
#'
#' `get_grna_matrix()` returns the gRNA matrix contained within a `sceptre_object`.
#'
#' `get_cell_covariates()` returns the cell covariate data frame contained within a `sceptre_object`.
#' @export
#' @examples
#' # 1. create a sceptre_object
#' directories <- paste0(system.file("extdata", package = "sceptre"),
#' "/highmoi_example/gem_group_", 1:2)
#' data(grna_target_data_frame_highmoi)
#' sceptre_object <- import_data_from_cellranger(
#'   directories = directories,
#'   moi = "high",
#'   grna_target_data_frame = grna_target_data_frame_highmoi
#' )
#'
#' # 2. extract data fields from the sceptre_object
#' response_matrix <- get_response_matrix(sceptre_object)
#' grna_matrix <- get_grna_matrix(sceptre_object)
#' cell_covariates <- get_cell_covariates(sceptre_object)
get_response_matrix <- function(sceptre_object) {
  return(sceptre_object@response_matrix[[1]])
}

#' Get gRNA matrix
#'
#' @rdname get_response_matrix
get_grna_matrix <- function(sceptre_object) {
  return(sceptre_object@grna_matrix[[1]])
}

#' Get cell covariates
#'
#' @rdname get_response_matrix
get_cell_covariates <- function(sceptre_object) {
  return(sceptre_object@covariate_data_frame)
}

set_response_matrix <- function(sceptre_object, response_matrix) {
  sceptre_object@response_matrix[[1]] <- response_matrix
  return(sceptre_object)
}

set_grna_matrix <- function(sceptre_object, grna_matrix) {
  sceptre_object@grna_matrix[[1]] <- grna_matrix
  return(sceptre_object)
}
