#' Basic getter functions
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
#' @export
get_grna_matrix <- function(sceptre_object) {
  return(sceptre_object@grna_matrix[[1]])
}

#' Get cell covariates
#'
#' @rdname get_response_matrix
#' @export
get_cell_covariates <- function(sceptre_object) {
  return(sceptre_object@covariate_data_frame)
}

#' Get gRNA assignments
#'
#' Obtains the gRNA-to-cell assignments from a `sceptre_object`. The output is a sparse logical matrix, with gRNAs in the rows and cells in the columns. A given entry of the matrix is TRUE if the given gRNA is assigned to the given cell.
#'
#' - The assignments correspond to the original gRNA expression matrix, i.e. the expression matrix that has not yet had cellwise QC performed on it.
#' - When using the "maximum" assignment strategy, exactly one gRNA is assigned to a given cell. In other words, each column of the gRNA-to-cell assignment matrix contains exactly one TRUE entry.
#'
#' @param sceptre_object a `sceptre_object` that has had `assign_grnas()` called on it
#'
#' @return a sparse logical matrix containing the gRNA-to-cell assignments
#' @export
get_grna_assignments <- function(sceptre_object) {
  if (!sceptre_object@functs_called[["assign_grnas"]]) {
    stop("`assign_grnas()` has not yet been called on the `sceptre_object`.")
  }
  initial_grna_assignment_list <- sceptre_object@initial_grna_assignment_list
  grna_ids <- names(initial_grna_assignment_list)
  names(initial_grna_assignment_list) <- NULL
  j <- unlist(initial_grna_assignment_list)
  increment_vector(j, -1L)
  l <- sapply(initial_grna_assignment_list, length)
  p <- c(0L, cumsum(l))
  mat <- Matrix::sparseMatrix(j = 1, p = c(0, 1), repr = "R")
  mat@j <- j
  mat@p <- p
  mat@Dim <- get_grna_matrix(sceptre_object) |> dim()
  rownames(mat) <- grna_ids
  return(mat)
}

set_response_matrix <- function(sceptre_object, response_matrix) {
  sceptre_object@response_matrix[[1]] <- response_matrix
  return(sceptre_object)
}

set_grna_matrix <- function(sceptre_object, grna_matrix) {
  sceptre_object@grna_matrix[[1]] <- grna_matrix
  return(sceptre_object)
}
