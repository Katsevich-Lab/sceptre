#' Data getter functions
#'
#' The data getter functions (i.e., `get_response_matrix()`, `get_grna_matrix()`, `get_cell_covariates()`) return of a specified data field from a `sceptre_object`.
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
#'
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
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
#' `get_grna_assignments()` returns the gRNA-to-cell assignments contained within a `sceptre_object`. The output is a sparse logical matrix, with gRNAs in the rows and cells in the columns. A given entry of the matrix is set to `TRUE` if the given gRNA is assigned to the given cell (and `FALSE` otherwise).
#'
#' When using the "maximum" assignment strategy, exactly one gRNA is assigned to a given cell. In other words, each column of the gRNA-to-cell assignment matrix contains exactly one TRUE entry.
#'
#' @param sceptre_object a `sceptre_object` that has had `assign_grnas()` called on it
#' @param apply_cellwise_qc a logical value (i.e., `TRUE` or `FALSE`) indicating whether to return the gRNA-to-cell assignment matrix after cellwise QC has been applied (default `FALSE`)
#'
#' @return a sparse logical matrix containing the gRNA-to-cell assignments
#' @export
#' @examples
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#' discovery_pairs <- construct_cis_pairs(sceptre_object)
#' sceptre_object <- sceptre_object |>
#'   set_analysis_parameters(
#'     discovery_pairs = discovery_pairs,
#'     side = "left"
#'   ) |>
#'   assign_grnas(
#'     method = "mixture", parallel = TRUE, n_processors = 2
#'   ) |>
#'   run_qc()
#'
#' grna_assignment_matrix <- get_grna_assignments(
#'   sceptre_object = sceptre_object
#' )
#' grna_assignment_matrix_with_qc <- get_grna_assignments(
#'   sceptre_object = sceptre_object,
#'   apply_cellwise_qc = TRUE
#' )
get_grna_assignments <- function(sceptre_object, apply_cellwise_qc = FALSE) {
  if (!sceptre_object@functs_called[["assign_grnas"]]) {
    stop("`assign_grnas()` has not yet been called on the `sceptre_object`.")
  }
  initial_grna_assignment_list <- sceptre_object@initial_grna_assignment_list
  grna_ids <- names(initial_grna_assignment_list)
  names(initial_grna_assignment_list) <- NULL
  j <- unlist(initial_grna_assignment_list)
  increment_vector(j, -1L)
  l <- vapply(initial_grna_assignment_list, length, FUN.VALUE = integer(1))
  p <- c(0L, cumsum(l))
  mat <- Matrix::sparseMatrix(j = 1, p = c(0, 1), repr = "R")
  mat@j <- j
  mat@p <- p
  mat@Dim <- c(
    length(initial_grna_assignment_list),
    get_grna_matrix(sceptre_object) |> ncol()
  )
  rownames(mat) <- grna_ids
  if (apply_cellwise_qc) {
    if (!sceptre_object@functs_called[["run_qc"]]) {
      stop("QC has not yet been called on this sceptre_object.")
    }
    mat <- mat[, sceptre_object@cells_in_use]
  }
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
