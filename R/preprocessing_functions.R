check_inputs <- function(response_matrix, grna_matrix, covariate_data_frame, grna_group_data_frame, formula_object, calibration_check, response_grna_group_pairs, regression_method) {
  # 1. check column names of grna_group_data_frame
  colnames_present <- all(c("grna_id", "grna_group") %in% colnames(grna_group_data_frame))
  if (!colnames_present) {
    stop("The data frame `grna_group_data_frame` must have columns `grna_id` and `grna_group`. The `grna_group` column should specify the group to which each `grna_id` belongs.")
  }

  # 2. check for the presence of "non-targeting" in the grna_group column
  nt_present <- "non-targeting" %in% grna_group_data_frame$grna_group
  if (!nt_present) {
    stop(paste0("The string 'non-targeting' must be present in the `grna_group` column of the `grna_group_data_frame`."))
  }
  # verify also that >= 2 NT gRNAs are present when running a calibration check
  n_nt_grnas <- grna_group_data_frame |>
    dplyr::filter(grna_group == "non-targeting") |> nrow()
  if (calibration_check && (n_nt_grnas <= 1)) {
    stop("Two or more non-targeting gRNAs must be present when running a calibration check.")
  }

  # 3. verify that the row names are unique for both response and grna modalities
  response_ids <- rownames(response_matrix)
  grna_ids <- rownames(grna_matrix)
  if (length(response_ids) != length(unique(response_ids))) stop("The rownames of the `response_matrix` must be unique.")
  if (length(grna_ids) != length(unique(grna_ids))) stop("The rownames of the `grna_matrix` must be unique.")

  # 4. ensure that the ampersand symbol (&) is absent from the grna ids; ensure that no gRNA is named "non-targeting"
  problematic_grna_ids <- grep(pattern = "&", x = grna_ids)
  if (length(problematic_grna_ids) >= 1) {
    stop(paste0("The ampersand character (&) cannot be present in the gRNA IDs. The following gRNA IDs contain an ampersand: ", paste0(grna_ids[problematic_grna_ids], collapse = ", ")))
  }
  if (any(grna_ids == "non-targeting")) {
    stop("No individual gRNA can have the ID `non-targeting`. The string `non-targeting` is reserved for the `grna_group` column of the `grna_group_data_frame`.")
  }

  # 5. if the pairs to analyze have been specified...
    # i. verify that `grna_group` and `response_id` are columns
    all(c("grna_group", "response_id") %in% colnames(response_grna_group_pairs))
    # ii. check that the response ids in the `response_grna_group_pairs` data frame are a subset of the response ids
    if (!all(response_grna_group_pairs$response_id %in% response_ids)) {
      stop("The column `response_id` of the `response_grna_group_pairs` data frame must be a subset of the row names of the response expression matrix.")
    }
    # iii. check that the grna ids in the  `grna_group_data_frame` data frame are a subset of the grna ids
    if (!all(grna_group_data_frame$grna_id %in% grna_ids)) {
      stop("The column `grna_id` of the `response_grna_group_pairs` data frame must be a subset of the row names of the grna expression matrix.")
    }
    # iv. check that the `grna_group` column of the `response_grna_group_pairs` data frame is a subset of the `grna_group` column of the `grna_group_data_frame`
    if (!all(response_grna_group_pairs$grna_group %in% grna_group_data_frame$grna_group)) {
      stop("The column `grna_group` of the `response_grna_group_pairs` data frame must be a subset of the colummn `grna_group` of the `grna_group_data_frame`.")
    }

  # 6. check that there are no offsets in the formula object
  if (grepl("offset", as.character(formula_object)[2])) stop("Offsets are not currently supported in formula objects.")

  # 7. check the covariate data frame has at least one column and that the variables in the formula object are a subset of the column names of the covariate data frame
  if (ncol(covariate_data_frame) == 0) {
    stop("The global cell covariate matrix must contain at least one column.")
  }
  formula_object_vars <- all.vars(formula_object)
  check_var <- formula_object_vars %in% colnames(covariate_data_frame)
  if (!all(check_var)) {
    stop(paste0("The variables in the `formula_object` must be a subset of the columns of the `covariate_data_frame`. Check the following variables: ", paste0(formula_object_vars[!check_var], collapse = ", ")))
  }

  # 8. check type of input matrices
  check_matrix_class <- function(input_matrix, input_matrix_name, allowed_matrix_classes) {
    ok_class <- sapply(X = allowed_matrix_classes, function(mat_class) methods::is(input_matrix, mat_class)) |> any()
    if (!ok_class) {
      stop(paste0("`", input_matrix_name, "` must be an object of class ", paste0(allowed_matrix_classes, collapse = ", "), "."))
    }
  }
  check_matrix_class(response_matrix, "response_matrix", c("matrix", "dgTMatrix", "dgCMatrix", "dgRMatrix"))
  check_matrix_class(grna_matrix, "grna_matrix", c("matrix", "dgTMatrix", "dgCMatrix", "dgRMatrix", "lgTMatrix", "lgCMatrix", "lgRMatrix"))

  # 9. check for agreement in number of cells
  check_ncells <- (ncol(response_matrix) == ncol(grna_matrix)) && (ncol(response_matrix) == nrow(covariate_data_frame))
  if (!check_ncells) {
    stop("The number of cells in the `response_matrix`, `grna_matrix`, and `covariate_data_frame` must coincide.")
  }

  # 10. ensure that "non-targeting" is not a group in the pairs to analyze data frame
  if ("non-targeting" %in% unique(response_grna_group_pairs$grna_group)) {
    stop("The `response_grna_group_pairs` data frame cannot contain the gRNA group `non-targeting`.")
  }

  # 11. ensure that regression_method is nb_glm or poisson_glm
  if (!(regression_method %in% c("nb_glm", "poisson_glm"))) {
    stop("`regression_method` should be either `nb_glm` or `poisson_glm`.")
  }

  # 12.

  return(NULL)
}


set_matrix_accessibility <- function(matrix_in, make_row_accessible = TRUE) {
  # if a logical matrix, convert to the corresponding numeric matrix; consider more efficient implementation later
  if (methods::is(matrix_in, "lgRMatrix")) {
    attr(matrix_in, "class") <- "dgRMatrix"
    matrix_in@x <- rep(1.0, length(matrix_in@j))
  }
  if (methods::is(matrix_in, "lgCMatrix")) {
    attr(matrix_in, "class") <- "dgCMatrix"
    matrix_in@x <- rep(1.0, length(matrix_in@i))
  }
  if (methods::is(matrix_in, "lgTMatrix")) {
    attr(matrix_in, "class") <- "dgTMatrix"
    matrix_in@x <- rep(1.0,  length(matrix_in@j))
  }

  if (methods::is(matrix_in, "dgRMatrix") && make_row_accessible) {
    out <- matrix_in
  } else if (methods::is(matrix_in, "dgCMatrix") && !make_row_accessible) {
    out <- matrix_in
  } else {
    # basic info: get matrix_in, n_responses, response_ids, and n_cells
    n_responses <- nrow(matrix_in)
    response_ids <- rownames(matrix_in)
    n_cells <- ncol(matrix_in)
    if (methods::is(matrix_in, "dgTMatrix")) {
      i <- matrix_in@i
      j <- matrix_in@j
    } else if (methods::is(matrix_in, "matrix")) {
      matrix_in <- methods::as(matrix_in, "TsparseMatrix")
      i <- matrix_in@i
      j <- matrix_in@j
    } else if (methods::is(matrix_in, "dgCMatrix")) {
      i <- matrix_in@i
      j <- convert_pointer_to_index_vector_v2(matrix_in@p)
    } else if (methods::is(matrix_in, "dgRMatrix")) {
      i <- convert_pointer_to_index_vector_v2(matrix_in@p)
      j <- matrix_in@j
    } else {
      stop("Class not recognized.")
    }
    x <- matrix_in@x

    # perform the sorting operation
    sort_in_place <- methods::is(matrix_in, "matrix") || methods::is(matrix_in, "dgTMatrix")
    if (sort_in_place) {
      l <- list(i = i, j = j, x = x)
      dt <- data.table::setDT(l) # assign by reference
    } else {
      dt <- data.table::data.table(i = i, j = j, x = x) # deep copy
    }

    # finally, initialize the output
    if (make_row_accessible) {
      data.table::setorderv(x = dt, cols = c("i", "j"))
      p <- obtain_pointer_vector(i = dt$i, dim = n_responses)
      out <- Matrix::sparseMatrix(j = 1, p = c(0, 1), x = 1, repr = "R")
      out@j <- dt$j
    } else {
      data.table::setorderv(x = dt, cols = c("j", "i"))
      p <- obtain_pointer_vector(i = dt$j, dim = n_cells)
      out <- Matrix::sparseMatrix(i = 1, p = c(0, 1), x = 1, repr = "C")
      out@i <- dt$i
    }
    out@p <- p
    out@x <- dt$x
    out@Dim <- c(n_responses, n_cells)
    rownames(out) <- response_ids
  }
  return(out)
}


convert_pointer_to_index_vector_v2 <- function(p) {
  l <- diff(p)
  rep(x = seq(0, length(l) - 1L), times = l)
}


#' Convert covariate DF to design matrix
#'
#' Converts the `covariate_data_frame` to a design matrix using the `formula_object`.
#'
#' @param covariate_data_frame the input covariate data frame
#' @param formula_object the input formula object
#'
#' @return the design matrix
#' @noRd
convert_covariate_df_to_design_matrix <- function(covariate_data_frame, formula_object) {
  global_cell_covariates_new <- stats::model.matrix(object = formula_object, data = covariate_data_frame)
  # verify that the global cell covariates are OK after transformation
  for (col_name in colnames(global_cell_covariates_new)) {
    vect <- global_cell_covariates_new[,col_name]
    if (any(vect == -Inf) || any(vect == Inf) || any(is.na(vect))) {
      stop(paste0("The column `", col_name, "` of the `covariate_data_frame` after the `formula object` has been applied contains entries that are -Inf, Inf, or NA. Remove these entries."))
    }
  }
  # verify that matrix is not rank-deficient
  rank_matrix <- Matrix::rankMatrix(global_cell_covariates_new)
  if (rank_matrix != ncol(global_cell_covariates_new)) {
    stop("The `formula_object` contains redundant information. This often occurs when one variable is `nested` inside another. For example, if the data frame contains a column `lane` with levels `lane_1`, `lane_2`, `lane_3`, and `lane_4` and a column `batch` with levels `batch_1` and `batch_2`, and if `lane_1` and `lane_2` are contained entirely within `batch_1` while `lane_3` and `lane_4` are contained entirely within `batch`2`, then `batch` is redundant. In this case `batch` should be removed from the formula object.")
  }
  return(global_cell_covariates_new)
}


assign_grnas_to_cells_lowmoi_v2 <- function(grna_matrix, grna_group_data_frame, calibration_check, n_calibration_pairs) {
  # get the idx of the maximum of each column
  grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = FALSE)
  grna_idx <- compute_colwise_max(i = grna_matrix@i, p = grna_matrix@p, x = grna_matrix@x, n_cells = ncol(grna_matrix))
  # grna_idx <- apply(X = grna_matrix, MARGIN = 2, which.max)
  grna_rownames <- factor(rownames(grna_matrix))
  indiv_grna_id_assignments <- grna_rownames[grna_idx]
  grna_group_assignments <- grna_group_data_frame$grna_group[match(x = indiv_grna_id_assignments,
                                                                   table = grna_group_data_frame$grna_id)]
  grna_groups <- grna_group_data_frame$grna_group

  # obtain all_nt_idxs
  out <- list()

  # add the targeting grna groups if doing a discovery analysis or n_calibration_pairs not specified
  if (!calibration_check || is.null(n_calibration_pairs)) {
    unique_grna_groups <- unique(grna_groups)
    unique_grna_groups <- unique_grna_groups[unique_grna_groups != "non-targeting"]
    grna_group_idxs <- lapply(unique_grna_groups, function(unique_grna_group) {
      which(grna_group_assignments == unique_grna_group)
    }) |> stats::setNames(unique_grna_groups)
    out$grna_group_idxs <- grna_group_idxs
  }

  # if running a discovery analysis, add the indices of all NT gRNAs
  if (!calibration_check) {
    out$all_nt_idxs <- which(grna_group_assignments == "non-targeting")
  }

  # if running a calibration check, add the indices of the individual NT gRNAs relative to all_nt_idxs
  if (calibration_check) {
    nt_grnas <- grna_group_data_frame |>
      dplyr::filter(grna_group == "non-targeting") |> dplyr::pull("grna_id")
    nt_grna_idxs <- lapply(nt_grnas, function(nt_grna) {
      which(nt_grna == indiv_grna_id_assignments)
    }) |> stats::setNames(nt_grnas)
    all_nt_idxs <- stats::setNames(unlist(nt_grna_idxs), NULL)
    n_cells_per_nt <- sapply(nt_grna_idxs, length)
    stop <- cumsum(n_cells_per_nt)
    start <- c(0L, stop[-length(stop)]) + 1L
    indiv_nt_grna_idxs <- lapply(seq(1, length(nt_grnas)), function(i) {
      seq(start[i], stop[i])
    }) |> stats::setNames(nt_grnas)
    out$all_nt_idxs <- all_nt_idxs
    out$indiv_nt_grna_idxs <- indiv_nt_grna_idxs
  }

  return(out)
}


get_synthetic_idxs_lowmoi <- function(grna_assignments, B, calibration_check, undercover_group_size = NULL) {
  if (calibration_check) {
    indiv_nt_sizes <- sapply(grna_assignments$indiv_nt_grna_idxs, length) |> sort(decreasing = TRUE)
    M <- sum(indiv_nt_sizes[seq(1, undercover_group_size)])
    n_control_cells <- length(grna_assignments$all_nt_idxs)
    out <- fisher_yates_samlper(n_tot = n_control_cells, M = M, B = B)
  } else { # discovery
    grna_group_sizes <- sapply(grna_assignments$grna_group_idxs, length)
    grna_group_sizes <- grna_group_sizes[grna_group_sizes != 0L]
    range_grna_group_sizes <- range(grna_group_sizes)
    n_control_cells <- length(grna_assignments$all_nt_idxs)
    out <- hybrid_fisher_iwor_sampler(N = n_control_cells,
                                      m = range_grna_group_sizes[1],
                                      M = range_grna_group_sizes[2],
                                      B = B)
  }
  return(out)
}


harmonize_arguments <- function(return_resampling_dist, fit_skew_normal) {
  if (return_resampling_dist) {
    assign(x = "B2", value = 0L, inherits = TRUE)
    assign(x = "B3", value = 0L, inherits = TRUE)
    assign(x = "fit_skew_normal", value = FALSE, inherits = TRUE)
  }
  if (!fit_skew_normal) {
    assign(x = "B2", value = 0L, inherits = TRUE)
  }

  return (NULL)
}


#' Compute cell covariates
#'
#' The \code{compute_cell_covariates()} function takes as input a feature-by-cell matrix and computes the cell covariates for that matrix.
#'
#' @param matrix_in a matrix with features (e.g., genes, gRNAs, or proteins) in the rows and cells in the columns. The matrix can be a standard (dense) matrix or a sparse matrix of class \code{dgCMatrix}, \code{dgRMatrix}, or \code{dgTMatrix}. Row names are optional.
#'
#' @return a data frame containing the following cell covariates.
#' \itemize{
#' \item{\code{n_nonzero}}: the number of features expressed in the cell
#' \item{\code{n_umi}}: the total number of UMIs sequenced in the cell
#' \item{\code{p_mito}}: the fraction of UMIs sequenced in the cell that maps to mitochondrial genes
#' }
#' @export
#' @note The variable \code{p_mito} is computed only if feature names are provided and some feature names begin with "MT-"; these features are assumed to be mitochondrial genes.
#'
#' @examples
#' data(response_matrix_lowmoi)
#' cell_covariates <- compute_cell_covariates(response_matrix_lowmoi)
compute_cell_covariates <- function(matrix_in) {
  # make response matrix column accessible
  matrix_in <- set_matrix_accessibility(matrix_in, make_row_accessible = FALSE)
  # get MT gene idxs
  mt_gene_idxs <- grep(pattern = "^MT-", x = rownames(matrix_in))
  compute_p_mito <- length(mt_gene_idxs) >= 1
  # call the low-level function
  out <- compute_cell_covariates_cpp(i = matrix_in@i,
                                     p = matrix_in@p,
                                     x = matrix_in@x,
                                     n_genes = nrow(matrix_in),
                                     n_cells = ncol(matrix_in),
                                     mt_gene_idxs = mt_gene_idxs,
                                     compute_p_mito = compute_p_mito)
  ret <- data.frame(n_nonzero = out$n_nonzero,
                    n_umis = out$n_umi)
  if (compute_p_mito) ret$p_mito <- out$p_mito
  return(ret)
}


compute_regression_ses <- function(covariate_matrix_nt, w) {
 info_mat <- t(covariate_matrix_nt * w) %*% covariate_matrix_nt
 ses <- sqrt(diag(solve(info_mat)))
 return(ses)
}
