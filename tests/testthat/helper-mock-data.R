#####################################################
## Functions for making mock data for unit testing ##
#####################################################

#' Makes a mock UMI count matrix for testing.
#'
#' This can replicate either a response matrix or a gRNA expression matrix.
#'
#' The specific patterns that appear:
#' * constant
#' * all but one value constant
#' * linear sequences
#' * random
#'
#' There are 13 deterministic patterns tested, so with at least one random pattern we
#' require at least 14 in the dimension in which the patterns go; in other words, if
#' \code{patterns_at_col_level = TRUE} then \code{num_cells >= 14} is required (since
#' the patterns are at the column / cell level), and if instead \code{patterns_at_col_level = FALSE}
#' then \code{num_genes >= 14} is required.
#'
#'
#' @param num_genes, @param num_cells : the number of rows and columns in the resulting data, respectively
#' @param rowname_base : the rownames are of the form `rowname_base` + "_" + (1:num_genes)
#' @param patterns_at_col_level : if TRUE the columns contain the above patterns, else the rows do
#' @param big : a constant representing a large data value
#' @param seed : for the random part of the data
#'
#' @return a sparse matrix of class \code{dgTMatrix} mocking a response matrix.
#' @export
#'
#' @examples
#' # mock response UMI matrix with 5 rows and 20 columns
#' # the structured patterns are for each cell, so e.g. cell 1 has UMI counts of all 0
#' make_mock_count_matrix(num_genes = 5, num_cells = 20, rowname_base = "gene")
#'
#' # mock response UMI matrix with 20 rows and 5 columns
#' # the structured patterns are now for each gene, so e.g. gene 1 has UMI counts of all 0
#' make_mock_count_matrix(num_genes = 20, num_cells = 5, rowname_base = "gene", patterns_at_col_level = FALSE)
#'
#' # we can also mock gRNA count matrices
#' make_mock_count_matrix(num_genes = 7, num_cells = 20, rowname_base = "grna")
make_mock_count_matrix <- function(num_genes, num_cells, rowname_base, patterns_at_col_level = TRUE, big = 1000, seed = NULL) {
  if(is.null(seed)) {
    seed = 101000
  }
  set.seed(seed)
  dims <- if(patterns_at_col_level) c(num_genes, num_cells) else c(num_cells, num_genes)
  if(dims[2] < 14) {
    stop("There are 14 required patterns and the provided dimensions are too small.")
  }
  data_systematic <- cbind(
    0, 1, big, # constant columns
    rep(c(0,1), c(dims[1]-1, 1)),  rep(c(0, 1), c(1, dims[1]-1)), # all but one element constant
    rep(c(big, 0), c(dims[1]-1, 1)) ,rep(c(big, 1), c(dims[1]-1, 1)),
    rep(c(0, big), c(dims[1]-1, 1)), rep(c(1, big), c(dims[1]-1, 1)),
    0:(dims[1] - 1), dims[1]:1, # sequence columns
    (0:(dims[1] - 1)) * big, (dims[1]:1) * big
  ) |>
    `colnames<-`(NULL) # the use of `big` adds a single column name
  data_random <- matrix(
    rpois(dims[1] * (dims[2] - ncol(data_systematic)), 1),
    nrow = dims[1]
  )
  mock_response_matrix <- cbind(data_systematic, data_random)

  if(!patterns_at_col_level) {
    mock_response_matrix <- t(mock_response_matrix)
  }
  rownames(mock_response_matrix) <- paste0(rowname_base, "_", 1:nrow(mock_response_matrix))  # rows are always genes
  return(as(mock_response_matrix, "TsparseMatrix"))
}

# `sample` has unwanted behavior here where e.g. `sample(2:2, 1)` can return 1
# even though only returning 2 is desired.
sample_as_vec <- function(x) {
  if(length(x) == 1)
    return(x)
  else {
    return(sample(x, 1))
  }
}

#' Helper function to make a factor with random numbers of entries per level.
#'
#' @param num_factor_levels : int, the number of distinct values for this factor to have
#' @param num_entries : int, the length of the resulting vector
#' @param level_name_base : string, all values start with this + "_".
#' @param return_factor : bool, if \code{FALSE} the result is a character vector while if \code{TRUE}
#' it is actually a factor.
#' @param shuffle : if TRUE then the values of the vector are shuffled
#'
#' @return
#' @export
#'
#' @examples
make_random_factor <- function(num_factor_levels, num_entries, level_name_base = "lev", return_factor = TRUE, shuffle = FALSE) {
  lev_counts <- numeric(num_factor_levels)
  lev_counts[1] <- sample_as_vec(2:(num_entries - 2 * (num_factor_levels - 1)))
  if(num_factor_levels > 2) { # middle levels get random counts
    for(j in 2:(num_factor_levels-1)) {
      # we need to save 2 * (# levels remaining) entries for the future levels, and
      # we've allocated `sum(lev_counts)` so far
      lev_counts[j] <- sample_as_vec(2:(num_entries - 2 * (num_factor_levels - j) - sum(lev_counts)))
    }
  }
  lev_counts[num_factor_levels] <- num_entries - sum(lev_counts) # needs to sum to `num_entries`
  lev_names <- paste0(level_name_base, "_", 1:num_factor_levels)
  values <- rep(lev_names, lev_counts)
  if(shuffle) {
    values <- sample(values, num_entries, FALSE)
  }
  if(return_factor) {
    values <- factor(values, levels = lev_names)
  }
  return(values)
}

#' Mock cell-level covariates
#'
#' @param rep_level_counts : this is for mocking a replicate or batch feature. \code{length(rep_level_counts)}
#' gives the number of replicates, and each entry is the number of cells within that replicate.
#' \code{sum(rep_level_counts)} is the total number of cells. For example,
#' \code{rep_level_counts = c(2,4,6)} will result in data for 12 cells total, with
#' replicate sizes of 2, 4, and 6 respectively.
#' @param add_numeric : integter, if positive then that many columns of
#' iid N(0,1) data are added, mocking continuous features.
#' @param add_count integer, if positive then that many columns of
#' iid Pois(1) data are added, mocking count features.
#' @param add_factor integer, if positive then that many random factor features
#' are added. For each random factor, the number of levels is sampled uniformly from
#' \code{2:max_num_factor_levels}, and then the number of cells per level are filled
#' in randomly, with the constraint that no level has fewer than 2 cells.
#'
#' @param max_num_factor_levels : integer, the maximum number of levels for each random
#' factor. Ignored if \code{add_factor = 0}.
#' @param seed
#'
#' @return a data.frame with the specified covariates, with each row corresponding to a cell.
#' @export
#'
#' @examples
#' # just doing replicate for 20 cells
#' make_mock_extra_covariates(rep_level_counts = c(4, 13, 3))
#' # simulating features of all types. With this seed, two of the factors ended up with
#' # just 2 levels, while the third ended up with 4 levels.
#' make_mock_extra_covariates(rep_level_counts = c(4, 13, 3), add_numeric = 2, add_count = 1,
#' add_factor = 3, max_num_factor_levels = 4, seed = 111)
make_mock_extra_covariates <- function(rep_level_counts, add_numeric = 0, add_count = 0, add_factor = 0,
                                       max_num_factor_levels = 6, seed = NULL) {
  if(is.null(seed)) {
    seed = 101001
  }
  set.seed(seed)

  num_cells <- sum(rep_level_counts)
  if(num_cells < 12) {
    stop("At least 12 cells are required.")
  }
  if(add_factor > 0 && max_num_factor_levels * 2 > num_cells) {
    stop("`max_num_factor_levels` must be decreased so that at least 2 cells per level are possible.")
  }
  if(max_num_factor_levels < 2) {
    stop("`max_num_factor_levels` must be at least 2.")
  }
  covariates <- data.frame(
    replicate = rep(paste0("rep_", 1:length(rep_level_counts)), rep_level_counts) |>
      factor(levels = paste0("rep_", 1:length(rep_level_counts)))
  )
  if(add_numeric > 0) {
    for(i in 1:add_numeric) {
      covariates[[paste0("numeric_", i)]] <- rnorm(num_cells)
    }
  }
  if(add_count > 0) {
    for(i in 1:add_count) {
      covariates[[paste0("count_", i)]] <- rpois(num_cells, 1)
    }
  }
  if(add_factor > 0) {
    for(i in 1:add_factor) {
      # at least 12 cells required so guaranteed to get at least 2 cells / level
      num_levs <- sample(2:max_num_factor_levels, 1)
      covariates[[paste0("factor_", i)]] <- make_random_factor(
        num_factor_levels = num_levs, num_entries = num_cells,
        level_name_base = paste0("factor_", i),
        return_factor = TRUE, shuffle = TRUE
      )
    }
  }
  return(covariates)
}

#' Mock grna_id to grna_target data.frame
#'
#' @param grna_labels : a vector of unique grna id's
#' @param num_targets : the number of groups to form out of those
#'
#' @return : a data.frame with two columns, one giving the provided \code{grna_labels},
#' and the other giving the newly created \code{grna_target}.
#' @export
#'
#' @examples
#' make_mock_grna_target_data(paste0("grna_", 1:12), 4)
make_mock_grna_target_data <- function(grna_labels, num_targets, seed = NULL) {
  if(is.null(seed)) {
    seed = 10102
  }
  set.seed(seed)

  data.frame(
    grna_id = grna_labels,
    grna_target = make_random_factor(
      num_factor_levels = num_targets, num_entries = length(grna_labels),
      level_name_base = "grna_target", return_factor = FALSE
    )
  )
}


#############
## Example ##
#############
# num_cells <- 15
# num_genes <- 18
# num_grnas <- 5
# num_distinct_targets <- 4
# rep_levels <- c(4, 8, 3); stopifnot(sum(rep_levels) == num_cells)
# (mock_response <- make_mock_count_matrix(num_genes = num_genes, num_cells = num_cells, rowname_base = "gene", patterns_at_col_level = FALSE))
# (mock_grna <- make_mock_count_matrix(num_genes = num_grnas, num_cells = num_cells, rowname_base = "grna", patterns_at_col_level = TRUE))
# (mock_covariates <- make_mock_extra_covariates(rep_level_counts = rep_levels))
# ## if we want a fancier one
# # make_mock_extra_covariates(rep_level_counts = rep_levels, add_numeric = 2, add_count = 1, add_factor = 3, max_num_factor_levels = 4)
# (mock_grna_target_data_frame <- make_mock_grna_target_data(rownames(mock_response), num_distinct_targets))

