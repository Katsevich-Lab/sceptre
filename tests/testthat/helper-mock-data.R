######################################################################################################
##                          Functions for making mock data for unit testing                         ##
######################################################################################################

# This function takes in a vector like c(1,5,3) and returns a data.frame with
# columns "grna_target" and "chr" with rows representing different genomic targets
# and what chromosome they are on.
# The input determines how many different guides there will be for each target, so
# in the example `num_guides_per_target = c(1, 5, 3)` we will have one target with one guide,
# one target with 5 guides, and one target with 3 guides. These will all be on different chromosomes.
# This is then repeated except all are on the same chromosome, so this example will have 4 chromosomes, with
# the final chromosome having length(c(1, 5, 3)) targets with 1+5+3 guides total for it.
.make_mock_target_and_chr_only <- function(num_guides_per_target) {
  targets <- rep(paste0("t", 1:length(num_guides_per_target)), num_guides_per_target)
  chrs <- rep(paste0("c", 1:(length(num_guides_per_target) + 1)), c(num_guides_per_target, sum(num_guides_per_target)))
  grna_target_data_simple <- data.frame(
    grna_target = c(
      paste0("t1_c", rep(1:length(num_guides_per_target), num_guides_per_target)),
      paste0(targets, "_c", length(num_guides_per_target) + 1)
    ),
    chr = chrs
  )
  return(grna_target_data_simple)
}

# This function takes the output of `.make_mock_target_and_chr_only` and adds in columns "start" and "end"
# representing where the target starts and stops on the given chromosome.
# This function is not meant to be directly called. It is assumed `length(chr_distances) == length(chr_starts)`.
# The result is a data.frame with columns "grna_target", "chr", "start", and "end", consisting of
# `target_and_chr_data` copied once per entry of `chr_distances` with each copy having the corresponding
# values of `chr_distances` and `chr_starts`.
# The input's `grna_target` and `chr` are modified to indicate which distance they correspond to.
# Note: in the semantic naming of targets and chromosomes, the index of the distance is used rather
# than the value of the distance, because there could be repeated distances.
.make_mock_cross_target_chr_with_distances <- function(target_and_chr_data, chr_distances, chr_starts) {
  target_chr_and_dist_data <- lapply(seq_along(chr_distances), function(i) {
    target_and_chr_data |>
      dplyr::transmute(
        grna_target = paste0(grna_target, "_d", i),
        chr = paste0(chr, "_d", i),
        start = chr_starts[i], end = chr_starts[i] + chr_distances[i]
      )
  }) |>
    do.call(what = rbind)
  return(target_chr_and_dist_data)
}

#' Make mock \code{grna_target_data_frame}s
#'
#' This function makes a mock grna_target_data_frame in the spirit of
#' \code{data(grna_target_data_frame_highmoi)} and \code{data(grna_target_data_frame_lowmoi)}.
#' It has columns \code{grna_id}, \code{grna_target}, \code{chr}, \code{start}, and \code{end}.
#'
#' It builds the mock data in three steps.
#'
#' The first step takes in \code{num_guides_per_target} and makes a data.frame with one grna target
#' per entry of \code{num_guides_per_target}, repeated according to the values in that vector. Each
#' of these grna targets is on its own chromosome. This data.frame is then repeated except all targets
#' are now on the same chromosome.
#'
#' The second step adds in \code{start} and \code{end}. For every unique value of \code{chr_distances}
#' (and the corresponding value of \code{chr_starts}), the previous step is duplicated and then this copy's
#' \code{start} and \code{end} are filled in with the corresponding values of \code{chr_distances} and
#' \code{chr_starts}.
#'
#' Finally, the \code{grna_id} column is added in and optionally NT guides are appended.
#'
#' @param num_guides_per_target : a vector (length >=1) of positive ints giving the number of distinct grnas per target,
#' before the doubling in the chromosome addition step and the multiplication in the \code{start} and \code{end}
#' steps as described above.
#' @param chr_distances : a vector (length >=1) of non-negative ints giving the length on the chromosome of the various
#' guides.
#' @param chr_starts : a vector of positive ints either of length \code{chr_distances} or length 1 (in which case that
#' value is recycled) determining the start locations on the chromosome of the various guides.
#' @param num_nt_guides : an int determining the number of different NT guides to have. It can be 0 in which case
#' no NT guides are added.
#'
#' @return A \code{data.frame} with columns \code{grna_id}, \code{grna_target}, \code{chr}, \code{start}, and
#' \code{end} as described above. The values of  \code{grna_id}, \code{grna_target}, and \code{chr} indicate
#' how they relate. For instance, a \code{grna_id} of \code{"g2_t1_c3_d2"} is the 2nd guide for the target
#' \code{"t1_c3_d2"}; that target is the first target on the chromosome \code{"c3_d2"}, which used
#' \code{chr_distances[2]} and \code{chr_starts[2]} (unless the latter is a scalar in which case it's just
#' \code{chr_starts}).
#'
#' @examples
#' # See `?make_mock_response_matrix_list` for a full example using all of these functions.
make_mock_grna_target_data <- function(num_guides_per_target, chr_distances, chr_starts, num_nt_guides) {
  if (length(chr_starts) == 1) {
    chr_starts <- rep(chr_starts, length(chr_distances))
  }
  if (length(chr_starts) != length(chr_distances)) {
    stop("`chr_starts` must be length 1 or have the same length as `chr_distances`.")
  }

  # first we make grna target data that is just doing different targets for different chromosomes
  # This data is formed by first having every target on a different chromosome,
  # and then repeating the guide counts but all on the same chromosome.
  target_and_chr_data <- .make_mock_target_and_chr_only(num_guides_per_target)

  # now we cross `target_and_chr_data` with each of the values in `chr_distance`
  # the result has the length of the final data and is just missing the grna_id column
  grna_target_data <- .make_mock_cross_target_chr_with_distances(
    target_and_chr_data = target_and_chr_data,
    chr_distances = chr_distances, chr_starts = chr_starts
  )

  # adding in unique names for each grna
  # this gets the unique target names, ensuring they stay in the order that they appear in the data
  target_labels <- with(grna_target_data, grna_target[!duplicated(grna_target)])
  # This vector enumerates each unique target: we get a sequence 1:(# times that target appears)
  grna_id_number <- table(grna_target_data$grna_target)[target_labels] |> sapply(function(i) 1:i) |> unlist()
  grna_target_data <- cbind(
    grna_id = paste0("g", grna_id_number, "_", grna_target_data$grna_target),
    grna_target_data
  )
  if (num_nt_guides > 0) {
    grna_target_data <- grna_target_data |>
      rbind(
        data.frame(
          grna_id = paste0("nt", 1:num_nt_guides),
          grna_target = "non-targeting",
          chr = NA_character_, start = NA, end = NA
        )
      )
  }
  return(grna_target_data)
}

# This function returns a matrix with the given dimension with patterns either across the rows or across the columns.
# The "patterns" are structured values that cover a wide range of cases. There are 13 patterns in all:
# constant: all 0, all 1, all `big`
# constant but for one element: 6 variations using values from c(0, 1, big)
# sequences: 4 variations, either increasing or decreasing
#
# To fill the specified number of dimensions, the patterns are repeated as many times as the whole set of 13 fits,
# and then if more dimensions needed the rest are filled in with iid Pois(1) noise.
# So for example, if we run `out <- .make_mock_patterned_matrix(12, 27, TRUE)`
# then we get a 12 x 27 matrix with the 1st column all 0, the second column all 1, and etc., and
# out[,1:13] is identical to out[,14:26] since the column patterns are repeated, and the final column out[,27] is
# Pois(1) noise.
.make_mock_patterned_matrix <- function(num_rows, num_cols, patterns_at_col_level = TRUE, big = 1000) {
  dims <- if (patterns_at_col_level) c(num_rows, num_cols) else c(num_cols, num_rows)
  num_patterns <- 13 # `patterned_matrix` will have 13 columns, one per pattern

  patterned_matrix <- cbind(
    0, 1, big, # constant
    rep(c(0,1), c(dims[1]-1, 1)),  rep(c(0, 1), c(1, dims[1]-1)), # all but one element constant
    rep(c(big, 0), c(dims[1]-1, 1)), rep(c(big, 1), c(dims[1]-1, 1)),
    rep(c(0, big), c(dims[1]-1, 1)), rep(c(1, big), c(dims[1]-1, 1)),
    0:(dims[1] - 1), dims[1]:1, # sequence columns
    seq(0, big, length = dims[1]) |> round(), seq(big, 1, length = dims[1]) |> round()
  ) |>
    magrittr::set_colnames(NULL) # the use of `big` adds a single column name

  # if we have more columns than the deterministic patterns, then
  # duplicate the patterned matrix as many times as we can and
  # fill the rest in with iid Pois(1) noice
  if (dims[2] > num_patterns) {
    num_copies <- dims[2] %/% num_patterns
    num_extra_cols <- dims[2] %% num_patterns
    patterned_matrix <- replicate(num_copies, patterned_matrix, simplify = FALSE) |>
      do.call(what = cbind)
    if (num_extra_cols > 0) {
      patterned_matrix <- cbind(
        patterned_matrix,
        matrix(rpois(num_extra_cols * dims[1], 1), ncol = num_extra_cols) # appending noise
      )
    }
  } else if (dims[2] < 13) {  # trim down to however many columns were asked for
    patterned_matrix <- patterned_matrix[,1:dims[2], drop = FALSE]
  }
  if (!patterns_at_col_level) {
    patterned_matrix <- t(patterned_matrix)
  }
  return(patterned_matrix)
}

#' Make a list of mock grna x cell expression matrices
#'
#' This function returns a list of matrices of grna expressions based on the provided \code{grna_target_data_frame}.
#'
#' These matrices are constructed by creating 5 response patterns for the targeting grnas, and separately
#' 5 response patterns for the non-targeting guides if any are present in \code{grna_target_data_frame}. The returned
#' list is all combinations of these: if there are no NT grnas provided, then the return is just a list of 5 with one
#' entry per response pattern; if instead there are NT grnas present in \code{grna_target_data_frame}, then the returned
#' list has all 25 combinations of each targeting grna response matrix on top of each NT grna response matrix. Note that
#' \code{grna_target_data_frame} will always have the non-NT on top of the NT since that is how
#' \code{make_mock_grna_target_data} structures it.
#'
#' @param grna_target_data_frame : the output of \code{make_mock_grna_target_data}
#' @param num_cells : int, the number of cells desired for the output
#' @param big : a large value to use in the expression matrices
#' @param return_as_sparse : return as a standard \code{matrix} or (if \code{TRUE}) as a \code{TsparseMatrix}.
#'
#' @return a list of matrices (sparse or not, as determined by \code{return_as_sparse}).
#'
#' @examples
#' # See `?make_mock_response_matrix_list` for a full example using all of these functions.
make_mock_grna_matrix_list <- function(grna_target_data_frame, num_cells, big = 10000, return_as_sparse = TRUE) {
  # first we make various data sets just for the targeting grnas
  non_nt_data <- dplyr::filter(grna_target_data_frame, grna_target != "non-targeting")
  num_targeting_guides <- nrow(non_nt_data)

  patterns <- list(
    all_zero = matrix(0, num_targeting_guides, num_cells),
    all_one = matrix(1, num_targeting_guides, num_cells),
    all_big = matrix(big, num_targeting_guides, num_cells),
    row_patterns = .make_mock_patterned_matrix(
      num_targeting_guides, num_cells, patterns_at_col_level = FALSE, big = big),
    col_patterns = .make_mock_patterned_matrix(
      num_targeting_guides, num_cells, patterns_at_col_level = TRUE, big = big)
  ) |>
    lapply(`dimnames<-`, list(non_nt_data$grna_id)) # , cell_names)) # not currently using cell names

  # we will do NT and non-NT separately and bind those together
  num_nt <- sum(grna_target_data_frame$grna_target == "non-targeting")
  if (num_nt > 0) {
    nt_names <- with(grna_target_data_frame, grna_id[grna_target == "non-targeting"])
    nt_patterns <- list(
      all_zero = matrix(0, num_nt, num_cells),
      all_one  = matrix(1, num_nt, num_cells),
      all_big  = matrix(big, num_nt, num_cells),
      # num_nt can realistically be 1, so just doing a pattern across the first row
      row_patterns = matrix(c(seq(big, 1, length = num_cells) |> round(), rep(1, num_nt * num_cells - num_cells)),
                            num_nt, num_cells, byrow = TRUE),
      col_patterns = .make_mock_patterned_matrix(num_nt, num_cells,
                                                 patterns_at_col_level = TRUE, big = big)
    ) |>
      lapply(`dimnames<-`, list(nt_names)) # , cell_names))  # not currently using cell names

    # combine by taking all combos
    new_patterns <- vector("list", length(patterns) * length(nt_patterns))
    k <- 1
    for (i in seq_along(patterns)) {
      for (j in seq_along(nt_patterns)) {
        new_patterns[[k]] <- rbind(patterns[[i]], nt_patterns[[j]])
        names(new_patterns)[k] <- paste0("non_nt_", names(patterns)[i], "_and_nt_", names(nt_patterns)[j])
        k <- k + 1
      }
    }
    patterns <- new_patterns
  }
  if (return_as_sparse) {
    return(lapply(patterns, as, "TsparseMatrix"))
  } else {
    return(patterns)
  }
}

#' Make a list of mock response expression matrices
#'
#' This function returns a list of 5 matrices of response expressions based on the provided dimensions.
#'
#' The matrices are: (1) all zero; (2) all 1; (3) all `big`;
#' (4) row patterns, as determined by \code{.make_mock_patterned_matrix}; and (5) column patterns, also as
#' determined by \code{.make_mock_patterned_matrix}.
#'
#' @param num_responses : int, the number of responses to make data for
#' @param num_cells : int, the number of cells to make data for
#' @param big : a large value to use in the expression matrices
#' @param return_as_sparse : return as a standard \code{matrix} or (if \code{TRUE}) as a \code{TsparseMatrix}.
#'
#' @return a list of matrices (sparse or not, as determined by \code{return_as_sparse}).
#'
#' @examples
#' set.seed(12321)
#' num_cells <- 50
#' mock_grna_target_data_frame <- make_mock_grna_target_data(
#'   num_guides_per_target = c(1, 5, 7, 2, 3), chr_distances = c(1, 3, 10),
#'   chr_starts = 1, num_nt_guides = 5
#' )
#' with(mock_grna_target_data_frame,
#'      cat("There are ", length(unique(grna_target)), " unique targets, ",
#'          length(unique(chr)), " unique chromosomes, and ", length(grna_id),
#'          " unique grna ids overall.", sep="")
#' )
#' mock_grna_matrix_list <- make_mock_grna_matrix_list(
#'   grna_target_data_frame = mock_grna_target_data_frame,
#'   num_cells = num_cells
#' )
#'
#' # looking at one particular example
#' mock_grna_matrix_list$non_nt_col_patterns_and_nt_row_patterns |> dim()
#' mock_grna_matrix_list$non_nt_col_patterns_and_nt_row_patterns[1:5,1:15]
#'
#' mock_response_matrix_list <- make_mock_response_matrix_list(
#'   num_responses = 30, num_cells = num_cells
#' )
#'
#' # looking at one particular example
#' mock_response_matrix_list$row_patterns |> dim()
#' mock_response_matrix_list$row_patterns[1:5,1:20]
#'
#' mock_extra_covariates  <- make_mock_extra_covariates_list(num_cells)
#' mock_extra_covariates$almost_constant_many_levels_two_vals |> dim()
#' mock_extra_covariates$almost_constant_many_levels_two_vals |> head()
#' mock_extra_covariates$almost_constant_many_levels_two_vals$batch |> table()
#'
make_mock_response_matrix_list <- function(num_responses, num_cells, big = 10000, return_as_sparse = TRUE) {
  response_names <- paste0("response_", 1:num_responses)
  patterns <- list(
    all_zero = matrix(0, num_responses, num_cells),
    all_one = matrix(1, num_responses, num_cells),
    all_big = matrix(big, num_responses, num_cells),
    row_patterns = .make_mock_patterned_matrix(
      num_responses, num_cells, patterns_at_col_level = FALSE, big = big),
    col_patterns = .make_mock_patterned_matrix(
      num_responses, num_cells, patterns_at_col_level = TRUE, big = big)
  ) |>
    lapply(`dimnames<-`, list(response_names)) #, cell_names))  # not currently using cell names

  if (return_as_sparse) {
    return(lapply(patterns, as, "TsparseMatrix"))
  } else {
    return(patterns)
  }
}

#' Make a list of mock extra covariates data.frames
#'
#' @param num_cells : the number of rows each data.frame will have.
#'
#' @return a list of 12 data.frames, each with \code{num_cells} rows. The
#' first 7 only have a single column, \code{batch}, and test constant, nearly constant,
#' many levels, and one missing factor level. The next 5 test adding additional
#' columns, which are either numeric, count, factor (two different options are considered),
#' or all of the above together.
#'
#' @examples
#' # See `?make_mock_response_matrix_list` for a full example using all of these functions.
#'
make_mock_extra_covariates_list <- function(num_cells) {
  # what replicate patterns do we want to have?
  if (num_cells < 20) {
    stop("`num_cells` must be >= 20 for this function.")
  }
  # this is mainly used for testing features other than batch
  three_level_batch <- data.frame(
    batch = rep(c("b1", "b2", "b3"), times = num_cells %/% 3 + c(0, 0, num_cells %% 3)) |>
      factor(levels = c("b1", "b2", "b3"))
  )
  patterns <- list(
    constant = data.frame(batch = rep("b1", num_cells) |> factor(levels = "b1")),
    almost_constant_two_levels_one_val = data.frame(batch = rep(c("b1", "b2"), c(1, num_cells - 1)) |>
                                                      factor(levels = c("b1", "b2"))),
    almost_constant_two_levels_two_vals = data.frame(batch = rep(c("b1", "b2"), c(2, num_cells - 2)) |>
                                                       factor(levels = c("b1", "b2"))),
    almost_constant_many_levels_one_val = data.frame(
      batch = rep(paste0("b", 1:10), c(rep(1, 9), num_cells - 9)) |>
        factor(levels = paste0("b", 1:10))
    ),
    almost_constant_many_levels_two_vals = data.frame(
      batch = rep(paste0("b", 1:10), c(rep(2, 9), num_cells - 18)) |>
        factor(levels = paste0("b", 1:10))
    ),
    # repeat b1:b10 as many times as they fit, padding the last few with b1
    many_levels = data.frame(
      batch = c(rep(paste0("b", 1:10), times = num_cells %/% 10 + (num_cells %% 10) * rep(0:1, c(9, 1)))) |>
        factor(levels = paste0("b", 1:10))
    ),
    # two levels appear but three are in the factor
    missing_level = data.frame(
      batch = rep(c("b1", "b2"), c(num_cells %/% 2, num_cells %/% 2 + num_cells %% 2)) |>
        factor(levels = c("b1", "b2", "b3"))
    ),
    numeric = three_level_batch |> dplyr::mutate(numeric = rnorm(num_cells)),
    count = three_level_batch |> dplyr::mutate(count = rpois(num_cells, 2)),
    factor_one_value_level = three_level_batch |> dplyr::mutate(
      factor = rep(c("fac1", "fac2"), c(1, num_cells - 1)) |>
        factor(levels = c("fac1", "fac2"))
    ),
    factor_many_levels = three_level_batch |> dplyr::mutate(
      factor = rep(paste0("fac", 1:5), num_cells %/% 5 + c(0, 0, 0, 0, num_cells %% 5)) |>
        factor(levels = paste0("fac", 1:5))
    ),
    many_columns = three_level_batch |> dplyr::mutate(
      numeric = rnorm(num_cells), count = rpois(num_cells, 2),
      factor = rep(paste0("fac", 1:5), num_cells %/% 5 + c(0, 0, 0, 0, num_cells %% 5)) |>
        factor(levels = paste0("fac", 1:5))
    )
  )
  return(patterns)
}
