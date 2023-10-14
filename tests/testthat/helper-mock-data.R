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
  attr(patterned_matrix, "dimnames") <- NULL
  return(patterned_matrix)
}

# this function takes a character vector `patterns_to_make`, which is a subset of
# `c("zero", "one", "big", "row", "column")` (note that the calling function enforces this),
# and returns a list of matrices with the indicated patterns and with `pattern_nrow` and
# `pattern_ncol` rows and columns, respectively. `big` is passed to `.make_mock_patterned_matrix`
# and the final rownames of each matrix are set to be `pattern_rownames`.
.make_mock_matrix_pattern_list <- function(patterns_to_make, pattern_nrow, pattern_ncol, big, pattern_rownames) {
  # ensuring they appear in this order if `patterns_to_make` doesn't list them in this order
  patterns_to_make <- intersect(c("zero", "one", "big", "row", "column"), patterns_to_make)
  pattern_list <- list()
  if ("zero" %in% patterns_to_make) { # matrix with all 0
    pattern_list <- c(pattern_list, matrix(0, pattern_nrow, pattern_ncol) |> list())
  }
  if ("one" %in% patterns_to_make ) { # matrix with all 1
    pattern_list <- c(pattern_list, matrix(1, pattern_nrow, pattern_ncol) |> list())
  }
  if ("big" %in% patterns_to_make ) { # matrix with all `big`
    pattern_list <- c(pattern_list, matrix(big, pattern_nrow, pattern_ncol) |> list())
  }
  if ("row" %in% patterns_to_make ) { # matrix with patterns in the rows
    pattern_list <- c(
      pattern_list,
      .make_mock_patterned_matrix(
        pattern_nrow, pattern_ncol, patterns_at_col_level = FALSE, big = big
      ) |> list()
    )
  }
  if ("column" %in% patterns_to_make ) { # matrix with patterns in the columns
    pattern_list <- c(
      pattern_list,
      .make_mock_patterned_matrix(
        pattern_nrow, pattern_ncol, patterns_at_col_level = TRUE, big = big
      ) |> list()
    )
  }

  pattern_list <- pattern_list |> lapply(`dimnames<-`, list(pattern_rownames))
  names(pattern_list) <- patterns_to_make
  return(pattern_list)
}

#' Makes either a single mock grna x cell expression matrix or a list of them,
#' based on the provided \code{grna_target_data_frame}
#'
#' These matrices are constructed using 5 possible response patterns for the targeting grnas, and separately
#' 5 possible response patterns for the non-targeting guides if any are present in \code{grna_target_data_frame}.
#' The returned matrices are all combinations of these: if there are no NT grnas provided, then the return is
#' just the chosen patterns for the non-NT guide RNAs; if instead there are NT grnas present in
#' \code{grna_target_data_frame}, then there are up to 25 matrices that can be returned, corresponding to
#' all possible combinations of the 5 non-NT patterns and the 5 NT patterns.
#'
#' Note that \code{grna_target_data_frame} will always have the non-NT on top of the NT since that is how
#' \code{make_mock_grna_target_data} structures it, so if both non-NT and NT guides are present, the
#' returned matrices will always have the non-NT rows on top.
#'
#' @param grna_target_data_frame : the output of \code{make_mock_grna_target_data}
#' @param num_cells : int, the number of cells desired for the output
#' @param non_nt_patterns : a character vector that is a subset of c("zero", "one", "big", "row", "column")
#' or is equal to "all". This determines which non-NT patterns are made. If "all" then all 5 patterns are made.
#' @param nt_patterns : a character vector that is a subset of c("zero", "one", "big", "row", "column")
#' or is equal to "all". This determines which NT patterns are made. If "all" then all 5 patterns are made.
#' Ignored if there are no NT guides present in \code{grna_target_data_frame}.
#' @param big : int, a large value to use in the expression matrices
#' @param return_as_sparse : logical, return as a standard \code{matrix} or (if \code{TRUE}) as a \code{TsparseMatrix}.
#' @param unlist_if_single_pattern : logical, if TRUE, then if only a single pattern is made (which will happen precisely
#' when \code{non_nt_patterns} has a single non-"all" value, and either no NTs are present or \code{nt_patterns}
#' also specifies just a single non-"all" value) then the actual matrix is returned rather than a list of
#' length 1 containing that pattern.
#'
#' @return either a matrix or a list of matrices (sparse or not, as determined by \code{return_as_sparse}),
#' depending on \code{unlist_if_single_pattern} and the number of patterns specified.
#'
#' @examples
#' # See `?make_mock_response_matrix_list` for a full example using all of these functions.
make_mock_grna_matrices <- function(grna_target_data_frame, num_cells, non_nt_patterns, nt_patterns = NULL,
                                    big = 10000, return_as_sparse = TRUE, unlist_if_single_pattern = TRUE) {
  allowed_patterns <- c("all", "zero", "one", "big", "row", "column")
  if (!all(non_nt_patterns %in% allowed_patterns) || length(non_nt_patterns) == 0) {
    stop(paste0('`non_nt_patterns` must be a non-empty subset of c("',
                paste0(allowed_patterns[-1], collapse = '", "'), '") or be equal to "',
                allowed_patterns[1], '".'))
  }
  if ("all" %in% non_nt_patterns) {
    non_nt_patterns <- allowed_patterns[-1] # taking every actual pattern
  }

  # non-NT and NT are handled separately, and combined at the end if any NT are present.
  num_nt <- sum(grna_target_data_frame$grna_target == "non-targeting")
  num_non_nt <- nrow(grna_target_data_frame) - num_nt
  non_nt_ids <- with(grna_target_data_frame, grna_id[grna_target != "non-targeting"])

  if (num_nt > 0) {
    if (is.null(nt_patterns)) {
      stop("`grna_target_data_frame` contains non-targeting gRNA but `nt_patterns` has not been set.")
    }
    if (!all(nt_patterns %in% allowed_patterns) || length(non_nt_patterns) == 0) {
      stop(paste0('`nt_patterns` must be a non-empty subset of c("',
                  paste0(allowed_patterns[-1], collapse = '", "'), '") or be equal to "',
                  allowed_patterns[1], '".'))
    }
    if ("all" %in% nt_patterns) {
      nt_patterns <- allowed_patterns[-1] # taking every actual pattern
    }
  }

  non_nt_pattern_list <- .make_mock_matrix_pattern_list(
    patterns_to_make = non_nt_patterns,
    pattern_nrow = num_non_nt,
    pattern_ncol = num_cells,
    big = big,
    pattern_rownames = non_nt_ids
  )
  names(non_nt_pattern_list) <- paste0("non_nt_", names(non_nt_pattern_list))


  if (num_nt == 0) {
    final_pattern_list <- non_nt_pattern_list  # no NTs means these are all the patterns there are
  } else {
    # if there are NTs, we need to take each combination of non-NT pattern with each
    # combination of NT pattern
    nt_ids <- with(grna_target_data_frame, grna_id[grna_target == "non-targeting"])
    nt_pattern_list <- .make_mock_matrix_pattern_list(
      patterns_to_make = nt_patterns,
      pattern_nrow = num_nt,
      pattern_ncol = num_cells,
      big = big,
      pattern_rownames = nt_ids
    )
    names(nt_pattern_list) <- paste0("nt_", names(nt_pattern_list))

    # loop through this list filling in each element with one element of
    # `non_nt_pattern_list` stacked on top of another element of `nt_pattern_list`
    final_pattern_list <- vector("list", length(non_nt_pattern_list) * length(nt_pattern_list))
    k <- 1
    for (i in seq_along(non_nt_pattern_list)) {
      for (j in seq_along(nt_pattern_list)) {
        final_pattern_list[[k]] <- rbind(non_nt_pattern_list[[i]], nt_pattern_list[[j]])
        names(final_pattern_list)[k] <- paste0(names(non_nt_pattern_list)[i], "_", names(nt_pattern_list)[j])
        k <- k + 1
      }
    }
  }
  if (return_as_sparse) {
    final_pattern_list <- lapply(final_pattern_list, as, "TsparseMatrix")
  }
  if (unlist_if_single_pattern && length(final_pattern_list) == 1) {
    return(final_pattern_list[[1]])
  } else {
    return(final_pattern_list)
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
#' @param patterns : a character vector that is a subset of c("zero", "one", "big", "row", "column")
#' or is equal to "all". This determines which patterns are made. If "all" then all 5 patterns are made.
#' @param big : int, a large value to use in the expression matrices
#' @param return_as_sparse : logical, return as a standard \code{matrix} or (if \code{TRUE}) as a \code{TsparseMatrix}.
#' @param unlist_if_single_pattern logical, if TRUE, then if only a single pattern is made (which will happen precisely
#' when \code{patterns} has a single non-"all" value) then the actual matrix is returned rather than a list of
#' length 1 containing that pattern.
#'
#' @return either a single matrix or a list of matrices: whether or not it is a list is determined by the
#' length of \code{patterns} and the value of \code{unlist_if_single_pattern}; whether or not it is sparse is
#' determined by \code{return_as_sparse}.
#'
#' @examples
#' set.seed(12321)
#' num_cells <- 27
#' num_responses <- 21
#' mock_grna_target_data_frame <- make_mock_grna_target_data(
#'   num_guides_per_target = c(1, 5, 7, 2, 3), chr_distances = c(1, 3, 10),
#'   chr_starts = 1, num_nt_guides = 5
#' )
#' with(mock_grna_target_data_frame,
#'      cat("There are ", length(unique(grna_target)), " unique targets, ",
#'          length(unique(chr)), " unique chromosomes, and ", length(grna_id),
#'          " unique grna ids overall.", sep="")
#' )
#' mock_grna_matrix_list <- make_mock_grna_matrices(
#'   grna_target_data_frame = mock_grna_target_data_frame,
#'   num_cells = num_cells,
#'   non_nt_patterns = "all",
#'   nt_patterns = "all"
#' )
#'
#' # looking at one particular example
#' mock_grna_matrix_list$non_nt_column_nt_row |> dim()
#' mock_grna_matrix_list$non_nt_column_nt_row[1:5,1:15]
#' # if we just wanted this one matrix, we could have gotten it this way
#' # (note that due to RNG some of the values in the full matrix aren't the
#' # same, but the pattern is the same and this subset is the same):
#' make_mock_grna_matrices(
#'   grna_target_data_frame = mock_grna_target_data_frame,
#'   num_cells = num_cells,
#'   non_nt_patterns = "column",
#'   nt_patterns = "row"
#' )[1:5, 1:15]
#'
#' mock_response_matrix_list <- make_mock_response_matrices(
#'   num_responses = num_responses,
#'   num_cells = num_cells,
#'   patterns = "all"
#' )
#'
#' # looking at one particular example
#' mock_response_matrix_list$row |> dim()
#' mock_response_matrix_list$row[1:5,1:20]
#'
#' # if we only wanted that one matrix, we could have gotten it this way:
#' make_mock_response_matrices(
#'   num_responses = num_responses,
#'   num_cells = num_cells,
#'   patterns = "row"
#' )[1:5, 1:15]
#'
#'
#' mock_extra_covariates  <- make_mock_extra_covariates_list(num_cells)
#' mock_extra_covariates$almost_constant_many_levels_two_vals |> dim()
#' mock_extra_covariates$almost_constant_many_levels_two_vals |> head()
#' mock_extra_covariates$almost_constant_many_levels_two_vals$batch |> table()
make_mock_response_matrices <- function(num_responses, num_cells, patterns,
                                        big = 10000, return_as_sparse = TRUE,
                                        unlist_if_single_pattern = TRUE) {

  allowed_patterns <- c("all", "zero", "one", "big", "row", "column")
  if (!all(patterns %in% allowed_patterns) || length(patterns) == 0) {
    stop(paste0('`patterns` must be a non-empty subset of c("',
                paste0(allowed_patterns[-1], collapse = '", "'), '") or be equal to "',
                allowed_patterns[1], '".'))
  }
  if ("all" %in% patterns) {
    patterns <- allowed_patterns[-1] # taking every actual pattern
  }

  pattern_list <- .make_mock_matrix_pattern_list(
    patterns_to_make = patterns,
    pattern_nrow = num_responses,
    pattern_ncol = num_cells,
    big = big,
    pattern_rownames = paste0("response_", 1:num_responses)
  )

  if (return_as_sparse) {
    pattern_list <- lapply(pattern_list, as, "TsparseMatrix")
  }
  if (unlist_if_single_pattern && length(pattern_list) == 1) {
    return(pattern_list[[1]])
  } else {
    return(pattern_list)
  }
}

#' Makes either a single data.frame or a list of data.frames, mocking extra_covariates
#'
#' @param num_cells : int, the number of rows per data.frame. Must be >= 20.
#' @param patterns : a character vector indicating which of the 12 patterns to include. See below for
#' details. If \code{patterns = "all"} then all 12 patterns are returned.
#' @param unlist_if_single_pattern : boolean, if TRUE then if only a single pattern is selected
#' then just that data.frame is returned, rather than a list containing that data.frame. Ignored if
#' \code{length(patterns) >= 2}.
#'
#' There are 12 different patterns possible.
#'
#' The first 7 patterns are all a data.frame with exactly one column, a factor named "batch".
#' (1) "constant" : "batch" is constant
#' (2) "almost_constant_two_levels_one_val" : "batch" has two levels. The first level appears
#' exactly one time and the rest of the values are the second level.
#' (3) "almost_constant_two_levels_two_vals" : "batch" has two levels. The first level appears
#' exactly two times and the rest of the values are the second level.
#' (4) "almost_constant_many_levels_one_val" : "batch" has 10 levels. The first 9 levels appears
#' exactly one time and the rest of the values are the 10th level.
#' (5) "almost_constant_many_levels_two_vals" :"batch" has 10 levels. The first 9 levels appears
#' exactly two times and the rest of the values are the 10th level.
#' (6) "many_levels" : "batch" has 10 levels, with each appearing an equal number of times
#' (or as close to equal as possible).
#' (7) "missing_level" :  "batch" has 3 levels, but only values for the first two levels
#' appear in the data
#'
#' The remaining 5 patterns all have at least two columns. In every case the first column is a factor named "batch"
#' with 3 levels which is as close to balanced as possible.
#'
#' (8) "numeric" : a data.frame with two columns: "batch" and "numeric", where the latter is \code{rnorm(num_cells)}
#' (9) "count" : a data.frame with two columns: "batch" and "count", where the latter is \code{rpois(num_cells, 2)}
#' (10) "factor_one_value_level" : a data.frame with two columns: "batch" and "factor", where "factor" has two levels,
#' with the first level appearing just once and the remaining values all being the second level
#' (11) "factor_many_levels" : a data.frame with two columns: "batch" and "factor", where "factor" has 5 levels and
#' is as close to balanced as possible.
#' (12) "many_columns" : a data.frame with 4 columns: "batch", "numeric" (same as in pattern 8),
#' "count" (same as in pattern 9), and "factor" (same as in pattern 11).
#'
#' @return either a data.frame or a list of data.frames (depending on \code{length(patterns)} and
#' \code{unlist_if_single_pattern}), each with \code{num_cells} many rows.
#'
#' @examples
#' See \code{make_mock_response_matrices} for a full example with all functions.
make_mock_extra_covariates_data_frames <- function(num_cells, patterns, unlist_if_single_pattern = TRUE) {
  # we need at least 20 obs to be able to do all the desired factor patterns.
  if (num_cells < 20) {
    stop("`num_cells` must be >= 20 for this function.")
  }
  allowed_patterns <- c(
    "all", "constant", "almost_constant_two_levels_one_val",
    "almost_constant_two_levels_two_vals", "almost_constant_many_levels_one_val",
    "almost_constant_many_levels_two_vals", "many_levels", "missing_level",
    "numeric", "count", "factor_one_value_level", "factor_many_levels", "many_columns"
  )

  if (!all(patterns %in% allowed_patterns) || length(patterns) == 0) {
    stop(paste0('`patterns` must be a non-empty subset of c("',
                paste0(allowed_patterns[-1], collapse = '", "'), '") or be equal to "',
                allowed_patterns[1], '".'))
  }
  if ("all" %in% patterns) {
    patterns <- allowed_patterns[-1] # taking every actual pattern
  }
  patterns <- intersect(allowed_patterns[-1], patterns)  # ensuring ordering is as expected

  # the idea of this is to have a well-behaved "batch" column to which other more
  # weird columns are added.
  three_level_batch <- data.frame(
    batch = rep(c("b1", "b2", "b3"), times = num_cells %/% 3 + c(0, 0, num_cells %% 3)) |>
      factor(levels = c("b1", "b2", "b3"))
  )

  pattern_list <- list()
  if ("constant" %in% patterns) {
    pattern_list <- c(
      pattern_list,
      data.frame(batch = rep("b1", num_cells) |> factor(levels = "b1")) |> list()
    )
  }
  if ("almost_constant_two_levels_one_val" %in% patterns) {
    pattern_list <- c(
      pattern_list,
      data.frame(batch = rep(c("b1", "b2"), c(1, num_cells - 1)) |>
                   factor(levels = c("b1", "b2"))) |> list()
    )
  }
  if ("almost_constant_two_levels_two_vals" %in% patterns) {
    pattern_list <- c(
      pattern_list,
      data.frame(batch = rep(c("b1", "b2"), c(2, num_cells - 2)) |>
                   factor(levels = c("b1", "b2"))) |> list()
    )
  }
  if ("almost_constant_many_levels_one_val" %in% patterns) {
    pattern_list <- c(
      pattern_list,
      data.frame(
        batch = rep(paste0("b", 1:10), c(rep(1, 9), num_cells - 9)) |>
          factor(levels = paste0("b", 1:10))
      ) |> list()
    )
  }
  if ("almost_constant_many_levels_two_vals" %in% patterns) {
    pattern_list <- c(
      pattern_list,
      data.frame(
        batch = rep(paste0("b", 1:10), c(rep(2, 9), num_cells - 18)) |>
          factor(levels = paste0("b", 1:10))
      ) |> list()
    )
  }
  if ("many_levels" %in% patterns) {
    # repeat b1:b10 as many times as they fit
    pattern_list <- c(
      pattern_list,
      data.frame(
        batch = c(rep(paste0("b", 1:10), times = num_cells %/% 10 + (num_cells %% 10) * rep(0:1, c(9, 1)))) |>
          factor(levels = paste0("b", 1:10))
      ) |> list()
    )
  }
  if ("missing_level" %in% patterns) {
    # two levels appear but three are in the factor
    pattern_list <- c(
      pattern_list,
      data.frame(
        batch = rep(c("b1", "b2"), c(num_cells %/% 2, num_cells %/% 2 + num_cells %% 2)) |>
          factor(levels = c("b1", "b2", "b3"))
      ) |> list()
    )
  }
  if ("numeric" %in% patterns) {
    pattern_list <- c(
      pattern_list,
      three_level_batch |> dplyr::mutate(numeric = rnorm(num_cells)) |> list()
    )
  }
  if ("count" %in% patterns) {
    pattern_list <- c(
      pattern_list,
      three_level_batch |> dplyr::mutate(count = rpois(num_cells, 2)) |> list()
    )
  }
  if ("factor_one_value_level" %in% patterns) {
    pattern_list <- c(
      pattern_list,
      three_level_batch |> dplyr::mutate(
        factor = rep(c("fac1", "fac2"), c(1, num_cells - 1)) |>
          factor(levels = c("fac1", "fac2"))
      ) |> list()
    )
  }
  if ("factor_many_levels" %in% patterns) {
    pattern_list <- c(
      pattern_list,
      three_level_batch |> dplyr::mutate(
        factor = rep(paste0("fac", 1:5), num_cells %/% 5 + c(0, 0, 0, 0, num_cells %% 5)) |>
          factor(levels = paste0("fac", 1:5))
      ) |> list()
    )
  }
  if ("many_columns" %in% patterns) {
    pattern_list <- c(
      pattern_list,
      three_level_batch |> dplyr::mutate(
        numeric = rnorm(num_cells), count = rpois(num_cells, 2),
        factor = rep(paste0("fac", 1:5), num_cells %/% 5 + c(0, 0, 0, 0, num_cells %% 5)) |>
          factor(levels = paste0("fac", 1:5))
      ) |> list()
    )
  }
  names(pattern_list) <- patterns

  if (unlist_if_single_pattern && length(pattern_list) == 1) {
    return(pattern_list[[1]])
  } else {
    return(pattern_list)
  }
}
