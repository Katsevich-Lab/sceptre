#' Construct *cis* pairs
#'
#' `construct_cis_pairs()` is a helper function to facilitate construction the *cis* pairs. See \href{https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#sec-set-analysis-parameters_construct_cis_pairs}{Section 2.2.2} of the manual for more detailed information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#' @param distance_threshold (optional) target-response pairs located within `distance_threshold` bases of one another and on the same chromosome are included in the *cis* discovery set.
#' @param positive_control_pairs (optional) a data frame with columns `grna_target` and `response_id` containing the positive control pairs; if supplied, the positive control targets are excluded from the *cis* pairs.
#' @param response_position_data_frame (optional) a data frame with columns `response_id`, `chr`, and `position` giving the genomic coordinate of each response; by default `response_position_data_frame` is set to a data frame containing the genomic coordinate of each gene in the human genome relative to reference genome GRCh38.
#'
#' @return a data frame with columns `grna_target` and `response_id` containing the *cis* discovery set
#' @export
construct_cis_pairs <- function(sceptre_object, positive_control_pairs = data.frame(), distance_threshold = 500000L,
                                response_position_data_frame = gene_position_data_frame_grch38) {
  if (!all(colnames(response_position_data_frame) %in% c("response_id", "chr", "position"))) {
    stop("`response_position_data_frame` must contain columns 'response_id', 'chr', and 'position'.")
  }
  if (nrow(sceptre_object@grna_target_data_frame_with_vector) >= 1L) {
    grna_target_data_frame <- sceptre_object@grna_target_data_frame_with_vector
  } else {
    grna_target_data_frame <- sceptre_object@grna_target_data_frame
  }
  grna_target_data_frame <- data.table::as.data.table(grna_target_data_frame)
  response_ids <- rownames(sceptre_object@response_matrix)
  distance_threshold <- as.integer(distance_threshold)
  grna_targets_to_exclude <- c("non-targeting", as.character(positive_control_pairs$grna_target))

  # 1. subset grna group data frame so as to exclude non-targeting gRNAs and gRNA groups in grna_targets_to_exclude
  grna_target_data_frame <- grna_target_data_frame |> dplyr::filter(!(grna_target %in% grna_targets_to_exclude))

  # 2. subset response_position_data_frame so that only genes within response matrix remain
  response_position_data_frame <- response_position_data_frame |> dplyr::filter(response_id %in% response_ids)

  # 3. loop over chromosomes in grna group df
  unique_chrs <- unique(grna_target_data_frame$chr)
  out_pairs <- lapply(X = unique_chrs, FUN = function(unique_chr) {
    response_position_data_frame_curr_chr <- response_position_data_frame[response_position_data_frame$chr == unique_chr,]
    response_posits <- response_position_data_frame_curr_chr$position
    response_ids <- response_position_data_frame_curr_chr$response_id
    grna_target_data_frame_curr_chr <- grna_target_data_frame[grna_target_data_frame$chr == unique_chr,]

    # 4. loop over unique grna groups in grna group df curr chr
    unique_grna_targets <- unique(grna_target_data_frame_curr_chr$grna_target)
    lapply(X = unique_grna_targets, FUN = function(unique_grna_target) {
      x <- grna_target_data_frame_curr_chr[grna_target_data_frame_curr_chr$grna_target == unique_grna_target,]
      min_posit <- min(x$start); max_posit <- max(x$end)
      midpoint <- as.integer(floor((min_posit + max_posit)/2))
      paired_responses <- response_ids[compute_genes_within_distance(midpoint, response_posits, distance_threshold)]
      if (length(paired_responses) >= 1L) {
        data.table::data.table(response_id = paired_responses, grna_target = unique_grna_target)
      } else {
        NULL
      }
    }) |> data.table::rbindlist()
  }) |> data.table::rbindlist() |> as.data.frame() |> dplyr::select(grna_target, response_id)

  return(out_pairs)
}


#' Construct *trans* pairs
#'
#' `construct_trans_pairs()` is a helper function to facilitate construction the set of *trans* pairs. See \href{https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#sec-set-analysis-parameters_construct_trans_pairs}{Section 2.2.2} of the manual for more detailed information about this function.
#'
#' Typically, in screens of genes (resp., noncoding regulatory elements), we set `pairs_to_exclude` to `"pc_pairs"` (resp., `"pairs_containing_pc_targets"`).
#'
#' @param sceptre_object a `sceptre_object`
#' @param positive_control_pairs (optional) the set of positive control pairs
#' @param pairs_to_exclude (optional; default `"none"`) a string specifying pairs to exclude from the *trans* pairs, one of `"none"`, `"pc_pairs"`, or `"pairs_containing_pc_targets"`
#'
#' @return a data frame with columns `grna_target` and `response_id` containing the *trans* discovery set
#' @export
construct_trans_pairs <- function(sceptre_object, positive_control_pairs = data.frame(), pairs_to_exclude = "none") {
  if (!(pairs_to_exclude %in% c("none", "pc_pairs", "pairs_containing_pc_targets"))) {
    stop("`pairs_to_exclude` must be set to 'none', 'pc_pairs', or 'pairs_containing_pc_targets'.")
  }
  response_ids <- rownames(sceptre_object@response_matrix)
  grna_targets_to_exclude <- c("non-targeting", if (pairs_to_exclude == "pairs_containing_pc_targets") as.character(positive_control_pairs$grna_target) else NULL)
  if (nrow(sceptre_object@grna_target_data_frame_with_vector) >= 1L) {
    grna_target_data_frame <- sceptre_object@grna_target_data_frame_with_vector
  } else {
    grna_target_data_frame <- sceptre_object@grna_target_data_frame
  }
  grna_targets <- grna_target_data_frame |>
    dplyr::filter(!(grna_target %in% grna_targets_to_exclude)) |>
    dplyr::pull(grna_target) |> unique()
  out_pairs <- expand.grid(grna_target = grna_targets, response_id = response_ids, stringsAsFactors = FALSE)
  data.table::setDT(out_pairs)
  if (pairs_to_exclude == "pc_pairs") {
    data.table::setDT(positive_control_pairs)
    if (nrow(positive_control_pairs) >= 1) {
      out_pairs <- data.table::fsetdiff(out_pairs, positive_control_pairs)
    }
  }
  return(out_pairs)
}


#' Construct positive control pairs
#'
#' `construct_positive_control_pairs()` is a helper function to facilitate construction of the positive control pairs. See \href{https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#positive-control-pairs}{Section 2.2 in the manual} for more detailed information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#'
#' @return a data frame with columns `grna_target` and `response_id` containing the positive control pairs
#' @export
construct_positive_control_pairs <- function(sceptre_object) {
  if (nrow(sceptre_object@grna_target_data_frame_with_vector) >= 1L) {
    grna_target_data_frame <- sceptre_object@grna_target_data_frame_with_vector
  } else {
    grna_target_data_frame <- sceptre_object@grna_target_data_frame
  }
  response_ids <- rownames(sceptre_object@response_matrix)
  pc_grna_targets <- grna_target_data_frame$grna_target[
    grna_target_data_frame$grna_target %in% response_ids] |> unique()
  df <- data.frame(grna_target = pc_grna_targets, response_id = pc_grna_targets)
  return(df)
}
