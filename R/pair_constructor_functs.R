#' Construct cis pairs
#'
#' @param sceptre_object TBD
#' @param positive_control_pairs TBD
#' @param distance_threshold TBD
#' @param exclude_positive_control_grna_targets TBD
#' @param ref_genome TBD
#'
#' @return TBD
#' @export
construct_cis_pairs <- function(sceptre_object, positive_control_pairs = NULL,
                                distance_threshold = 500000L, ref_genome = "10X_GRCh38_2020") {
  if (ref_genome != "10X_GRCh38_2020") {
    stop("The only reference genome currently available is the GRCh38 (2020) reference genome provided by 10X.")
  }
  grna_target_data_frame <- data.table::as.data.table(sceptre_object@grna_target_data_frame)
  gene_ids <- rownames(sceptre_object@response_matrix)
  distance_threshold <- as.integer(distance_threshold)
  grna_targets_to_exclude <- c("non-targeting", if (!is.null(positive_control_pairs)) as.character(positive_control_pairs$grna_target) else NULL)

  # 1. subset grna group data frame so as to exclude non-targeting gRNAs and gRNA groups in grna_targets_to_exclude
  grna_target_data_frame <- grna_target_data_frame |> dplyr::filter(!(grna_target %in% grna_targets_to_exclude))

  # 2. subset gene_table so that only genes within response matrix remain
  gene_table <- gene_table |> dplyr::filter(gene_id %in% gene_ids)

  # 3. loop over chromosomes in grna group df
  unique_chrs <- unique(grna_target_data_frame$chr)
  out_pairs <- lapply(X = unique_chrs, FUN = function(unique_chr) {
    gene_table_curr_chr <- gene_table[gene_table$chr == unique_chr,]
    gene_tss_posits <- gene_table_curr_chr$tss_position
    gene_ids <- gene_table_curr_chr$gene_id
    grna_target_data_frame_curr_chr <- grna_target_data_frame[grna_target_data_frame$chr == unique_chr,]

    # 4. loop over unique grna groups in grna group df curr chr
    unique_grna_targets <- unique(grna_target_data_frame_curr_chr$grna_target)
    lapply(X = unique_grna_targets, FUN = function(unique_grna_target) {
      x <- grna_target_data_frame_curr_chr[grna_target_data_frame_curr_chr$grna_target == unique_grna_target,]
      min_posit <- min(x$start); max_posit <- max(x$end)
      midpoint <- as.integer(floor((min_posit + max_posit)/2))
      paired_genes <- gene_ids[compute_genes_within_distance(midpoint, gene_tss_posits, distance_threshold)]
      if (length(paired_genes) >= 1L) {
        data.table::data.table(response_id = paired_genes, grna_target = unique_grna_target)
      } else {
        NULL
      }
    }) |> data.table::rbindlist()
  }) |> data.table::rbindlist() |> as.data.frame() |> dplyr::select(grna_target, response_id)

  return(out_pairs)
}


#' Construct trans pairs
#'
#' A helper function to construct the trans pairs
#'
#' @param sceptre_object TBD
#' @param positive_control_pairs TBD
#' @param exclude_positive_control_pairs TBD
#'
#' @return TBD
#' @export
construct_trans_pairs <- function(sceptre_object, positive_control_pairs = data.frame(), exclude_positive_control_pairs = TRUE) {
  response_ids <- rownames(sceptre_object@response_matrix)
  grna_targets <- sceptre_object@grna_target_data_frame |>
    dplyr::filter(!(grna_target %in% "non-targeting")) |>
    dplyr::pull(grna_target) |> unique()
  out_pairs <- expand.grid(response_id = response_ids, grna_target = grna_targets)
  if (exclude_positive_control_pairs) out_pairs <- exclude_pairs(out_pairs, positive_control_pairs)
  return(out_pairs)
}


#' Construct positive control pairs
#'
#' Helper function to construct the positive control pairs.
#'
#' @param sceptre_object TBD
#'
#' @return TBD
#' @export
construct_positive_control_pairs <- function(sceptre_object) {
  grna_target_data_frame <- sceptre_object@grna_target_data_frame
  response_ids <- rownames(sceptre_object@response_matrix)
  pc_grna_targets <- grna_target_data_frame$grna_target[
    grna_target_data_frame$grna_target %in% response_ids] |> unique()
  df <- data.frame(grna_target = pc_grna_targets, response_id = pc_grna_targets)
  return(df)
}


exclude_pairs <- function(pairs_df, pairs_to_exclude_df) {
  if (nrow(pairs_to_exclude_df) >= 1) {
    pairs_to_exclude_str <- paste0(pairs_to_exclude_df$response_id,
                                   pairs_to_exclude_df$grna_target)
    out <- pairs_df |> dplyr::mutate(pair_str = paste0(response_id, grna_target)) |>
      dplyr::filter(!(pair_str %in% pairs_to_exclude_str)) |>
      dplyr::select(-pair_str)
  } else {
    out <- pairs_df
  }
  return(out)
}
