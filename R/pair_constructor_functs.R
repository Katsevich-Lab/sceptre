#' Construct cis pairs
#'
#' @param sceptre_object TBD
#' @param distance_threshold TBD
#' @param grna_groups_to_exclude TBD
#' @param response_grna_group_pairs_to_exclude TBD
#' @param ref_genome TBD
#'
#' @return TBD
#' @export
construct_cis_pairs <- function(sceptre_object, positive_control_pairs = data.frame(), distance_threshold = 500000L,
                                exclude_pc_grna_groups = TRUE, ref_genome = "10X_GRCh38_2020") {
  if (ref_genome != "10X_GRCh38_2020") {
    stop("The only reference genome currently available is the GRCh38 (2020) reference genome provided by 10X.")
  }
  grna_group_data_frame <- data.table::as.data.table(sceptre_object@grna_group_data_frame)
  gene_ids <- rownames(sceptre_object@response_matrix)
  distance_threshold <- as.integer(distance_threshold)
  grna_groups_to_exclude <- c("non-targeting", if (exclude_pc_grna_groups) positive_control_pairs$grna_group else NULL)

  # 1. subset grna group data frame so as to exclude non-targeting gRNAs and gRNA groups in grna_groups_to_exclude
  grna_group_data_frame <- grna_group_data_frame |> dplyr::filter(!(grna_group %in% grna_groups_to_exclude))

  # 2. subset gene_table so that only genes within response matrix remain
  gene_table <- gene_table |> dplyr::filter(gene_id %in% gene_ids)

  # 3. loop over chromosomes in grna group df
  unique_chrs <- unique(grna_group_data_frame$chr)
  out_pairs <- lapply(X = unique_chrs, FUN = function(unique_chr) {
    gene_table_curr_chr <- gene_table[gene_table$chr == unique_chr,]
    gene_tss_posits <- gene_table_curr_chr$tss_position
    gene_ids <- factor(gene_table_curr_chr$gene_id)
    grna_group_data_frame_curr_chr <- grna_group_data_frame[grna_group_data_frame$chr == unique_chr,]

    # 4. loop over unique grna groups in grna group df curr chr
    unique_grna_groups <- factor(unique(grna_group_data_frame_curr_chr$grna_group))
    lapply(X = unique_grna_groups, FUN = function(unique_grna_group) {
      x <- grna_group_data_frame_curr_chr[grna_group_data_frame_curr_chr$grna_group == unique_grna_group,]
      min_posit <- min(x$start); max_posit <- max(x$end)
      midpoint <- as.integer(floor((min_posit + max_posit)/2))
      paired_genes <- gene_ids[compute_genes_within_distance(midpoint, gene_tss_posits, distance_threshold)]
      if (length(paired_genes) >= 1L) {
        data.table::data.table(response_id = paired_genes, grna_group = unique_grna_group)
      } else {
        NULL
      }
    }) |> data.table::rbindlist()
  }) |> data.table::rbindlist()

  return(out_pairs)
}


#' Construct trans pairs
#'
#' A helper function to construct the trans pairs
#'
#' @param sceptre_object TBD
#' @param grna_groups_to_exclude TBD
#' @param response_grna_group_pairs_to_exclude TBD
#'
#' @return TBD
#' @export
construct_trans_pairs <- function(sceptre_object, positive_control_pairs = data.frame(), exclude_positive_control_pairs = TRUE) {
  response_ids <- rownames(sceptre_object@response_matrix) |> factor()
  grna_groups <- sceptre_object@grna_group_data_frame |>
    dplyr::filter(!(grna_group %in% "non-targeting")) |>
    dplyr::pull(grna_group) |> unique() |> factor()
  out_pairs <- expand.grid(response_id = response_ids, grna_group = grna_groups)
  if (exclude_positive_control_pairs) out_pairs <- exclude_pairs(out_pairs, positive_control_pairs)
  return(out_pairs)
}


#' Construct positive control pairs
#'
#' Helper function to construct the positive control pairs.
#'
#' @param sceptre_object TBD
#' @param grna_groups_to_exclude TBD
#' @param response_grna_group_pairs_to_exclude TBD
#'
#' @return TBD
#' @export
construct_positive_control_pairs <- function(sceptre_object) {
  grna_group_data_frame <- sceptre_object@grna_group_data_frame
  response_ids <- rownames(sceptre_object@response_matrix)
  pc_grna_groups <- grna_group_data_frame$grna_group[
    grna_group_data_frame$grna_group %in% response_ids] |> unique()
  df <- data.frame(grna_group = pc_grna_groups, response_id = pc_grna_groups)
  return(df)
}


exclude_pairs <- function(pairs_df, pairs_to_exclude_df) {
  if (nrow(pairs_to_exclude_df) >= 1) {
    pairs_to_exclude_str <- paste0(pairs_to_exclude_df$response_id,
                                   pairs_to_exclude_df$grna_group)
    out <- pairs_df |> dplyr::mutate(pair_str = paste0(response_id, grna_group)) |>
      dplyr::filter(!(pair_str %in% pairs_to_exclude_str)) |>
      dplyr::select(-pair_str)
  } else {
    out <- pairs_df
  }
  return(out)
}
