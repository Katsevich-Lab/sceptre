#' Return cis pairs
#'
#' @param sceptre_object
#' @param grna_groups_to_exclude
#' @param distance_threshold
#'
#' @return
#' @export
return_cis_pairs <- function(sceptre_object, grna_groups_to_exclude = character(), distance_threshold = 500000L, ref_genome = "10X_GRCh38_2020") {
  if (ref_genome != "10X_GRCh38_2020") {
    stop("The only reference genome currently available is the GRCh38 (2020) reference genome provided by 10X.")
  }
  grna_group_data_frame <- data.table::as.data.table(sceptre_object@grna_group_data_frame)
  gene_ids <- rownames(sceptre_object@response_matrix)
  distance_threshold <- as.integer(distance_threshold)

  # 1. subset grna group data frame so as to exclude non-targeting gRNAs and gRNA groups in grna_groups_to_exclude
  grna_group_data_frame <- grna_group_data_frame |>
    dplyr::filter(!(grna_group %in% c(grna_groups_to_exclude, "non-targeting")))

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


#' Return all pairs
#'
#' @param sceptre_object
#' @param grna_groups_to_exclude
#'
#' @return
#' @export
return_all_pairs <- function(sceptre_object, grna_groups_to_exclude = NULL) {
  response_ids <- rownames(sceptre_object@response_matrix) |> factor()
  grna_groups <- sceptre_object@grna_group_data_frame |>
    dplyr::filter(grna_group != "non-targeting") |>
    dplyr::pull(grna_group) |> unique() |> factor()
  expand.grid(response_id = response_ids, grna_group = grna_groups)
}
