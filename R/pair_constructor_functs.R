#' Return cis pairs
#'
#' @param sceptre_object
#' @param grna_groups_to_exclude
#' @param distance_threshold
#'
#' @return
#' @export
#'
#' @examples
return_cis_pairs <- function(sceptre_object, grna_groups_to_exclude = NULL, distance_threshold = 500000, ref_genome = "10X_GRCh38") {
  if (ref_genome != "10X_GRCh38") {
    stop("The only reference genome currently available is the GRCh38 reference genome provided by 10X.")
  }
  grna_group_data_frame <- sceptre_object@grna_group_data_frame
  gene_ids <- rownames(sceptre_object@response_matrix)
  gene_table |> head()
  # etc...
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
