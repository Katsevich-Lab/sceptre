#' Construct *cis* pairs
#'
#' `construct_cis_pairs()` is a helper function to facilitate construction the *cis* pairs. `construct_cis_pairs()` returns the set of target-response pairs for which the target and response are located on the same chromosome and in close physical proximity to one another. `construct_cis_pairs()` is a useful pair constructor function for screens that aim to map noncoding regulatory elements (e.g., enhancers or noncoding GWAS variants) to target genes in *cis*. `construct_cis_pairs()` assumes that the columns `chr`, `start`, and `stop` are present in the `grna_target_data_frame`, giving the chromosome, start position, and end position, respectively, of the region that each gRNA targets. `construct_cis_pairs()` takes several arguments: `sceptre_object` (required), `distance_threshold` (optional), `positive_control_pairs` (optional), and `response_position_data_frame` (optional). By default, `construct_cis_pairs()` pairs each gRNA target to the set of responses on the same chromosome as that target and within `distance_threshold` bases of that target. (The default value of `distance_threshold` is 500,000 bases, or half a megabase.) The `positive_control_pairs` data frame optionally can be passed to `construct_cis_pairs()`, in which case the positive control targets (i.e., the entries within the `grna_target` column of `positive_control_pairs`) are excluded from the *cis* pairs. One may want to exclude these from the discovery analysis if these targets are intended for positive control purposes only. See \href{https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#sec-set-analysis-parameters_construct_cis_pairs}{Section 2.2.2 of the manual} for more detailed information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#' @param distance_threshold (optional) target-response pairs located within `distance_threshold` bases of one another and on the same chromosome are included in the *cis* discovery set.
#' @param positive_control_pairs (optional) a data frame with columns `grna_target` and `response_id` containing the positive control pairs; if supplied, the positive control targets are excluded from the *cis* pairs.
#' @param response_position_data_frame (optional) a data frame with columns `response_id`, `chr`, and `position` giving the genomic coordinate of each response; by default `response_position_data_frame` is set to a data frame containing the genomic coordinate of each gene in the human genome relative to reference genome GRCh38.
#'
#' @return a data frame with columns `grna_target` and `response_id` containing the *cis* pairs
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
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' discovery_pairs <- construct_cis_pairs(sceptre_object,
#'   positive_control_pairs = positive_control_pairs,
#'   distance_threshold = 5e6
#' )
construct_cis_pairs <- function(sceptre_object, positive_control_pairs = data.frame(), distance_threshold = 500000L,
                                response_position_data_frame = gene_position_data_frame_grch38) {
  if (!all(c("response_id", "chr", "position") %in% colnames(response_position_data_frame))) {
    stop("`response_position_data_frame` must contain columns 'response_id', 'chr', and 'position'.")
  }
  grna_target_data_frame <- sceptre_object@grna_target_data_frame
  grna_target_data_frame <- data.table::as.data.table(grna_target_data_frame)
  response_ids <- rownames(get_response_matrix(sceptre_object))
  distance_threshold <- as.integer(distance_threshold)
  grna_targets_to_exclude <- c("non-targeting", as.character(positive_control_pairs$grna_target))

  # 1. subset grna group data frame so as to exclude non-targeting gRNAs and gRNA groups in grna_targets_to_exclude
  grna_target_data_frame <- grna_target_data_frame |> dplyr::filter(!(grna_target %in% grna_targets_to_exclude))

  # 2. subset response_position_data_frame so that only genes within response matrix remain
  response_position_data_frame <- response_position_data_frame |> dplyr::filter(response_id %in% response_ids)
  if (nrow(response_position_data_frame) == 0L) {
    stop("The response IDs (i.e., the rownames of `response_matrix`) must be a subset of the `response_id` column of `gene_position_data_frame_grch38`. You can view the first few rows of `gene_position_data_frame_grch38` via `head(gene_position_data_frame_grch38)`.")
  }

  # 3. loop over chromosomes in grna group df
  unique_chrs <- unique(grna_target_data_frame$chr)
  out_pairs <- lapply(X = unique_chrs, FUN = function(unique_chr) {
    response_position_data_frame_curr_chr <- response_position_data_frame[response_position_data_frame$chr == unique_chr, ]
    response_posits <- response_position_data_frame_curr_chr$position
    response_ids <- response_position_data_frame_curr_chr$response_id
    grna_target_data_frame_curr_chr <- grna_target_data_frame[grna_target_data_frame$chr == unique_chr, ]

    # 4. loop over unique grna groups in grna group df curr chr
    unique_grna_targets <- unique(grna_target_data_frame_curr_chr$grna_target)
    lapply(X = unique_grna_targets, FUN = function(unique_grna_target) {
      x <- grna_target_data_frame_curr_chr[grna_target_data_frame_curr_chr$grna_target == unique_grna_target, ]
      min_posit <- min(x$start)
      max_posit <- max(x$end)
      midpoint <- as.integer(floor((min_posit + max_posit) / 2))
      paired_responses <- response_ids[compute_genes_within_distance(midpoint, response_posits, distance_threshold)]
      if (length(paired_responses) >= 1L) {
        data.table::data.table(response_id = paired_responses, grna_target = unique_grna_target)
      } else {
        NULL
      }
    }) |> data.table::rbindlist()
  }) |>
    data.table::rbindlist() |>
    as.data.frame() |>
    dplyr::select(grna_target, response_id)

  return(out_pairs)
}


#' Construct *trans* pairs
#'
#' `construct_trans_pairs()` is a helper function to facilitate construction the set of *trans* pairs. `construct_trans_pairs()` returns the entire set of possible target-response pairs. `construct_trans_pairs()` is a useful pair constructor function for analyses in which we seek to conduct a *trans* analysis, testing each target against each response. `construct_trans_pairs()` takes as arguments `sceptre_object` (required), `positive_control_pairs` (optional), and `pairs_to_exclude` (optional). By default `construct_trans_pairs()` returns a data frame with columns `grna_target` and `response_id`, where each gRNA target is mapped to each response ID.
#'
#' The optional argument `pairs_to_exclude` enables the user to remove specific pairs from the *trans* set and takes values `"none"`, `"pc_pairs"`, or `"pairs_containing_pc_targets"`. If `pairs_to_exclude` is set to `"none"` (the default), then no pairs are removed from the *trans* set. Next, if `pairs_to_exclude` is set to `"pc_pairs"` (and the `positive_control_pairs` data frame is passed), then then the positive control target-response pairs are excluded from the *trans* set. Finally, if `pairs_to_exclude` is set to `"pairs_containing_pc_targets"` (and `positive_control_pairs` is passed), then *all* pairs containing a positive control gRNA target are excluded from the *trans* pairs. (In this sense setting `pairs_to_exclude` to `"pairs_containing_pc_targets"` is stronger than setting `pairs_to_exclude` to `"pc_pairs"`.) Typically, in gene-targeting (resp., noncoding-regulatory-element-targeting) screens, we set `pairs_to_exclude` to `"pc_pairs"` (resp., `"pairs_containing_pc_targets"`). See \href{https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#sec-set-analysis-parameters_construct_trans_pairs}{Section 2.2.2 of the manual} for more detailed information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#' @param positive_control_pairs (optional) the set of positive control pairs
#' @param pairs_to_exclude (optional; default `"none"`) a string specifying pairs to exclude from the *trans* pairs, one of `"none"`, `"pc_pairs"`, or `"pairs_containing_pc_targets"`
#'
#' @return a data frame with columns `grna_target` and `response_id` containing the *trans* discovery set
#' @export
#' @examples
#' # 1. low-moi, gene-targeting screen
#' data("lowmoi_example_data")
#' sceptre_object <- import_data(
#'   response_matrix = lowmoi_example_data$response_matrix,
#'   grna_matrix = lowmoi_example_data$grna_matrix,
#'   extra_covariates = lowmoi_example_data$extra_covariates,
#'   grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
#'   moi = "low"
#' )
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' discovery_pairs <- construct_trans_pairs(
#'   sceptre_object = sceptre_object,
#'   positive_control_pairs = positive_control_pairs,
#'   pairs_to_exclude = "pc_pairs"
#' )
#'
#' # 2. high-moi, enhancer-targeting screen
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
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' discovery_pairs <- construct_trans_pairs(
#'   sceptre_object = sceptre_object,
#'   positive_control_pairs = positive_control_pairs,
#'   pairs_to_exclude = "pairs_containing_pc_targets"
#' )
construct_trans_pairs <- function(sceptre_object, positive_control_pairs = data.frame(), pairs_to_exclude = "none") {
  if (!(pairs_to_exclude %in% c("none", "pc_pairs", "pairs_containing_pc_targets"))) {
    stop("`pairs_to_exclude` must be set to 'none', 'pc_pairs', or 'pairs_containing_pc_targets'.")
  }
  response_ids <- rownames(get_response_matrix(sceptre_object))
  grna_targets_to_exclude <- c("non-targeting", if (pairs_to_exclude == "pairs_containing_pc_targets") as.character(positive_control_pairs$grna_target) else NULL)
  grna_target_data_frame <- sceptre_object@grna_target_data_frame
  grna_targets <- grna_target_data_frame |>
    dplyr::filter(!(grna_target %in% grna_targets_to_exclude)) |>
    dplyr::pull(grna_target) |>
    unique()
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
#' `construct_positive_control_pairs()` is a helper function to facilitate construction of the positive control pairs. Positive control pairs are target-response pairs for which we know (or have strong reason to believe) that there is a regulatory relationship between the target and the response. We can use positive control pairs to verify that `sceptre` is sensitive (i.e., capable of detecting true associations) on the dataset under analysis. `construct_positive_control_pairs()` takes as an argument a `sceptre_object` and returns a data frame with columns `grna_target` and `response_id`, where gRNA targets and response IDs with matching names are paired. Typically, the positive control set consists of transcription start sites paired to the gene regulated by those transcription start sites. See \href{https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#positive-control-pairs}{Section 2.2 in the manual} for more detailed information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#'
#' @return a data frame with columns `grna_target` and `response_id` containing the positive control pairs
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
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
construct_positive_control_pairs <- function(sceptre_object) {
  grna_target_data_frame <- sceptre_object@grna_target_data_frame
  response_ids <- rownames(get_response_matrix(sceptre_object))
  pc_grna_targets <- grna_target_data_frame$grna_target[
    grna_target_data_frame$grna_target %in% response_ids
  ] |> unique()
  df <- data.frame(grna_target = pc_grna_targets, response_id = pc_grna_targets)
  return(df)
}
