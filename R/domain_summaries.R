# Canonical guide metadata has one owner. Preparation stores only a row-aligned
# derived group label, preserving repeated metadata rows and their original order.
guide_group_metadata <- function(object) {
    metadata <- object@data@grna_target_data_frame
    groups <- object@preparation@guide_groups
    if (!length(groups)) groups <- metadata$grna_target
    if (length(groups) != nrow(metadata)) stop("Guide grouping does not align with canonical metadata.")
    metadata$grna_group <- groups
    metadata
}

# Totals count requested target-response rows, not compiled guide rows. Passing
# Bonferroni hypotheses count target-response families with at least one eligible
# guide; singleton/union retain the existing passing-row definition.
pair_count <- function(object, type = "discovery", passing = FALSE) {
    stopifnot(type %in% c("discovery", "positive_control"))
    if (!passing) return(as.integer(nrow(methods::slot(object@spec, paste0(type, "_pairs")))))
    pairs <- methods::slot(object@preparation, paste0(type, "_pairs_with_info"))
    if (!nrow(pairs) || !"pass_qc" %in% names(pairs)) return(0L)
    if (object@spec@grna_integration_strategy == "bonferroni") {
        metadata <- object@data@grna_target_data_frame[c("grna_id", "grna_target")]
        pairs <- dplyr::left_join(pairs, metadata, by = c("grna_group" = "grna_id"))
        eligible <- dplyr::summarize(dplyr::group_by(pairs, response_id, grna_target),
            any_pass_qc = any(pass_qc), .groups = "drop")
        return(as.integer(sum(eligible$any_pass_qc)))
    }
    as.integer(sum(pairs$pass_qc))
}

mean_assignment_size <- function(object) mean(lengths(object@preparation@guide_assignments))
