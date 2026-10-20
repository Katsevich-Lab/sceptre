# A cell can have three different position numbers:
# * source: its column in the original count matrix;
# * analysis: its position among the cells kept after QC;
# * fit: its position among the cells selected, in order, for one model fit.
# For example, if QC keeps original columns c(2, 4, 5), original cell 5 is
# analysis cell 3. If a fit uses original cells c(5, 2), it is fit cell 1.
# The helpers below translate between these lists. All R positions start at 1.

#' Check the original matrix columns selected for analysis
#' @param source_cells Original count-matrix column numbers, in the order to
#'   use.
#'   They must be whole numbers, without missing values or duplicates.
#' @param n_source Number of columns in the original count matrix.
#' @return The same column numbers as an integer vector. This helper accepts an
#'   empty vector; the QC workflow separately rejects an analysis with no cells.
#' @noRd
prepare_cell_map <- function(source_cells, n_source) {
    if (
        !is.numeric(source_cells) ||
            anyNA(source_cells) ||
            any(source_cells != as.integer(source_cells)) ||
            any(source_cells < 1L | source_cells > n_source) ||
            anyDuplicated(source_cells)
    ) {
        stop(
            "Analysis source cells must be unique integer positions within ",
            "the response store."
        )
    }
    as.integer(source_cells)
}

#' Find a group's cells in the list of cells kept after QC
#' @param source_members Group members given as original count-matrix columns.
#' @param analysis_to_source Original matrix columns of the cells kept after QC,
#'   in analysis order. For example, `c(2, 4, 5)` means analysis cell 1 is
#'   original column 2, analysis cell 2 is original column 4, and so on.
#' @return Positions in `analysis_to_source`, in that list's order, with each
#'   position included once. For example, members `c(5, 2, 5)` give `c(1, 3)`
#'   for the list above. Members removed by QC are omitted; no matches gives
#'   `integer()`. The original order of `source_members` does not affect this.
#' @noRd
source_membership_to_analysis <- function(source_members, analysis_to_source) {
    which(analysis_to_source %in% source_members)
}

#' Store nontargeting guide groups as positions within one NT-cell list
#' @param nt_groups_analysis Named list of nontargeting (NT) guide groups. Each
#'   group contains positions among the cells kept after QC. Low-MOI groups
#'   must not overlap.
#' @return A list with `all_nt_idxs`, the combined NT cells in first-seen order,
#'   and `indiv_nt_grna_idxs`, each group's positions within that combined list.
#'   For example, groups `A = c(4, 1)` and `B = 5` give the combined list
#'   `c(4, 1, 5)` and within-list positions `A = c(1, 2)` and `B = 3`.
#'   Empty input gives an empty cell vector and group list.
#' @details Matching cell numbers to the combined list also handles overlapping
#'   groups without inventing positions beyond the end of that list.
#' @noRd
encode_nt_fit_assignments <- function(nt_groups_analysis) {
    nt_to_analysis <- as.integer(unique(unlist(
        nt_groups_analysis,
        use.names = FALSE
    )))
    groups_fit <- lapply(nt_groups_analysis, match, table = nt_to_analysis)
    list(indiv_nt_grna_idxs = groups_fit, all_nt_idxs = nt_to_analysis)
}

#' Collect cells from the NT guides named in a calibration group
#' @param group_id NT guide names joined by `&`, for example `"NT1&NT2"`.
#' @param nt_groups Named list giving each NT guide's cells. These are usually
#'   positions among cells kept after QC; this function keeps those numbers.
#' @param low_moi If TRUE, concatenate the guides' cell lists. If FALSE, remove
#'   repeated cells, keeping the order in which they first appear.
#' @return One unnamed vector of cell positions. For `"NT1&NT2"`, NT1's cells
#'   come before any additional NT2 cells. The position numbers are not
#'   converted
#'   to original matrix columns or to positions within a model fit.
#' @noRd
calibration_group_membership <- function(group_id, nt_groups, low_moi) {
    nt_ids <- strsplit(group_id, "&", fixed = TRUE)[[1L]]
    members <- stats::setNames(unlist(nt_groups[nt_ids]), NULL)
    if (!low_moi) {
        members <- unique(members)
    }
    members
}

#' Match each requested gene-guide test to its gene and guide-group positions
#' @param pairs Requested tests, with `response_id` and `grna_group` columns.
#' @param response_ids,group_ids Lists of unique gene and guide-group names to
#'   match against. The caller has already checked that the requested names
#'   exist.
#' @return A table with the original request row (`pair_row_id`), gene position
#'   (`response_idx`) and guide-group position (`group_idx`). All start at 1.
#'   Request order and repeated requests are kept. This helper only matches
#'   names;
#'   checks for missing names or forbidden duplicates belong to the caller.
#' @noRd
resolve_pair_rows <- function(pairs, response_ids, group_ids) {
    data.frame(
        pair_row_id = seq_len(nrow(pairs)),
        response_idx = match(pairs$response_id, response_ids),
        group_idx = match(pairs$grna_group, group_ids)
    )
}

#' Prepare the cell and pair lookups used by an analysis
#' @param cells_in_use,n_source Original matrix columns kept after QC, in order,
#'   and the original number of matrix columns.
#' @param grna_assignments Saved guide assignments after QC. Targeting groups
#'   already use positions among the kept cells. When using NT controls, saved
#'   NT groups instead use positions within the separate NT-cell list.
#' @param pairs,response_ids Requested tests and the original matrix's gene
#'   names.
#' @param calibration_check,control_group_complement,low_moi Logical settings
#'   identifying calibration, complement controls and low-MOI data.
#' @return A list connecting original matrix columns to cells kept after QC,
#'   including the ordered NT-cell list and guide memberships. It also matches
#'   each requested pair to a gene and guide group and lists request rows for
#'   each gene and group, preserving request order within each list.
#' @details Saved NT-group positions are translated back to positions among all
#'   kept cells here, once. Subsequent helpers can then use the same numbering.
#'   This prepares lookups; it does not fit models or change the saved object.
#' @noRd
prepare_index_context <- function(
    cells_in_use,
    n_source,
    grna_assignments,
    pairs,
    response_ids,
    calibration_check,
    control_group_complement,
    low_moi
) {
    analysis_to_source <- prepare_cell_map(cells_in_use, n_source)
    nt_to_analysis <- if (control_group_complement) {
        integer()
    } else {
        grna_assignments$all_nt_idxs
    }
    nt_groups_analysis <- grna_assignments$indiv_nt_grna_idxs
    if (!control_group_complement) {
        nt_groups_analysis <- lapply(nt_groups_analysis, function(nt_fit_pos) {
            nt_to_analysis[nt_fit_pos]
        })
    }
    group_ids <- if (calibration_check) {
        unique(as.character(pairs$grna_group))
    } else {
        names(grna_assignments$grna_group_idxs)
    }
    pair_rows <- resolve_pair_rows(pairs, response_ids, group_ids)
    list(
        analysis_to_source = analysis_to_source,
        nt_to_analysis = nt_to_analysis,
        groups_analysis = grna_assignments$grna_group_idxs,
        nt_groups_analysis = nt_groups_analysis,
        pair_rows = pair_rows,
        response_ids = response_ids,
        group_ids = group_ids,
        rows_by_response = split(pair_rows$pair_row_id, pair_rows$response_idx),
        rows_by_group = split(pair_rows$pair_row_id, pair_rows$group_idx),
        calibration_check = calibration_check,
        control_group_complement = control_group_complement,
        low_moi = low_moi
    )
}

#' Get the treated cells for one guide group
#' @param index Cell and pair lookups returned by `prepare_index_context()`.
#' @param group_id A targeting guide group, or joined NT guide names for
#'   calibration.
#' @return Positions among the cells kept after QC. Calibration combines the NT
#'   guides using `calibration_group_membership()`; other analyses use the saved
#'   targeting group directly. These are not original matrix column numbers.
#' @noRd
group_treatment_analysis <- function(index, group_id) {
    if (index$calibration_check) {
        calibration_group_membership(
            group_id,
            index$nt_groups_analysis,
            index$low_moi
        )
    } else {
        index$groups_analysis[[group_id]]
    }
}

#' Choose the cells and their order for a model fit
#' @param index Cell and pair lookups returned by `prepare_index_context()`.
#' @param group_id Guide-group name, or NULL to prepare cells shared across
#'   groups.
#' @return A list containing `fit_to_analysis` (positions among cells kept after
#'   QC), `source_cells` (the corresponding original matrix columns),
#'   `trt_fit_pos` (treated cells' positions within this fit) and `n_trt`.
#' @details With complement controls, use all kept cells. With NT controls,
#'   calibration uses the ordered NT-cell list, while discovery uses the chosen
#'   treatment cells followed by the NT cells. For NT discovery, NULL returns
#'   all
#'   kept cells so each group worker can later select its own subset.
#'
#'   Example: QC keeps original columns `c(2, 4, 5, 8, 11)`. Treatment occupies
#'   positions `c(2, 3)` in that list and NT cells occupy `c(4, 1, 5)`.
#'   The discovery fit therefore reads original columns `c(4, 5, 8, 2, 11)`;
#'   treatment occupies positions `1:2` within this fit. Counts and
#'   design-matrix
#'   rows must follow this same order. Membership order and repeats are kept.
#'
#'   Each selection is built when needed, rather than keeping a cell list for
#'   every guide at once. This helper only selects cells; it does not read
#'   counts
#'   or fit a model.
#' @noRd
make_fit_view <- function(index, group_id = NULL) {
    treatment <- if (is.null(group_id)) {
        integer()
    } else {
        group_treatment_analysis(index, group_id)
    }
    if (index$control_group_complement) {
        fit_to_analysis <- seq_along(index$analysis_to_source)
        trt_fit_pos <- treatment
    } else if (index$calibration_check) {
        fit_to_analysis <- index$nt_to_analysis
        trt_fit_pos <- match(treatment, fit_to_analysis)
    } else if (is.null(group_id)) {
        fit_to_analysis <- seq_along(index$analysis_to_source)
        trt_fit_pos <- integer()
    } else {
        fit_to_analysis <- c(treatment, index$nt_to_analysis)
        trt_fit_pos <- seq_along(treatment)
    }
    view <- new(
        "CellView",
        fit_to_analysis = as.integer(fit_to_analysis),
        source_cells = as.integer(index$analysis_to_source[fit_to_analysis]),
        trt_fit_pos = as.integer(trt_fit_pos),
        n_trt = as.integer(length(trt_fit_pos))
    )
    validObject(view)
    view
}

#' Make a function that finds treated cells in the shared model-fitting data
#' @param index Cell and pair lookups returned by `prepare_index_context()`.
#' @return A function taking a guide-group name and returning `trt_idxs`, the
#'   treated cells' positions starting at 1, and `n_trt`, their count.
#' @details For NT calibration, positions refer to the ordered NT-cell list.
#'   Otherwise they refer to all cells kept after QC. NT discovery workers use
#'   their own treatment-plus-NT selections instead. The returned function keeps
#'   one set of lookups, rather than storing cell lists or design matrices for
#'   every requested pair.
#' @noRd
make_treatment_index_lookup <- function(index) {
    force(index)
    function(group_id) {
        trt_analysis_pos <- group_treatment_analysis(index, group_id)
        trt_fit_pos <- if (
            index$calibration_check && !index$control_group_complement
        ) {
            match(trt_analysis_pos, index$nt_to_analysis)
        } else {
            trt_analysis_pos
        }
        list(trt_idxs = trt_fit_pos, n_trt = length(trt_fit_pos))
    }
}

#' Collect the data, cell selections and settings needed to run the tests
#' @param object Sceptre object containing the original design matrix and any
#'   response fits that are still valid for its current settings.
#' @param analysis One of `calibration_check`, `power_check`,
#'   `discovery_analysis`.
#' @param pairs Validated gene-guide pairs that passed QC, in request order.
#' @param plan Optional validated transient run plan.
#' @return A temporary list containing the count matrix (`response_store`), cell
#'   and pair lookups (`index`), shared cell selection (`common_view`), and
#'   design
#'   rows selected in the same order (`design`). It also contains the requested
#'   pairs, analysis name, test settings and reusable response fits.
#' @details `shared_response_fit` says whether groups use the same response fit.
#'   NT discovery instead fits each group's treatment-plus-NT cells separately.
#'   `test_options` includes B1, B2, B3, their checked sum, test direction and
#'   the
#'   tail-approximation setting. Keeping `drop = FALSE` preserves a matrix even
#'   when the formula contains only an intercept. The saved object and its full
#'   design matrix are left unchanged.
#' @noRd
prepare_analysis_context <- function(object, analysis, pairs, plan = NULL) {
    plan <- validate_analysis_context_plan(object, analysis, pairs, plan)
    analysis <- sub("^run_", "", analysis)
    total_resamples <- checked_resampling_total(plan@B1, plan@B2, plan@B3)
    store <- get_response_matrix(object)
    index <- prepare_index_context(
        object@preparation@analysis_cells,
        ncol(store),
        object@preparation@analysis_groups,
        pairs,
        rownames(store),
        identical(analysis, "calibration_check"),
        object@spec@control_group_complement,
        object@data@low_moi
    )
    common_view <- make_fit_view(index)
    identity <- response_fit_identity(object, common_view@source_cells)
    perturbation_identity <- perturbation_cache_identity(object)
    cache <- if (identical(object@fits@response_identity, identity)) {
        object@fits@response_fits
    } else {
        list()
    }
    list(
        fit_identity = identity,
        response_store = store,
        index = index,
        common_view = common_view,
        design = object@preparation@design_matrix[
            common_view@source_cells,
            ,
            drop = FALSE
        ],
        pairs = pairs,
        analysis = analysis,
        shared_response_fit = index$calibration_check ||
            index$control_group_complement,
        test_options = analysis_test_options(object, plan, total_resamples),
        response_fits = cache,
        perturbation_identity = perturbation_identity,
        perturbation_fits = analysis_perturbation_cache(
            object,
            index,
            perturbation_identity
        )
    )
}

# Keep numerical settings and cache eligibility distinct from cell mapping.
analysis_test_options <- function(object, plan, total_resamples) {
    list(
        B1 = plan@B1,
        B2 = plan@B2,
        B3 = plan@B3,
        total_resamples = total_resamples,
        side_code = object@spec@side_code,
        fit_parametric_curve = object@spec@resampling_approximation ==
            "skew_normal"
    )
}

analysis_perturbation_cache <- function(object, index, perturbation_identity) {
    if (
        !index$calibration_check &&
            identical(
                object@fits@perturbation_identity,
                perturbation_identity
            )
    ) {
        object@fits@perturbation_fits
    } else {
        list()
    }
}

validate_analysis_context_plan <- function(object, analysis, pairs, plan) {
    if (is.null(plan)) {
        plan <- make_run_plan(object, analysis, pairs)
    }
    plan <- validate_run_plan(object, plan)
    if (!identical(plan@analysis, analysis_record_name(analysis))) {
        stop("Analysis type does not match run plan.")
    }
    expected <- dplyr::filter(plan@pairs, pass_qc)
    if (!identical(lapply(pairs, identity), lapply(expected, identity))) {
        stop("Analysis pairs do not match the run plan's eligible family.")
    }
    plan
}
