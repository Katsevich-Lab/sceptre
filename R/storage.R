# Internal matrix boundary. ondisc stays optional: its namespace is only used
# for an actual backed operation. The selected route centralizes its existing
# private APIs here rather than requiring a companion package interface change.

#' Resolve a feature ID in source-row coordinates
#'
#' @param store Row-accessible sparse matrix or initialized `odm` with row names.
#' @param id One feature ID, already validated by the caller.
#' @return Integer R 1-based row positions; an absent ID gives `integer()`.
#'   Duplicate row names are not resolved or rejected here.
#' @noRd
storage_response_index <- function(store, id) {
    which(id == rownames(store))
}

#' Read one response in an ordered source-cell view
#'
#' @param store A `dgRMatrix` or initialized, process-local `odm`.
#' @param response_idx One valid R 1-based row index in `store`.
#' @param source_cells Valid R 1-based source-column indices, or `NULL` for all
#'   columns. Their order and repetitions are retained; `integer()` reads no cells.
#' @return Numeric vector in the requested cell order, including for one cell.
#' @details Materializes one complete source row before subsetting, never the full
#'   assay. Does not mutate inputs, draw randomness, or write backing files.
#' @noRd
storage_read_response <- function(store, response_idx, source_cells = NULL) {
    y <- if (methods::is(store, "odm")) {
        as.numeric(store[response_idx, ])
    } else {
        load_csr_row(store@j, store@p, store@x, response_idx, ncol(store))
    }
    if (is.null(source_cells)) y else y[source_cells]
}

#' Find stored source positions passing a feature threshold
#'
#' @param store A `dgRMatrix` or initialized `odm`.
#' @param feature_idx One valid R 1-based source-row index.
#' @param threshold Scalar count threshold.
#' @return R 1-based source-column indices in native stored-column order.
#' @details Only stored entries are tested: implicit zeros are absent even at
#'   threshold zero, while explicit stored zeros can pass. Sparse values and the
#'   threshold remain doubles; the existing backed representation and native
#'   threshold use integer coercion. Does not mutate inputs or use RNG/file writes.
#' @noRd
storage_threshold_row <- function(store, feature_idx, threshold) {
    if (methods::is(store, "odm")) {
        ondisc:::threshold_count_matrix_cpp(
            file_name_in = store@h5_file, f_row_ptr = store@ptr,
            row_idx = feature_idx, threshold = threshold
        )
    } else {
        threshold_count_matrix(store@j, store@p, store@x, feature_idx,
                               as.numeric(threshold))
    }
}

#' Own disposable copies of indices passed to external ondisc kernels
#'
#' @param groups List of index vectors; coordinates are unchanged here.
#' @return List of freshly allocated integer vectors, retaining names, order and
#'   repetitions. The copy protects callers from ondisc kernels that decrement
#'   indices in place; integer coercion alone can alias an existing vector.
#' @noRd
storage_owned_groups <- function(groups) {
    lapply(groups, function(group) as.integer(group) + 0L)
}

#' Count requested pairs and optional NT summaries in analysis coordinates
#'
#' @param store A `dgRMatrix` or initialized `odm` with responses in rows.
#' @param pair_index Data frame with R 1-based `response_idx` source rows and
#'   `group_idx` positions in `groups`. Request rows may repeat or be unordered.
#' @param groups List of treatment groups, each containing R 1-based positions in
#'   the ordered analysis view `source_cells`, not original source columns.
#' @param source_cells Valid R 1-based source-column indices defining that view.
#' @param controls Either `"complement"` or R 1-based analysis positions of NT cells.
#' @param nt_groups Optional list of NT groups in analysis coordinates. The
#'   low-MOI caller supplies its current NT groups; no union/deduplication occurs.
#'   `NULL` requests pair counts without retaining the all-response summary.
#' @param nuclear Logical; backed storage skips the summary and scans requested
#'   response rows only. Sparse storage retains its fused scan when summarizing.
#' @return List with integer `n_nonzero_trt` and `n_nonzero_cntrl` vectors in exact
#'   request order. A retained summary contains `n_nonzero_mat` (NT groups by all
#'   source responses) and `n_nonzero_tot` (one total per source response); otherwise
#'   they are the legacy `matrix()` and `integer()` placeholders.
#' @details Counts refer to stored column positions, including explicit zeros,
#'   rather than positive loaded values. Complement totals use the analysis view;
#'   NT totals sum the group counts, so overlapping groups would count repeatedly.
#'   The adapter sorts requests for the scan and restores order. Sparse kernels
#'   own coordinate conversion; only ondisc gets localized NT groups and disposable
#'   R copies for its mutating private interface. No RNG draws or file writes.
#' @noRd
storage_pair_counts <- function(
    store, pair_index, groups, source_cells, controls,
    nt_groups = NULL, nuclear = FALSE
) {
    complement <- is.character(controls) && identical(controls, "complement")
    all_nt_idxs <- if (complement) integer() else as.integer(controls)
    summarize <- !is.null(nt_groups)
    if (is.null(nt_groups)) {
        nt_groups <- if (complement) list() else list(all_nt_idxs)
    }
    scan_order <- order(pair_index$response_idx, method = "radix")
    response_idxs <- as.integer(pair_index$response_idx[scan_order])
    group_idxs <- as.integer(pair_index$group_idx[scan_order])
    if (methods::is(store, "odm")) {
        if (!complement) {
            nt_groups <- lapply(nt_groups, function(group) match(group, all_nt_idxs))
        }
        groups <- storage_owned_groups(groups)
        nt_groups <- storage_owned_groups(nt_groups)
        common <- list(
            file_name_in = store@h5_file, f_row_ptr = store@ptr,
            n_genes = nrow(store), n_cells_orig = ncol(store),
            n_cells_sub = length(source_cells), grna_group_idxs = groups,
            all_nt_idxs = all_nt_idxs, to_analyze_response_idxs = response_idxs,
            to_analyze_grna_idxs = group_idxs,
            control_group_complement = complement, cells_in_use = source_cells
        )
        if (nuclear || !summarize) {
            out <- do.call(ondisc:::compute_n_ok_pairs_ondisc, c(common, list(
                unique_response_idxs = sort(unique(response_idxs))
            )))
            out$n_nonzero_mat <- matrix()
            out$n_nonzero_tot <- integer()
        } else {
            out <- do.call(ondisc:::compute_nt_nonzero_matrix_and_n_ok_pairs_ondisc,
                           c(common, list(indiv_nt_grna_idxs = nt_groups)))
        }
    } else {
        out <- sparse_pair_counts(
            j = store@j, p = store@p, n_cells_orig = ncol(store),
            grna_group_idxs = groups, nt_groups_analysis = nt_groups,
            to_analyze_response_idxs = response_idxs,
            to_analyze_grna_idxs = group_idxs,
            control_group_complement = complement, cells_in_use = source_cells
        )
        if (!summarize) {
            out$n_nonzero_mat <- matrix()
            out$n_nonzero_tot <- integer()
        }
    }
    restore_order <- order(scan_order)
    out$n_nonzero_trt <- out$n_nonzero_trt[restore_order]
    out$n_nonzero_cntrl <- out$n_nonzero_cntrl[restore_order]
    out
}

#' Count analysis groups across requested source responses
#'
#' @inheritParams storage_pair_counts
#' @param response_idxs Valid R 1-based source-row indices in requested order.
#' @return Integer matrix with `length(groups)` rows and `length(response_idxs)`
#'   columns. Request order, repetitions, and zero/one-sized dimensions are kept.
#' @details The group kernel scans all responses, then this boundary
#'   selects columns. Counts use stored positions, including explicit zeros.
#'   Sparse kernels own base conversion; only ondisc requires disposable R copies.
#'   Inputs, RNG state, and backing files are otherwise unchanged.
#' @noRd
storage_group_counts <- function(store, response_idxs, groups, source_cells) {
    if (methods::is(store, "odm")) {
        groups <- storage_owned_groups(groups)
        out <- ondisc:::compute_n_trt_cells_matrix_ondisc(
            file_name_in = store@h5_file, f_row_ptr = store@ptr,
            n_cells_orig = ncol(store), n_cells_sub = length(source_cells),
            n_genes = nrow(store), nt_grna_group_idxs = groups,
            cells_in_use = source_cells
        )
    } else {
        out <- sparse_group_counts(
            j = store@j, p = store@p, n_cells_orig = ncol(store),
            groups = groups, cells_in_use = source_cells
        )
    }
    out[, response_idxs, drop = FALSE]
}

#' Create backing files with the existing shared import identity
#'
#' @param response_matrix,grna_matrix Sparse count matrices prepared by import.
#' @param directory_to_write Existing output directory prepared by the caller.
#' @param integer_id Scalar integer import ID, already drawn by the caller.
#' @return Named list `response_matrix`, `grna_matrix` of initialized `odm`s.
#' @details Writes `response.odm` then `grna.odm`; creation is not transactional.
#'   The private constructor preserves the supplied shared ID, unlike the public
#'   constructor's fixed ID. Existing ondisc count coercion, chunk-size restrictions
#'   and errors are retained. No RNG draws occur here; handles are process-local.
#' @noRd
storage_create_backed_pair <- function(
    response_matrix, grna_matrix, directory_to_write, integer_id
) {
    # Preserve the shared import identity and existing construction/RNG order.
    # The public ondisc constructor currently hardcodes ID zero.
    response_matrix <- ondisc:::create_odm_from_r_matrix_internal(
        mat = response_matrix,
        file_to_write = file.path(directory_to_write, "response.odm"),
        integer_id = integer_id
    )
    grna_matrix <- ondisc:::create_odm_from_r_matrix_internal(
        mat = grna_matrix,
        file_to_write = file.path(directory_to_write, "grna.odm"),
        integer_id = integer_id
    )
    list(response_matrix = response_matrix, grna_matrix = grna_matrix)
}

#' Reattach backing files and check their saved import identity
#'
#' @param response_file,grna_file Existing backing-file paths prepared by caller.
#' @param integer_id Scalar integer ID saved in the sceptre object.
#' @return Named list `response_matrix`, `grna_matrix` of initialized `odm`s.
#' @details Opens both files before checking their IDs against the object, retaining
#'   the existing mismatch errors. Reconstructs process-local handles through
#'   ondisc; reads files without writing them or drawing randomness.
#' @noRd
storage_reopen_backed_pair <- function(response_file, grna_file, integer_id) {
    response <- ondisc::initialize_odm_from_backing_file(response_file)
    grna <- ondisc::initialize_odm_from_backing_file(grna_file)
    if (integer_id != response@integer_id) {
        stop("The `sceptre_object` and `response_odm` have distinct IDs. The `sceptre_object` likely is associated with a different backing .odm file.")
    }
    if (integer_id != grna@integer_id) {
        stop("The `sceptre_object` and `grna_odm` have distinct IDs. The `sceptre_object` likely is associated with a different backing .odm file.")
    }
    list(response_matrix = response, grna_matrix = grna)
}

#' Compute existing cell summaries from column-accessible counts
#'
#' @param store A `dgCMatrix`, with features in rows and cells in columns.
#' @param mt_gene_idxs R 1-based mitochondrial feature-row indices; may be empty.
#' @param compute_p_mito,compute_max_feature Logical summary switches.
#' @return List of cell-ordered vectors: integer `n_nonzero` (loaded values > 0.5),
#'   numeric `n_umi`, `p_mito`, `frac_umis_max_feature`, and integer `max_feature`.
#'   The latter is deliberately a native zero-based row index, with first maximum
#'   winning ties. Disabled optional summaries retain native zero-filled vectors;
#'   zero-total ratios retain native non-finite results when computed.
#' @details The native kernel owns mitochondrial-index conversion. Does not
#'   convert storage, mutate caller inputs, or use RNG/files.
#' @noRd
storage_cell_covariates <- function(
    store, mt_gene_idxs, compute_p_mito, compute_max_feature
) {
    compute_cell_covariates_cpp(
        i = store@i, p = store@p, x = store@x,
        n_genes = nrow(store), n_cells = ncol(store),
        mt_gene_idxs = mt_gene_idxs,
        compute_p_mito = compute_p_mito,
        compute_max_feature = compute_max_feature
    )
}

#' Convert the existing supported matrix representations for row/column access
#'
#' @param matrix_in Numeric dense matrix, numeric or logical sparse R/C/T matrix,
#'   or initialized `odm`.
#' @param make_row_accessible Logical; request `dgRMatrix` rather than `dgCMatrix`.
#' @return Sparse matrix with original dimensions and dimnames, or the input `odm`.
#'   A numeric sparse input already in the requested orientation is returned as-is.
#' @details Sparse numeric values, including stored zeros and fractions, are kept;
#'   stored logical entries become one, including explicitly stored `FALSE`.
#'   Triplet duplicates are not normalized here. Sorting owns coordinate/value
#'   vectors and leaves caller inputs unchanged. This incurs sparse ownership
#'   copies and sorting workspace. Dense inputs are first coerced to general
#'   structure, retaining both triangles and near-symmetric entries; this may
#'   add one dense value copy. Sparse/backed inputs are not densified.
#'   No RNG draws or file writes occur. The local conversion
#'   keeps in-memory workflows independent of the optional ondisc dependency.
#' @noRd
set_matrix_accessibility <- function(matrix_in, make_row_accessible = TRUE) {
    owns_i <- owns_j <- owns_x <- FALSE
    # Preserve legacy presence semantics: every stored logical entry becomes 1.
    if (methods::is(matrix_in, "lgRMatrix")) {
        attr(matrix_in, "class") <- "dgRMatrix"
        matrix_in@x <- rep(1.0, length(matrix_in@j))
        owns_x <- TRUE
    }
    if (methods::is(matrix_in, "lgCMatrix")) {
        attr(matrix_in, "class") <- "dgCMatrix"
        matrix_in@x <- rep(1.0, length(matrix_in@i))
        owns_x <- TRUE
    }
    if (methods::is(matrix_in, "lgTMatrix")) {
        attr(matrix_in, "class") <- "dgTMatrix"
        matrix_in@x <- rep(1.0, length(matrix_in@j))
        owns_x <- TRUE
    }

    if (methods::is(matrix_in, "dgRMatrix") && make_row_accessible) {
        out <- matrix_in
    } else if (methods::is(matrix_in, "dgCMatrix") && !make_row_accessible) {
        out <- matrix_in
    } else if (methods::is(matrix_in, "odm")) {
        out <- matrix_in
    } else {
        n_responses <- nrow(matrix_in)
        response_ids <- rownames(matrix_in)
        cell_ids <- colnames(matrix_in)
        n_cells <- ncol(matrix_in)
        if (methods::is(matrix_in, "dgTMatrix")) {
            i <- matrix_in@i
            j <- matrix_in@j
        } else if (methods::is(matrix_in, "matrix")) {
            # Auto-detecting symmetric structure would store only one triangle;
            # a general output must retain every original dense entry instead.
            matrix_in <- methods::as(
                methods::as(matrix_in, "generalMatrix"), "TsparseMatrix"
            )
            i <- matrix_in@i
            j <- matrix_in@j
            owns_i <- owns_j <- owns_x <- TRUE
        } else if (methods::is(matrix_in, "dgCMatrix")) {
            i <- matrix_in@i
            j <- sparse_pointer_to_indices(matrix_in@p)
            owns_j <- TRUE
        } else if (methods::is(matrix_in, "dgRMatrix")) {
            i <- sparse_pointer_to_indices(matrix_in@p)
            j <- matrix_in@j
            owns_i <- TRUE
        } else {
            stop("Class not recognized.")
        }
        x <- matrix_in@x

        # Sort only owned columns. Sparse slots may alias caller vectors; dense
        # coercion, pointer expansion and logical replacement create fresh data.
        # Explicit duplication preserves value bits and avoids recopying these
        # fresh columns. Ownership follows constructors, never NAMED/refcounts.
        dt <- list(i = if (owns_i) i else data.table::copy(i),
                   j = if (owns_j) j else data.table::copy(j),
                   x = if (owns_x) x else data.table::copy(x))
        data.table::setDT(dt)

        if (make_row_accessible) {
            data.table::setorderv(x = dt, cols = c("i", "j"))
            p <- obtain_pointer_vector(i = dt$i, dim = n_responses)
            out <- Matrix::sparseMatrix(j = 1, p = c(0, 1), x = 1, repr = "R")
            out@j <- dt$j
        } else {
            data.table::setorderv(x = dt, cols = c("j", "i"))
            p <- obtain_pointer_vector(i = dt$j, dim = n_cells)
            out <- Matrix::sparseMatrix(i = 1, p = c(0, 1), x = 1, repr = "C")
            out@i <- dt$i
        }
        out@p <- p
        out@x <- dt$x
        out@Dim <- c(n_responses, n_cells)
        rownames(out) <- response_ids
        colnames(out) <- cell_ids
    }
    return(out)
}


#' Expand a compressed sparse pointer into native coordinates
#'
#' @param p Valid zero-based compressed pointer, length dimension plus one, with
#'   at least one dimension entry; differences are nonnegative stored-entry counts.
#' @return Zero-based coordinate vector of length `tail(p, 1)`, with empty runs
#'   omitted. Does not mutate `p` or change entry order.
#' @noRd
sparse_pointer_to_indices <- function(p) {
    l <- diff(p)
    rep(x = seq(0, length(l) - 1L), times = l)
}

