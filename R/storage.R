# Read and convert count matrices here, whether counts are in R memory or in
# ondisc files. Private ondisc calls are kept in this file so the rest of the
# analysis need not know how those files work. ondisc is loaded only when used.

#' Find a gene's row in the count matrix
#' @param store Sparse matrix or open ondisc matrix (`odm`) with row names.
#' @param id Gene or other feature name, already checked by the caller.
#' @return Matching row numbers, starting at 1; `integer()` if none match.
#'   If row names repeat, all matching rows are returned.
#' @noRd
storage_response_index <- function(store, id) {
    which(id == rownames(store))
}

#' Read one gene's counts for the requested cells
#' @param store Row-accessible sparse matrix (`dgRMatrix`) or an `odm` opened
#'   in the current R process.
#' @param response_idx Gene's row number in `store`, starting at 1.
#' @param source_cells Original count-matrix column numbers, in the desired
#'   order.
#'   NULL selects all cells; `integer()` selects none. Repeats are kept.
#' @return Numeric count vector in the requested order, including for one cell.
#' @details Read the complete gene row before selecting its cells. This loads
#'   one
#'   row, rather than the whole matrix, and leaves the original counts
#'   unchanged.
#' @noRd
storage_read_response <- function(store, response_idx, source_cells = NULL) {
    y <- if (methods::is(store, "odm")) {
        require_ondisc_capabilities()
        as.numeric(store[response_idx, ])
    } else {
        load_csr_row(store@j, store@p, store@x, response_idx, ncol(store))
    }
    if (is.null(source_cells)) y else y[source_cells]
}

#' Find cells whose stored count reaches a feature's assignment threshold
#' @param store Row-accessible sparse matrix (`dgRMatrix`) or open `odm`.
#' @param feature_idx Feature's row number in `store`, starting at 1.
#' @param threshold Count threshold.
#' @return Original matrix column numbers, starting at 1, in storage order.
#' @details Test only entries actually stored in the sparse matrix. An unstored
#'   zero is never returned, even for threshold zero; an explicitly stored zero
#'   can pass. In-memory values and thresholds use doubles; ondisc uses integer
#'   conversion for both. The input matrix is unchanged.
#' @noRd
storage_threshold_row <- function(store, feature_idx, threshold) {
    if (methods::is(store, "odm")) {
        require_ondisc_capabilities()
        ondisc:::threshold_count_matrix_cpp(
            file_name_in = store@h5_file,
            f_row_ptr = store@ptr,
            row_idx = feature_idx,
            threshold = threshold
        )
    } else {
        threshold_count_matrix(
            store@j,
            store@p,
            store@x,
            feature_idx,
            as.numeric(threshold)
        )
    }
}

#' Copy cell-position vectors before passing them to ondisc
#' @param groups List of cell-position vectors, using the numbering the caller
#'   needs. This helper copies the numbers without changing their meaning.
#' @return Fresh integer vectors with names, order and repeats preserved.
#'   Some private ondisc functions subtract one from their input vectors in
#'   place.
#'   Give them copies so the original guide assignments remain unchanged.
#'   `as.integer()` alone may reuse the original vector and is not a safe copy.
#' @noRd
storage_owned_groups <- function(groups) {
    lapply(groups, function(group) as.integer(group) + 0L)
}

#' Count stored entries for each requested gene-guide pair
#' @param store Row-accessible sparse matrix (`dgRMatrix`) or open `odm`;
#'   genes are rows and original cells are columns.
#' @param pair_index Table giving each request's gene row (`response_idx`) and
#'   position in `groups` (`group_idx`), both starting at 1. Repeated requests
#'   and any request order are allowed.
#' @param groups Treatment-cell lists. Numbers are positions within
#'   `source_cells`,
#'   not original matrix columns. For `source_cells = c(2, 4, 5)`, group `c(1,
#'   3)`
#'   means original cells 2 and 5.
#' @param source_cells Original matrix columns kept after QC, in analysis order.
#' @param controls `"complement"` for all other kept cells, or NT-control
#'   positions
#'   within `source_cells`.
#' @param nt_groups Optional NT-guide cell lists, also as positions within
#'   `source_cells`. Low-MOI callers supply their NT groups for a summary across
#'   all genes. These lists are counted as supplied, without removing overlaps.
#'   NULL requests pair counts without retaining that summary.
#' @return `n_nonzero_trt` and `n_nonzero_cntrl`, integer vectors in request
#'   order.
#'   When retained, `n_nonzero_mat` has NT groups in rows and all genes in
#'   columns;
#'   `n_nonzero_tot` has a total per gene. Otherwise these two fields contain
#'   the
#'   existing empty placeholders `matrix()` and `integer()`.
#' @details These counts count stored entries, including explicitly stored
#'   zeros,
#'   rather than testing whether loaded values are positive. Complement totals
#'   count kept cells; NT totals sum the supplied NT groups, so an overlap would
#'   be counted again. Requests are sorted for reading, then restored to their
#'   original order. In-memory C++ code copies and converts its own positions.
#'   Private ondisc calls instead receive copied vectors and, where needed,
#'   positions within the NT-cell list. Caller inputs remain unchanged.
#' @noRd
storage_pair_counts <- function(
    store,
    pair_index,
    groups,
    source_cells,
    controls,
    nt_groups = NULL
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
        out <- storage_backed_pair_counts(
            store,
            response_idxs,
            group_idxs,
            groups,
            source_cells,
            all_nt_idxs,
            nt_groups,
            complement,
            summarize
        )
    } else {
        out <- storage_sparse_pair_counts(
            store,
            groups,
            source_cells,
            nt_groups,
            response_idxs,
            group_idxs,
            complement,
            summarize
        )
    }
    restore_order <- order(scan_order)
    out$n_nonzero_trt <- out$n_nonzero_trt[restore_order]
    out$n_nonzero_cntrl <- out$n_nonzero_cntrl[restore_order]
    out
}

#' Count stored entries for guide groups across a requested list of genes
#' @inheritParams storage_pair_counts
#' @param response_idxs Original gene-row numbers, starting at 1, in requested
#'   order.
#' @return Integer matrix with one row per group and one column per requested
#'   gene. Keep request order, repeats and matrix dimensions, even for zero or
#'   one row/column.
#' @details The C++ function counts all gene rows, then this helper selects the
#'   requested columns. Counts include explicitly stored zeros. Copy positions
#'   before private ondisc calls so the caller's vectors are not changed.
#' @noRd
storage_group_counts <- function(store, response_idxs, groups, source_cells) {
    if (methods::is(store, "odm")) {
        require_ondisc_capabilities()
        groups <- storage_owned_groups(groups)
        out <- ondisc:::compute_n_trt_cells_matrix_ondisc(
            file_name_in = store@h5_file,
            f_row_ptr = store@ptr,
            n_cells_orig = ncol(store),
            n_cells_sub = length(source_cells),
            n_genes = nrow(store),
            nt_grna_group_idxs = groups,
            cells_in_use = source_cells
        )
    } else {
        out <- sparse_group_counts(
            j = store@j,
            p = store@p,
            n_cells_orig = ncol(store),
            groups = groups,
            cells_in_use = source_cells
        )
    }
    out[, response_idxs, drop = FALSE]
}

#' Write response and guide matrices to ondisc files with the same import ID
#' @param response_matrix,grna_matrix Sparse count matrices prepared during
#'   import.
#' @param directory_to_write Output directory already created by the caller.
#' @param integer_id Import ID already randomly chosen by the caller.
#' @return Open ondisc matrices named `response_matrix` and `grna_matrix`.
#' @details Write `response.odm`, then `grna.odm`. If the second write fails,
#'   the
#'   first file may already exist. Use the private ondisc constructor to
#'   preserve
#'   the supplied shared ID; the public constructor fixes its ID at zero.
#'   ondisc's integer count conversion, chunk-size limits and errors still
#'   apply.
#'   Open file handles belong to the current R process.
#' @noRd
storage_create_backed_pair <- function(
    response_matrix,
    grna_matrix,
    directory_to_write,
    integer_id
) {
    require_ondisc_capabilities()
    constructor <- get(
        "create_odm_from_r_matrix_internal",
        asNamespace("ondisc")
    )
    # ondisc rejects its default 1000-entry HDF5 chunk for otherwise valid
    # small matrices. Admit both stores before creating either file, retain the
    # default for large stores, and explicitly lower only storage chunk size.
    # This does not subset counts, change values or recompute covariates.
    small_chunk <- function(matrix) {
        entries <- length(matrix@j)
        if (entries > 1000L) {
            return(NULL)
        }
        if (entries < 2L) {
            stop(
                "This ondisc constructor cannot create a store with fewer than two stored entries; use in-memory import for this screen."
            )
        }
        if (!"chunk_size" %in% names(formals(constructor))) {
            stop(
                "This ondisc constructor lacks the chunk_size argument required for a small count store."
            )
        }
        as.integer(entries - 1L)
    }
    response_chunk <- small_chunk(response_matrix)
    guide_chunk <- small_chunk(grna_matrix)
    create <- function(matrix, name, chunk) {
        args <- list(
            mat = matrix,
            file_to_write = file.path(directory_to_write, name),
            integer_id = integer_id
        )
        if (!is.null(chunk)) {
            args$chunk_size <- chunk
        }
        do.call(constructor, args)
    }
    # Both files must retain the caller's shared import ID so reopening can
    # check that they belong to this object. The public constructor uses zero.
    response_matrix <- create(response_matrix, "response.odm", response_chunk)
    grna_matrix <- create(grna_matrix, "grna.odm", guide_chunk)
    list(response_matrix = response_matrix, grna_matrix = grna_matrix)
}

#' Reopen both count files and check that they belong to the saved object
#' @param response_file,grna_file Paths to the response and guide ondisc files.
#' @param integer_id Import ID saved in the sceptre object.
#' @return Open ondisc matrices named `response_matrix` and `grna_matrix`.
#' @details Open both files in this R process, then compare their IDs with the
#'   saved ID. A mismatch raises the existing error. No files are rewritten.
#' @noRd
storage_reopen_backed_pair <- function(response_file, grna_file, integer_id) {
    require_ondisc_capabilities()
    response <- ondisc::initialize_odm_from_backing_file(response_file)
    grna <- ondisc::initialize_odm_from_backing_file(grna_file)
    if (integer_id != response@integer_id) {
        stop(
            "The `sceptre_object` and `response_odm` have distinct IDs. The `sceptre_object` likely is associated with a different backing .odm file."
        )
    }
    if (integer_id != grna@integer_id) {
        stop(
            "The `sceptre_object` and `grna_odm` have distinct IDs. The `sceptre_object` likely is associated with a different backing .odm file."
        )
    }
    list(response_matrix = response, grna_matrix = grna)
}

#' Calculate count summaries for each cell
#' @param store Column-accessible sparse matrix (`dgCMatrix`), with features in
#'   rows and cells in columns.
#' @param mt_gene_idxs Mitochondrial feature-row numbers, starting at 1; may be
#'   empty.
#' @param compute_p_mito,compute_max_feature Whether to calculate mitochondrial
#'   fractions and the largest feature's contribution.
#' @return One value per cell for `n_nonzero` (counts greater than 0.5),
#'   `n_umi`,
#'   `p_mito`, `frac_umis_max_feature` and `max_feature`. Unlike R row numbers,
#'   `max_feature` starts at 0; ties use the first maximum. Disabled optional
#'   summaries remain zero. A computed ratio with zero total counts can be
#'   non-finite. C++ copies mitochondrial positions before converting them,
#'   leaving the caller's vector unchanged.
#' @noRd
storage_cell_covariates <- function(
    store,
    mt_gene_idxs,
    compute_p_mito,
    compute_max_feature
) {
    compute_cell_covariates_cpp(
        i = store@i,
        p = store@p,
        x = store@x,
        n_genes = nrow(store),
        n_cells = ncol(store),
        mt_gene_idxs = mt_gene_idxs,
        compute_p_mito = compute_p_mito,
        compute_max_feature = compute_max_feature
    )
}

#' Convert a count matrix for efficient access by row or column
#' @param matrix_in Numeric dense matrix, numeric or logical sparse matrix, or
#'   open ondisc matrix (`odm`). Sparse input may use row, column or triplet
#'   storage.
#' @param make_row_accessible TRUE requests row storage (`dgRMatrix`); FALSE
#'   requests column storage (`dgCMatrix`).
#' @return Sparse matrix with the same dimensions and names, or the supplied
#'   `odm`. A numeric sparse matrix already in the requested form is returned as
#'   is.
#' @details Keep stored numeric values, including zeros and fractions. For
#'   logical
#'   sparse input, every stored entry becomes one, even an explicitly stored
#'   FALSE;
#'   this preserves the package's existing convention. Duplicate triplet entries
#'   are not combined here. Sorting uses copies so callers' sparse vectors are
#'   unchanged. Dense input is first treated as a general matrix to retain both
#'   triangles of symmetric-looking data; this can require an extra dense copy.
#'   Sparse and ondisc input are not converted to dense matrices. In-memory
#'   conversion works without installing ondisc.
#' @noRd
set_matrix_accessibility <- function(matrix_in, make_row_accessible = TRUE) {
    owns_i <- owns_j <- owns_x <- FALSE
    # For logical sparse input, count every stored entry as 1, even stored
    # FALSE.
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
                methods::as(matrix_in, "generalMatrix"),
                "TsparseMatrix"
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

        # Sorting changes vectors in place. Copy vectors taken from the input
        # sparse matrix so its slots stay unchanged. Dense conversion, pointer
        # expansion and replacing logical values already made fresh vectors;
        # those can be sorted directly. Do not guess whether copying is needed
        # from R's reference counts.
        dt <- list(
            i = if (owns_i) i else data.table::copy(i),
            j = if (owns_j) j else data.table::copy(j),
            x = if (owns_x) x else data.table::copy(x)
        )
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


#' Expand sparse-matrix entry counts into a row or column number per entry
#' @param p Compressed sparse pointer: cumulative entry counts starting at zero,
#'   with length equal to the number of rows/columns plus one. Counts must not
#'   decrease, and the vector must describe at least one row/column.
#' @return Row/column numbers starting at zero, one per stored entry. For
#'   example,
#'   `p = c(0, 2, 2, 3)` gives `c(0, 0, 2)`: two entries in the first
#'   row/column,
#'   none in the second and one in the third. The input and entry order are
#'   kept.
#' @noRd
sparse_pointer_to_indices <- function(p) {
    l <- diff(p)
    rep(x = seq(0, length(l) - 1L), times = l)
}

# Retain the sentinel request required by the qualified ondisc ABI.
storage_backed_pair_counts <- function(
    store,
    response_idxs,
    group_idxs,
    groups,
    source_cells,
    all_nt_idxs,
    nt_groups,
    complement,
    summarize
) {
    require_ondisc_capabilities()
    if (!complement) {
        nt_groups <- lapply(nt_groups, function(group) {
            match(group, all_nt_idxs)
        })
    }
    groups <- storage_owned_groups(groups)
    nt_groups <- storage_owned_groups(nt_groups)
    # Some installed ondisc binaries inspect the next requested row before
    # checking whether the final pair was consumed. Keep an inert request
    # beyond the gene domain so that lookahead remains within the vector.
    # Neither native loop visits this row; discard its unused output below.
    # The original unique rows are still the only rows requested for reading.
    n_requested <- length(response_idxs)
    native_response_idxs <- c(response_idxs, as.integer(nrow(store) + 1L))
    native_group_idxs <- c(group_idxs, 1L)
    common <- storage_backed_pair_arguments(
        store,
        groups,
        source_cells,
        all_nt_idxs,
        native_response_idxs,
        native_group_idxs,
        complement
    )
    out <- invoke_backed_pair_count(common, response_idxs, nt_groups, summarize)
    out$n_nonzero_trt <- out$n_nonzero_trt[seq_len(n_requested)]
    out$n_nonzero_cntrl <- out$n_nonzero_cntrl[seq_len(n_requested)]
    out
}

storage_backed_pair_arguments <- function(
    store,
    groups,
    source_cells,
    all_nt_idxs,
    native_response_idxs,
    native_group_idxs,
    complement
) {
    list(
        file_name_in = store@h5_file,
        f_row_ptr = store@ptr,
        n_genes = nrow(store),
        n_cells_orig = ncol(store),
        n_cells_sub = length(source_cells),
        grna_group_idxs = groups,
        all_nt_idxs = all_nt_idxs,
        to_analyze_response_idxs = native_response_idxs,
        to_analyze_grna_idxs = native_group_idxs,
        control_group_complement = complement,
        cells_in_use = source_cells
    )
}

storage_sparse_pair_counts <- function(
    store,
    groups,
    source_cells,
    nt_groups,
    response_idxs,
    group_idxs,
    complement,
    summarize
) {
    out <- sparse_pair_counts(
        j = store@j,
        p = store@p,
        n_cells_orig = ncol(store),
        grna_group_idxs = groups,
        nt_groups_analysis = nt_groups,
        to_analyze_response_idxs = response_idxs,
        to_analyze_grna_idxs = group_idxs,
        control_group_complement = complement,
        cells_in_use = source_cells
    )
    if (!summarize) {
        out$n_nonzero_mat <- matrix()
        out$n_nonzero_tot <- integer()
    }
    out
}

invoke_backed_pair_count <- function(
    common,
    response_idxs,
    nt_groups,
    summarize
) {
    if (!summarize) {
        out <- do.call(
            ondisc:::compute_n_ok_pairs_ondisc,
            c(
                common,
                list(
                    unique_response_idxs = sort(unique(response_idxs))
                )
            )
        )
        out$n_nonzero_mat <- matrix()
        out$n_nonzero_tot <- integer()
    } else {
        out <- do.call(
            ondisc:::compute_nt_nonzero_matrix_and_n_ok_pairs_ondisc,
            c(common, list(indiv_nt_grna_idxs = nt_groups))
        )
    }
    out
}
