# Internal matrix boundary. ondisc stays optional: its namespace is only used
# for an actual backed operation. The selected route centralizes its existing
# private APIs here rather than requiring a companion package interface change.

storage_response_index <- function(store, id) {
    which(id == rownames(store))
}

storage_read_response <- function(store, response_idx, source_cells = NULL) {
    y <- if (methods::is(store, "odm")) {
        as.numeric(store[response_idx, ])
    } else {
        load_csr_row(store@j, store@p, store@x, response_idx, ncol(store))
    }
    if (is.null(source_cells)) y else y[source_cells]
}

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

# Some legacy native kernels decrement their R input vectors in place. Integer
# coercion alone may alias an existing integer vector; addition owns fresh data.
storage_owned_groups <- function(groups) {
    lapply(groups, function(group) as.integer(group) + 0L)
}

# groups and controls use positions in source_cells, never original columns.
# Counts intentionally mean stored sparse positions, including explicit zeros.
# Supplying nt_groups retains the existing fused all-response NT summary scan.
# Pair rows can be unordered or repeated; the scan receives response-sorted rows
# and its results are restored to precisely the original request order.
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
    if (!complement) {
        nt_groups <- lapply(nt_groups, function(group) match(group, all_nt_idxs))
    }
    scan_order <- order(pair_index$response_idx, method = "radix")
    response_idxs <- as.integer(pair_index$response_idx[scan_order])
    group_idxs <- as.integer(pair_index$group_idx[scan_order])
    groups <- storage_owned_groups(groups)
    nt_groups <- storage_owned_groups(nt_groups)
    if (methods::is(store, "odm")) {
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
        out <- compute_nt_nonzero_matrix_and_n_ok_pairs_v3(
            j = store@j, p = store@p, n_cells_orig = ncol(store),
            n_cells_sub = length(source_cells), grna_group_idxs = groups,
            indiv_nt_grna_idxs = nt_groups, all_nt_idxs = all_nt_idxs,
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

# The current calibration consumer needs the legacy groups-by-all-responses
# integer matrix. Keep that batched native scan and subset its columns only at
# this boundary; do not replace it with per-cell R loops or eager dense counts.
storage_group_counts <- function(store, response_idxs, groups, source_cells) {
    groups <- storage_owned_groups(groups)
    if (methods::is(store, "odm")) {
        out <- ondisc:::compute_n_trt_cells_matrix_ondisc(
            file_name_in = store@h5_file, f_row_ptr = store@ptr,
            n_cells_orig = ncol(store), n_cells_sub = length(source_cells),
            n_genes = nrow(store), nt_grna_group_idxs = groups,
            cells_in_use = source_cells
        )
    } else {
        out <- compute_n_trt_cells_matrix(
            j = store@j, p = store@p, n_cells_orig = ncol(store),
            n_cells_sub = length(source_cells), n_genes = nrow(store),
            nt_grna_group_idxs = groups, cells_in_use = source_cells
        )
    }
    out[, response_idxs, drop = FALSE]
}

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

storage_cell_covariates <- function(
    store, mt_gene_idxs, compute_p_mito, compute_max_feature
) {
    compute_cell_covariates_cpp(
        i = store@i, p = store@p, x = store@x,
        n_genes = nrow(store), n_cells = ncol(store),
        mt_gene_idxs = as.integer(mt_gene_idxs) + 0L,
        compute_p_mito = compute_p_mito,
        compute_max_feature = compute_max_feature
    )
}

# This conversion remains package-local so sparse workflows never require
# ondisc. Its legacy ondisc counterpart has the same characterization tests;
# sharing a binary implementation would create an inappropriate dependency.
set_matrix_accessibility <- function(matrix_in, make_row_accessible = TRUE) {
    # if a logical matrix, convert to the corresponding numeric matrix; consider
    # more efficient implementation later
    if (methods::is(matrix_in, "lgRMatrix")) {
        attr(matrix_in, "class") <- "dgRMatrix"
        matrix_in@x <- rep(1.0, length(matrix_in@j))
    }
    if (methods::is(matrix_in, "lgCMatrix")) {
        attr(matrix_in, "class") <- "dgCMatrix"
        matrix_in@x <- rep(1.0, length(matrix_in@i))
    }
    if (methods::is(matrix_in, "lgTMatrix")) {
        attr(matrix_in, "class") <- "dgTMatrix"
        matrix_in@x <- rep(1.0, length(matrix_in@j))
    }

    if (methods::is(matrix_in, "dgRMatrix") && make_row_accessible) {
        out <- matrix_in
    } else if (methods::is(matrix_in, "dgCMatrix") && !make_row_accessible) {
        out <- matrix_in
    } else if (methods::is(matrix_in, "odm")) {
        out <- matrix_in
    } else {
        # basic info: get matrix_in, n_responses, response_ids, and n_cells
        n_responses <- nrow(matrix_in)
        response_ids <- rownames(matrix_in)
        cell_ids <- colnames(matrix_in)
        n_cells <- ncol(matrix_in)
        if (methods::is(matrix_in, "dgTMatrix")) {
            i <- matrix_in@i
            j <- matrix_in@j
        } else if (methods::is(matrix_in, "matrix")) {
            matrix_in <- methods::as(matrix_in, "TsparseMatrix")
            i <- matrix_in@i
            j <- matrix_in@j
        } else if (methods::is(matrix_in, "dgCMatrix")) {
            i <- matrix_in@i
            j <- convert_pointer_to_index_vector_v2(matrix_in@p)
        } else if (methods::is(matrix_in, "dgRMatrix")) {
            i <- convert_pointer_to_index_vector_v2(matrix_in@p)
            j <- matrix_in@j
        } else {
            stop("Class not recognized.")
        }
        x <- matrix_in@x

        # perform the sorting operation
        sort_in_place <- methods::is(matrix_in, "matrix") ||
            methods::is(matrix_in, "dgTMatrix")
        if (sort_in_place) {
            l <- list(i = i, j = j, x = x)
            dt <- data.table::setDT(l) # assign by reference
        } else {
            dt <- data.table::data.table(i = i, j = j, x = x) # deep copy
        }

        # finally, initialize the output
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


convert_pointer_to_index_vector_v2 <- function(p) {
    l <- diff(p)
    rep(x = seq(0, length(l) - 1L), times = l)
}

