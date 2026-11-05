# Response counts kept on disk ("backed"): the counts stay in screen.h5mu and
# are read through a DelayedMatrix (HDF5Array) as they are needed, instead of
# being loaded into memory. The gRNA counts, which are much smaller and are
# scanned whole by every gRNA assignment method, are always held in memory.
#
# Layout. AnnData stores a modality's X as cells x features. sceptre reads one
# response's counts across all cells at a time, so X is stored as
# `csc_matrix`: each feature's nonzero entries are contiguous. MuData's
# backed reader opens this as a features x cells CSR_H5ADMatrixSeed, whose row
# reads touch only that row's entries, and Python's anndata and mudata read the
# same file, in memory or backed. MuData writes this layout from a dgRMatrix;
# the response counts are written into it block by block with rhdf5
# (write_backed_counts()), as 32-bit integer counts.
#
# Analysis code reads backed counts in two ways: one response at a time
# (load_row()), and in blocks of consecutive responses realized as in-memory
# dgRMatrix (row_blocks()), so that code written for in-memory CSR counts can
# scan every response with bounded memory. Block sizes follow DelayedArray's
# automatic block size (DelayedArray::setAutoBlockSize()). Blocks of counts
# read straight from screen.h5mu hold as many nonzero entries as a block holds
# elements, and are read from the file's index and value arrays directly:
# realizing them through DelayedArray's sparse coercions costs several times
# the read itself.
#
# The blocks are assembled from their compressed sparse components with
# methods::new() on the slots that Matrix documents for dgRMatrix and
# dgCMatrix (?"dgRMatrix-class"), rather than with Matrix::sparseMatrix(),
# which would expand and re-validate every block; the components are read
# back only through compressed_sparse_components().

# the group of screen.h5mu holding the response counts
response_counts_group <- "/mod/response/X"

is_backed <- function(x) {
    methods::is(x, "DelayedMatrix")
}

# The counts of one response across all cells, as a numeric vector.
backed_row <- function(mat, id) {
    i <- match(id, rownames(mat))
    as.numeric(DelayedArray::extract_array(mat, list(i, NULL)))
}

# Counts as an in-memory dgRMatrix. A DelayedMatrix is realized through its
# transpose, whose column-compressed form is the CSR form of the counts.
as_row_counts <- function(x) {
    if (!is_backed(x)) {
        return(set_matrix_accessibility(x, make_row_accessible = TRUE))
    }
    tx <- compressed_sparse_components(
        methods::as(methods::as(DelayedArray::t(x), "CsparseMatrix"), "dMatrix")
    )
    methods::new(
        "dgRMatrix",
        j = tx$i,
        p = tx$p,
        x = tx$x,
        Dim = dim(x),
        Dimnames = if (is.null(dimnames(x))) list(NULL, NULL) else dimnames(x)
    )
}

# The indices that a block (viewport) of a grid spans along dimension `along`.
viewport_range <- function(viewport, along) {
    seq(
        DelayedArray::start(viewport)[along],
        DelayedArray::end(viewport)[along]
    )
}

# Rows of a DelayedMatrix, or of an in-memory dgRMatrix, as a dgRMatrix.
row_slice <- function(mat, rows) {
    as_row_counts(mat[rows, , drop = FALSE])
}

# Calls f(block, offset) on consecutive blocks of rows of `mat`, each an
# in-memory dgRMatrix holding rows offset + 1, ..., offset + nrow(block), and
# returns the results as a list. An in-memory dgRMatrix is passed whole.
row_blocks <- function(mat, f) {
    if (!is_backed(mat)) {
        return(list(f(mat, 0L)))
    }
    location <- h5_counts_location(mat)
    if (!is.null(location)) {
        return(h5_row_blocks(mat, location, f))
    }
    grid <- DelayedArray::rowAutoGrid(mat)
    lapply(seq_len(length(grid)), function(k) {
        rows <- viewport_range(grid[[k, 1L]], 1L)
        f(row_slice(mat, rows), rows[1L] - 1L)
    })
}

# Where in screen.h5mu backed counts read from it are: the response modality's
# X, an AnnData csc_matrix that MuData's backed reader opens as a features x
# cells CSR_H5ADMatrixSeed, under no delayed operation but names; NULL for
# other counts. The group's shape and number of entries must be the counts'.
h5_counts_location <- function(mat) {
    if (
        !is_backed(mat) ||
            !DelayedArray::isPristine(mat, ignore.dimnames = TRUE) ||
            !methods::is(DelayedArray::seed(mat), "CSR_H5ADMatrixSeed")
    ) {
        return(NULL)
    }
    file <- DelayedArray::path(mat)
    group <- response_counts_group
    indptr <- tryCatch(
        {
            shape <- rhdf5::h5readAttributes(file, group)$shape
            stopifnot(identical(as.numeric(shape), as.numeric(rev(dim(mat)))))
            as.numeric(rhdf5::h5read(
                file,
                paste0(group, "/indptr"),
                bit64conversion = "double"
            ))
        },
        error = function(e) NULL
    )
    if (
        length(indptr) != nrow(mat) + 1L ||
            indptr[length(indptr)] != HDF5Array::nzcount(DelayedArray::seed(mat))
    ) {
        return(NULL)
    }
    list(file = file, group = group, indptr = indptr)
}

# Blocks of consecutive features, given the number of entries of each
# (`widths`): each block holds about as many entries as a block of
# DelayedArray's automatic size holds elements, or a single feature.
entry_blocks <- function(widths) {
    if (length(widths) == 0L) {
        return(list())
    }
    budget <- max(1, DelayedArray::getAutoBlockSize() / 8)
    first <- c(0, cumsum(as.numeric(widths)))[seq_along(widths)]
    block_of_row <- c(0, cumsum(diff(floor(first / budget)) > 0))
    unname(split(seq_along(widths), block_of_row))
}

# Entries start, ..., start + count - 1 of an open HDF5 dataset, as a vector
# of `type`; read with H5Dread(), since h5read() copies the values to drop
# their dimension and scans them for missing values.
h5_read_entries <- function(dataset, start, count, type) {
    if (count == 0) {
        return(vector(type))
    }
    space <- rhdf5::H5Dget_space(dataset)
    on.exit(rhdf5::H5Sclose(space), add = TRUE)
    rhdf5::H5Sselect_hyperslab(space, start = start, count = count)
    memory <- rhdf5::H5Screate_simple(count)
    on.exit(rhdf5::H5Sclose(memory), add = TRUE)
    values <- rhdf5::H5Dread(
        dataset,
        h5spaceFile = space,
        h5spaceMem = memory,
        drop = TRUE
    )
    if (typeof(values) == type) values else methods::as(values, type)
}

# row_blocks() for counts read from screen.h5mu (location:
# h5_counts_location(mat)): blocks of consecutive rows (entry_blocks()), each
# read with one read of the indices and one of the values.
h5_row_blocks <- function(mat, location, f) {
    widths <- diff(location$indptr)
    if (length(widths) == 0L) {
        return(list())
    }
    first <- location$indptr[-length(location$indptr)]
    file <- rhdf5::H5Fopen(location$file, flags = "H5F_ACC_RDONLY")
    on.exit(rhdf5::H5Fclose(file), add = TRUE)
    datasets <- lapply(
        c(indices = "indices", data = "data"),
        function(name) rhdf5::H5Dopen(file, paste0(location$group, "/", name))
    )
    on.exit(lapply(datasets, rhdf5::H5Dclose), add = TRUE, after = FALSE)
    lapply(entry_blocks(widths), function(rows) {
        start <- first[rows[1L]] + 1
        count <- sum(widths[rows])
        block <- methods::new(
            "dgRMatrix",
            j = h5_read_entries(datasets$indices, start, count, "integer"),
            p = as.integer(c(0, cumsum(widths[rows]))),
            x = h5_read_entries(datasets$data, start, count, "double"),
            Dim = c(length(rows), ncol(mat)),
            Dimnames = list(rownames(mat)[rows], colnames(mat))
        )
        f(block, rows[1L] - 1L)
    })
}

# Counts to write with write_backed_counts(), other than a matrix: `dim`
# (features x cells), the number of entries of each feature (`widths`), and
# `block(rows)`, the counts of the consecutive features `rows` as a cells x
# features dgCMatrix.
counts_source <- function(dim, widths, block) {
    structure(
        list(dim = dim, widths = widths, block = block),
        class = "sceptre_counts_source"
    )
}

# A matrix as a counts source. In-memory counts are sliced as a dgRMatrix,
# whose rows are the columns of the transpose, several times faster than
# through DelayedArray; other counts are realized through DelayedArray, in
# blocks of its automatic grid (`widths` NULL).
matrix_counts_source <- function(counts) {
    if (is_backed(counts)) {
        return(counts_source(
            dim = dim(counts),
            widths = NULL,
            block = function(rows) {
                block <- compressed_sparse_components(row_slice(counts, rows))
                methods::new(
                    "dgCMatrix",
                    i = block$j,
                    p = block$p,
                    x = block$x,
                    Dim = c(ncol(counts), length(rows))
                )
            }
        ))
    }
    counts <- set_matrix_accessibility(counts, make_row_accessible = TRUE)
    csr <- compressed_sparse_components(counts)
    counts_source(
        dim = dim(counts),
        widths = diff(csr$p),
        block = function(rows) {
            first <- csr$p[rows[1L]]
            entries <- seq_len(csr$p[rows[length(rows)] + 1L] - first) + first
            methods::new(
                "dgCMatrix",
                i = csr$j[entries],
                p = csr$p[c(rows, rows[length(rows)] + 1L)] - first,
                x = csr$x[entries],
                Dim = c(ncol(counts), length(rows))
            )
        }
    )
}

# The columns `feature_cols` (increasing) of in-memory counts held as cells x
# features (a dgCMatrix, `t_counts`) as a counts source of their features.
# With `check_name`, each block is checked to hold counts
# (check_integer_counts()) under that name.
t_counts_source <- function(t_counts, feature_cols, check_name = NULL) {
    csc <- compressed_sparse_components(t_counts)
    counts_source(
        dim = c(length(feature_cols), nrow(t_counts)),
        widths = diff(csc$p)[feature_cols],
        block = function(rows) {
            block <- t_counts[, feature_cols[rows], drop = FALSE]
            if (!is.null(check_name)) {
                check_integer_counts(block, check_name)
            }
            block
        }
    )
}

# Counts held in HDF5 files written by write_backed_counts() (group `group`),
# one per set of cells, as a counts source of all the cells, the cells of each
# file following those of the files before it. Each block is read from the
# files' index and value arrays directly.
gathered_counts_source <- function(files, group) {
    indptrs <- lapply(files, function(file) {
        as.numeric(rhdf5::h5read(
            file,
            paste0(group, "/indptr"),
            bit64conversion = "double"
        ))
    })
    n_cells <- vapply(
        files,
        function(file) as.integer(rhdf5::h5readAttributes(file, group)$shape[1L]),
        integer(1)
    )
    cell_offsets <- c(0L, cumsum(n_cells))[seq_along(files)]
    n_features <- length(indptrs[[1L]]) - 1L
    counts_source(
        dim = c(n_features, sum(n_cells)),
        widths = Reduce(`+`, lapply(indptrs, diff)),
        block = function(rows) {
            pieces <- Map(
                function(file, indptr) {
                    h5 <- rhdf5::H5Fopen(file, flags = "H5F_ACC_RDONLY")
                    on.exit(rhdf5::H5Fclose(h5), add = TRUE)
                    first <- indptr[rows[1L]]
                    count <- indptr[rows[length(rows)] + 1L] - first
                    read <- function(name, type) {
                        dataset <- rhdf5::H5Dopen(h5, paste0(group, "/", name))
                        on.exit(rhdf5::H5Dclose(dataset), add = TRUE)
                        h5_read_entries(dataset, first + 1, count, type)
                    }
                    list(
                        i = read("indices", "integer"),
                        x = read("data", "double"),
                        p = as.integer(indptr[c(rows, rows[length(rows)] + 1L)] - first)
                    )
                },
                files,
                indptrs
            )
            # within a feature, the cells of each file in turn
            stacked <- stack_csc_rows(
                i_list = lapply(pieces, `[[`, "i"),
                x_list = lapply(pieces, `[[`, "x"),
                p_list = lapply(pieces, `[[`, "p"),
                row_offsets = cell_offsets,
                n_cols = length(rows)
            )
            rm(pieces)
            methods::new(
                "dgCMatrix",
                i = stacked$i,
                p = stacked$p,
                x = stacked$x,
                Dim = c(sum(n_cells), length(rows))
            )
        }
    )
}

# Writes counts (features x cells: any matrix that DelayedArray can read, or a
# counts source) into an existing HDF5 file as the AnnData csc_matrix group
# `group`, block by block, so that the counts are never all in memory at once.
# Blocks of consecutive features hold a bounded number of entries
# (entry_blocks()) when the number of entries of each feature is known, and
# follow DelayedArray's automatic grid otherwise. Each block's indices and
# values are appended to the group's datasets with one write each; the values
# are counts, stored as 32-bit integers (as Cell Ranger stores them). `level`
# is the deflate level (HDF5Array's default; 0 for a temporary file).
write_backed_counts <- function(counts, file, group, level = NULL) {
    source <- if (inherits(counts, "sceptre_counts_source")) {
        counts
    } else {
        matrix_counts_source(counts)
    }
    dim <- source$dim
    level <- if (is.null(level)) HDF5Array::getHDF5DumpCompressionLevel() else level
    blocks <- if (is.null(source$widths)) {
        grid <- DelayedArray::rowAutoGrid(counts)
        lapply(seq_len(nrow(grid)), function(k) viewport_range(grid[[k, 1L]], 1L))
    } else {
        entry_blocks(source$widths)
    }
    if (!file.exists(file)) {
        rhdf5::h5createFile(file)
    }
    h5 <- rhdf5::H5Fopen(file)
    on.exit(rhdf5::H5Fclose(h5), add = TRUE)
    rhdf5::h5createGroup(h5, group)
    datasets <- lapply(
        c(indices = "indices", data = "data"),
        function(name) {
            h5_create_extendible(h5, paste0(group, "/", name), "H5T_NATIVE_INT32", level)
        }
    )
    on.exit(lapply(datasets, rhdf5::H5Dclose), add = TRUE, after = FALSE)
    widths <- numeric(dim[1L])
    n_entries <- 0
    for (rows in blocks) {
        block <- compressed_sparse_components(source$block(rows))
        values <- block$x
        if (any(values != round(values)) || any(abs(values) > .Machine$integer.max)) {
            stop("Only integer counts can be written to screen.h5mu.")
        }
        n <- length(values)
        h5_append_entries(datasets$indices, block$i, n_entries, n)
        h5_append_entries(datasets$data, as.integer(values), n_entries, n)
        widths[rows] <- diff(block$p)
        n_entries <- n_entries + n
    }
    rhdf5::h5createDataset(
        h5,
        paste0(group, "/indptr"),
        dims = dim[1L] + 1,
        H5type = "H5T_NATIVE_INT64",
        chunk = min(4096, dim[1L] + 1),
        level = level
    )
    rhdf5::h5write(c(0, cumsum(widths)), h5, paste0(group, "/indptr"))
    x <- rhdf5::H5Gopen(h5, group)
    on.exit(rhdf5::H5Gclose(x), add = TRUE, after = FALSE)
    rhdf5::h5writeAttribute(
        "csc_matrix",
        x,
        "encoding-type",
        variableLengthString = TRUE,
        asScalar = TRUE
    )
    rhdf5::h5writeAttribute(
        "0.1.0",
        x,
        "encoding-version",
        variableLengthString = TRUE,
        asScalar = TRUE
    )
    rhdf5::h5writeAttribute(rev(dim), x, "shape")
    invisible(file)
}

# An empty, extendible, chunked one-dimensional dataset of `H5type`, open.
h5_create_extendible <- function(h5, name, H5type, level, chunk = 16384) {
    rhdf5::h5createDataset(
        h5,
        name,
        dims = 0,
        maxdims = rhdf5::H5Sunlimited(),
        H5type = H5type,
        chunk = chunk,
        level = level
    )
    rhdf5::H5Dopen(h5, name)
}

# Appends `count` values to an open one-dimensional dataset holding `offset`.
h5_append_entries <- function(dataset, values, offset, count) {
    if (count == 0) {
        return(invisible())
    }
    rhdf5::H5Dset_extent(dataset, offset + count)
    space <- rhdf5::H5Dget_space(dataset)
    on.exit(rhdf5::H5Sclose(space), add = TRUE)
    rhdf5::H5Sselect_hyperslab(space, start = offset + 1, count = count)
    memory <- rhdf5::H5Screate_simple(count)
    on.exit(rhdf5::H5Sclose(memory), add = TRUE)
    rhdf5::H5Dwrite(dataset, values, h5spaceMem = memory, h5spaceFile = space)
    invisible()
}

# An empty dgRMatrix of the given dimensions and names: the placeholder that
# stands for the response counts while they are written separately
# (write_screen()), kept in sceptre_object.rds in place of the counts in
# screen.h5mu, or combined with the other matrices of a backed import.
empty_row_counts <- function(dim, dimnames = NULL) {
    Matrix::sparseMatrix(
        i = integer(),
        j = integer(),
        x = numeric(),
        dims = dim,
        dimnames = dimnames,
        repr = "R"
    )
}

empty_counts <- function(counts) {
    empty_row_counts(dim(counts), dimnames(counts))
}

# The HDF5 files that back `x`, if any.
backing_files <- function(x) {
    if (!is_backed(x)) {
        return(character())
    }
    files <- tryCatch(
        unlist(DelayedArray::seedApply(x, DelayedArray::path)),
        error = function(e) character()
    )
    unique(normalizePath(files, winslash = "/", mustWork = FALSE))
}

