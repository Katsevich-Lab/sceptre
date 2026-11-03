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
# backed counts are copied into it block by block with HDF5Array's 10x
# writer (TENxRealizationSink), whose column-compressed groups are AnnData's
# csc_matrix once the AnnData encoding attributes are added.
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

# row_blocks() for counts read from screen.h5mu (location:
# h5_counts_location(mat)): blocks of consecutive rows holding about as many
# entries as a block of DelayedArray's automatic size holds elements (or one
# row), each read with one read of the indices and one of the values.
h5_row_blocks <- function(mat, location, f) {
    widths <- diff(location$indptr)
    if (length(widths) == 0L) {
        return(list())
    }
    budget <- max(1, DelayedArray::getAutoBlockSize() / 8)
    first <- location$indptr[-length(location$indptr)]
    block_of_row <- c(0, cumsum(diff(floor(first / budget)) > 0))
    file <- rhdf5::H5Fopen(location$file, flags = "H5F_ACC_RDONLY")
    on.exit(rhdf5::H5Fclose(file), add = TRUE)
    datasets <- lapply(
        c(indices = "indices", data = "data"),
        function(name) rhdf5::H5Dopen(file, paste0(location$group, "/", name))
    )
    on.exit(lapply(datasets, rhdf5::H5Dclose), add = TRUE, after = FALSE)
    # entries start, ..., start + count - 1 of a dataset, as a vector of
    # `type`; read with H5Dread(), since h5read() copies the values to drop
    # their dimension and scans them for missing values
    read <- function(dataset, start, count, type) {
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
    lapply(unname(split(seq_along(widths), block_of_row)), function(rows) {
        start <- first[rows[1L]] + 1
        count <- sum(widths[rows])
        block <- methods::new(
            "dgRMatrix",
            j = read(datasets$indices, start, count, "integer"),
            p = as.integer(c(0, cumsum(widths[rows]))),
            x = read(datasets$data, start, count, "double"),
            Dim = c(length(rows), ncol(mat)),
            Dimnames = list(rownames(mat)[rows], colnames(mat))
        )
        f(block, rows[1L] - 1L)
    })
}

# Writes counts (features x cells; any matrix that DelayedArray can read) into
# an existing HDF5 file as the AnnData csc_matrix group `group`, block by
# block, so that the counts are never all in memory at once. In-memory counts
# are sliced as a dgRMatrix, several times faster than through DelayedArray.
# `level` is the deflate level (HDF5Array's default; 0 for a temporary file).
write_backed_counts <- function(counts, file, group, level = NULL) {
    if (!is_backed(counts)) {
        counts <- set_matrix_accessibility(counts, make_row_accessible = TRUE)
    }
    template <- DelayedArray::DelayedArray(counts)
    sink <- HDF5Array::TENxRealizationSink(
        dim = rev(dim(counts)),
        type = DelayedArray::type(template),
        filepath = file,
        group = group,
        level = if (is.null(level)) HDF5Array::getHDF5DumpCompressionLevel() else level
    )
    # blocks of consecutive features; each block's transpose (cells x
    # features) is the same CSR data read as a dgCMatrix
    grid <- DelayedArray::colAutoGrid(sink)
    for (k in seq_len(ncol(grid))) {
        viewport <- grid[[1L, k]]
        rows <- viewport_range(viewport, 2L)
        block <- compressed_sparse_components(row_slice(counts, rows))
        DelayedArray::write_block(
            sink,
            viewport,
            methods::new(
                "dgCMatrix",
                i = block$j,
                p = block$p,
                x = block$x,
                Dim = c(ncol(counts), length(rows))
            )
        )
    }
    DelayedArray::close(sink)
    # the 10x layout keeps the shape in a dataset; AnnData in attributes
    rhdf5::h5delete(file, paste0(group, "/shape"))
    h5 <- rhdf5::H5Fopen(file)
    on.exit(rhdf5::H5Fclose(h5), add = TRUE)
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
    rhdf5::h5writeAttribute(rev(dim(counts)), x, "shape")
    invisible(file)
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

