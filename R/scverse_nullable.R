# MuData writes R logical masks as integer datasets with storage.mode metadata.
# Native AnnData nullable arrays require actual HDF5 booleans. Rewrite only
# nullable masks and boolean values as the standard FALSE=0/TRUE=1 enum; raw
# count datasets are never selected or read. Integer payloads keep their type.
# Remove R's storage.mode marker on those enums: rhdf5 otherwise attempts to
# change the storage mode of a factor before MuData can decode the boolean.
write_nullable_boolean_arrays <- function(file) {
    entries <- rhdf5::h5ls(file, recursive = TRUE)
    groups <- entries[entries$otype == "H5I_GROUP", , drop = FALSE]
    # rhdf5 exports no H5Tclose(); the datatype handle is released with the
    # other handles when rhdf5 closes them (h5closeAll())
    boolean_type <- NULL
    for (i in seq_len(nrow(groups))) {
        path <- paste0(sub("/$", "", groups$group[[i]]), "/", groups$name[[i]])
        encoding <- rhdf5::h5readAttributes(file, path)[["encoding-type"]]
        if (
            is.null(encoding) ||
                !encoding %in% c("nullable-integer", "nullable-boolean")
        ) {
            next
        }
        members <- if (encoding == "nullable-boolean") {
            c("mask", "values")
        } else {
            "mask"
        }
        for (member in members) {
            target <- paste0(path, "/", member)
            value <- rhdf5::h5read(file, target)
            dimensions <- dim(value)
            if (is.null(dimensions)) {
                dimensions <- length(value)
            }
            value <- as.logical(value)
            if (anyNA(value)) {
                stop("A nullable boolean payload may not contain NA")
            }
            attributes <- rhdf5::h5readAttributes(file, target)
            attributes[["storage.mode"]] <- NULL
            if (is.null(boolean_type)) {
                boolean_type <- rhdf5::H5Tenum_create("H5T_NATIVE_UCHAR")
                rhdf5::H5Tenum_insert(boolean_type, "FALSE", 0L)
                rhdf5::H5Tenum_insert(boolean_type, "TRUE", 1L)
            }
            # the handles are closed in reverse order of opening (`after =
            # FALSE`), so that the file is closed last
            replace_dataset <- function() {
                fid <- rhdf5::H5Fopen(file)
                on.exit(rhdf5::H5Fclose(fid), add = TRUE)
                rhdf5::H5Ldelete(fid, target)
                space <- rhdf5::H5Screate_simple(dimensions)
                on.exit(rhdf5::H5Sclose(space), add = TRUE, after = FALSE)
                settings <- rhdf5::H5Pcreate("H5P_DATASET_CREATE")
                on.exit(
                    rhdf5::H5Pclose(settings),
                    add = TRUE,
                    after = FALSE
                )
                if (length(value)) {
                    chunk <- rep.int(1L, length(dimensions))
                    chunk[[1L]] <- min(dimensions[[1L]], 65536L)
                    rhdf5::H5Pset_chunk(settings, chunk)
                    rhdf5::H5Pset_deflate(settings, 6L)
                }
                dataset <- rhdf5::H5Dcreate(
                    fid,
                    target,
                    boolean_type,
                    space,
                    dcpl = settings
                )
                on.exit(rhdf5::H5Dclose(dataset), add = TRUE, after = FALSE)
                rhdf5::H5Dwrite(
                    dataset,
                    as.raw(as.integer(value)),
                    h5type = boolean_type
                )
                for (name in names(attributes)) {
                    rhdf5::h5writeAttribute(
                        attributes[[name]],
                        dataset,
                        name,
                        asScalar = length(attributes[[name]]) == 1L
                    )
                }
            }
            replace_dataset()
        }
    }
    invisible(file)
}
