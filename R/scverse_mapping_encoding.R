# MuData 1.16 omits the standard encoding attributes on these global mapping
# elements. Add only the metadata of our newly owned two-modality screen;
# references and in-place analysis updates never pass through this writer.
write_h5mu_mapping_encodings <- function(file) {
    groups <- c("/obsm", "/obsmap")
    arrays <- c(
        "/obsm/grna",
        "/obsm/response",
        "/obsmap/grna",
        "/obsmap/response"
    )
    paths <- c(groups, arrays)
    fid <- rhdf5::H5Fopen(file)
    objects <- list()
    on.exit(
        {
            for (object in rev(objects)) {
                rhdf5::H5Oclose(object)
            }
            rhdf5::H5Fclose(fid)
        },
        add = TRUE
    )

    entries <- rhdf5::h5ls(fid, recursive = 2L)
    entry_paths <- paste0(
        sub("/$", "", entries$group),
        "/",
        entries$name
    )
    check_array <- function(object, path) {
        datatype <- entries$dclass[match(path, entry_paths)]
        space <- rhdf5::H5Dget_space(object)
        on.exit(rhdf5::H5Sclose(space), add = TRUE)
        if (
            !datatype %in% c("INTEGER", "FLOAT", "ENUM") ||
                rhdf5::H5Sget_simple_extent_dims(space)$rank < 1L
        ) {
            stop(
                sprintf("Invalid H5MU mapping array at %s.", path),
                call. = FALSE
            )
        }
    }
    check_attribute <- function(object, path, name, expected) {
        if (!rhdf5::H5Aexists(object, name)) {
            return(FALSE)
        }
        attribute <- rhdf5::H5Aopen(object, name)
        space <- rhdf5::H5Aget_space(attribute)
        on.exit(
            {
                rhdf5::H5Sclose(space)
                rhdf5::H5Aclose(attribute)
            },
            add = TRUE
        )
        if (
            rhdf5::H5Sget_simple_extent_dims(space)$rank != 0L ||
                !identical(rhdf5::H5Aread(attribute), expected)
        ) {
            stop(
                sprintf("Contradictory H5MU %s metadata at %s.", name, path),
                call. = FALSE
            )
        }
        TRUE
    }

    # Validate all six elements before changing any of them. Unknown groups
    # and datasets are deliberately outside this format-specific repair.
    pending <- list()
    for (path in paths) {
        if (!rhdf5::H5Lexists(fid, path)) {
            stop(
                sprintf("Missing H5MU mapping element %s.", path),
                call. = FALSE
            )
        }
        object <- rhdf5::H5Oopen(fid, path)
        objects[[path]] <- object
        group <- path %in% groups
        expected_type <- if (group) "H5I_GROUP" else "H5I_DATASET"
        if (!identical(rhdf5::H5Iget_type(object), expected_type)) {
            stop(
                sprintf("Invalid H5MU mapping element %s.", path),
                call. = FALSE
            )
        }
        if (!group) {
            check_array(object, path)
        }
        expected <- list(
            "encoding-type" = if (group) "dict" else "array",
            "encoding-version" = if (group) "0.1.0" else "0.2.0"
        )
        present <- vapply(
            names(expected),
            function(name) {
                check_attribute(object, path, name, expected[[name]])
            },
            logical(1)
        )
        pending[[path]] <- expected[!present]
    }
    for (path in names(pending)) {
        for (name in names(pending[[path]])) {
            rhdf5::h5writeAttribute(
                pending[[path]][[name]],
                objects[[path]],
                name,
                asScalar = TRUE
            )
        }
    }
    invisible(file)
}
