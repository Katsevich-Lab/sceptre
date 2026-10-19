##########################
# READ AND WRITE FUNCTIONS
##########################
#' Write or read an `ondisc`-backed `sceptre_object`
#'
#' `write_ondisc_backed_sceptre_object()` and
#' `read_ondisc_backed_sceptre_object()` enable the writing and reading of
#' `ondisc`-backed `sceptre_object`s, respectively. First,
#' `write_ondisc_backed_sceptre_object()` writes an `ondisc`-backed
#' `sceptre_object` to disk, creating a file `sceptre_object.rds` in the
#' specified directory. Next, `read_ondisc_backed_sceptre_object()` reads and
#' initializes a `sceptre_object` from a `sceptre_object.rds` file,
#' `response.odm` file, and `grna.odm` file stored on disk.
#'
#' @details
#' The saved RDS contains anchored backing descriptors and closed handles, rather
#' than live process-local file handles. Reading it with `readRDS()` permits
#' inspection of saved metadata and historical results, but does not reopen counts.
#' Use `read_ondisc_backed_sceptre_object()` before continuing backed analyses or
#' exporting reports that need counts. The optional ondisc installation is checked
#' for the required API when backed counts are opened; in-memory analyses do not
#' load ondisc.
#'
#' The response and guide files must be the original matching pair. Reopening
#' checks their shared import ID, anchored file checksums, dimensions and feature
#' order. Identical files may be moved or copied: pass their new paths to the
#' reader. Verified relocation preserves compatible fits and saved result
#' provenance. An import ID alone is not a content fingerprint.
#'
#' Active backing files are immutable: do not replace or edit them while an object
#' is in use. Saving and reopening verify the existing content anchors; analysis
#' getters do not repeatedly hash files. If count content or its imported summaries
#' must change, create a new screen with `import_data()` instead of overwriting the
#' backing files. Old saved objects whose original feature labels are unknown may
#' likewise require reimport before count replacement.
#'
#' @param sceptre_object a `sceptre_object`
#' @param directory_to_write the directory in which to write the
#' `sceptre_object.rds` file
#' @param sceptre_object_fp file path to a `sceptre_object.rds` file
#' @param response_odm_file_fp file path to the backing `.odm` file for the
#' response modality
#' @param grna_odm_file_fp file path to the backing `.odm` file for the gRNA
#' modality
#'
#' @return `write_ondisc_backed_sceptre_object()` returns NULL, and
#' `read_ondisc_backed_sceptre_object()` returns an `ondisc`-backed
#' `sceptre_object`
#' @export
#' @examples
#' if (requireNamespace("ondisc", quietly = TRUE)) {
#'   data(lowmoi_example_data)
#'   # 1. create ondisc-backed sceptre_object
#'   sceptre_object <- import_data(
#'     response_matrix = lowmoi_example_data$response_matrix,
#'     grna_matrix = lowmoi_example_data$grna_matrix,
#'     grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
#'     extra_covariates = lowmoi_example_data$extra_covariates,
#'     moi = "low",
#'     use_ondisc = TRUE,
#'     directory_to_write = tempdir()
#'   )
#'
#'   # 2. write
#'   write_ondisc_backed_sceptre_object(
#'     sceptre_object = sceptre_object,
#'     directory_to_write = tempdir()
#'   )
#'
#'   # 3. read
#'   rm(sceptre_object)
#'   sceptre_object <- read_ondisc_backed_sceptre_object(
#'     sceptre_object_fp = paste0(tempdir(), "/sceptre_object.rds"),
#'     response_odm_file_fp = paste0(tempdir(), "/response.odm"),
#'     grna_odm_file_fp = paste0(tempdir(), "/grna.odm")
#'   )
#' }
read_ondisc_backed_sceptre_object <- function(
    sceptre_object_fp,
    response_odm_file_fp,
    grna_odm_file_fp
) {
    # read in objects
    sceptre_object_fp <- prepare_ondisc_file_path(
        sceptre_object_fp,
        "`sceptre_object_fp`"
    )
    response_odm_file_fp <- prepare_ondisc_file_path(
        response_odm_file_fp,
        "`response_odm_file_fp`"
    )
    grna_odm_file_fp <- prepare_ondisc_file_path(
        grna_odm_file_fp,
        "`grna_odm_file_fp`"
    )
    sceptre_object <- migrate_sceptre_object(readRDS(sceptre_object_fp))
    response_descriptor <- count_store_descriptor(
        sceptre_object@data@response_matrix
    )
    grna_descriptor <- count_store_descriptor(sceptre_object@data@grna_matrix)
    expected_id <- response_descriptor$integer_id
    if (is.null(expected_id)) {
        expected_id <- sceptre_object@results@import_metadata$legacy_backend_id
    }
    if (!length(expected_id)) {
        stop("Backing store has no saved import ID; reimport the screen.")
    }
    matrices <- storage_reopen_backed_pair(
        response_odm_file_fp,
        grna_odm_file_fp,
        expected_id
    )
    if (
        (is.null(response_descriptor) || is.null(grna_descriptor)) &&
            is.null(sceptre_object@results@import_metadata$legacy_backend_id)
    ) {
        stop(
            "Backing store has no anchored content descriptor; reimport the screen."
        )
    }
    # Verified relocation changes only the descriptor path and process handle.
    sceptre_object@data@response_matrix <- as_count_store(
        matrices$response_matrix,
        response_descriptor
    )
    sceptre_object@data@grna_matrix <- as_count_store(
        matrices$grna_matrix,
        grna_descriptor
    )
    if (is.null(response_descriptor) || is.null(grna_descriptor)) {
        sceptre_object <- initialize_screen_identity(sceptre_object)
        for (stage in names(sceptre_object@preparation@stages)) {
            sceptre_object@preparation@stages[[
                stage
            ]]$identity <- stage_input_identity(sceptre_object, stage)
        }
    }
    sceptre_object@results@import_metadata$legacy_backend_id <- NULL
    return(sceptre_object)
}


#' @rdname read_ondisc_backed_sceptre_object
#' @export
write_ondisc_backed_sceptre_object <- function(
    sceptre_object,
    directory_to_write
) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    validate_backing_for_save(sceptre_object)
    # Keep the imported descriptor identity. Only the location may change after
    # verified reattachment; never bless changed content by writing fresh hashes.
    # Persist descriptors, never process-local live HDF5 handles.
    sceptre_object@data@grna_matrix <- close_count_store(
        sceptre_object@data@grna_matrix
    )
    sceptre_object@data@response_matrix <- close_count_store(
        sceptre_object@data@response_matrix
    )
    # check that dir exists
    directory_to_write <- prepare_directory_to_write(directory_to_write)
    saveRDS(sceptre_object, file.path(directory_to_write, "sceptre_object.rds"))
}
