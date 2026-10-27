# Frozen behavior suite: fixture files, digests and provenance.
#
# Sourced together with harness.R. Writes only to tempdir().

consistency_sha256 <- function(files) {
    if (exists("sha256sum", envir = asNamespace("tools"))) {
        return(unname(tools::sha256sum(files)))
    }
    vapply(
        files,
        function(file) digest::digest(file = file, algo = "sha256"),
        character(1),
        USE.NAMES = FALSE
    )
}

# Digests of the inputs that come from the package under test.
consistency_input_digests <- function() {
    example <- consistency_example_data()
    # Serialization format 2 without its 14-byte header (which records the
    # writing R version), so the digest depends on the object only.
    object_md5 <- function(x) {
        bytes <- serialize(x, NULL, xdr = TRUE, version = 2L)
        file <- tempfile("sceptre-consistency-input-")
        on.exit(unlink(file), add = TRUE)
        writeBin(bytes[-seq_len(14L)], file)
        unname(tools::md5sum(file))
    }
    extdata <- c(
        "highmoi_example/gem_group_1/barcodes.tsv.gz",
        "highmoi_example/gem_group_1/features.tsv.gz",
        "highmoi_example/gem_group_1/matrix.mtx.gz",
        "highmoi_example/gem_group_2/barcodes.tsv.gz",
        "highmoi_example/gem_group_2/features.tsv.gz",
        "highmoi_example/gem_group_2/matrix.mtx.gz",
        "parse_example/all_genes.csv",
        "parse_example/all_grnas.csv",
        "parse_example/gene_mat.mtx",
        "parse_example/grna_mat.mtx"
    )
    list(
        synthetic_md5 = object_md5(consistency_synthetic_inputs()),
        example_data_md5 = lapply(example, object_md5),
        extdata_sha256 = stats::setNames(
            as.list(consistency_sha256(vapply(
                extdata,
                function(f) consistency_extdata(f),
                character(1)
            ))),
            extdata
        )
    )
}
