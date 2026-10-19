schema_fixture <- function(version) {
    file <- if (version == 1L) "sceptre-0.99.0-21f9ea0.rds" else "sceptre-schema2-4aad80b.rds"
    readRDS(test_path("fixtures", "saved-objects", file))
}

schema_physical_components <- function(object, version) {
    a <- attributes(object)
    if (version == 1L) stats::setNames(rep(list(a), 5L), c("data", "spec", "preparation", "fits", "results")) else
        lapply(a[c("data", "spec", "preparation", "fits", "results")], attributes)
}

test_that("real schema one and two fixtures migrate to one current authority", {
    for (version in c(1L, 2L)) {
        old <- schema_fixture(version)
        physical <- schema_physical_components(old, version)
        before <- serialize(old, NULL)
        migrated <- sceptre:::migrate_sceptre_object(old)
        expect_identical(migrated@schema_version, 3L)
        expect_setequal(slotNames(migrated), c("data", "spec", "preparation", "fits", "results", "schema_version"))
        expect_true(validObject(migrated))
        for (field in c("response_matrix", "grna_matrix")) {
            expect_identical(slot(migrated@data, field), physical$data[[field]][[1L]])
            expect_false(is.list(slot(migrated@data, field)))
        }
        aliases <- c(
            design_matrix = "covariate_matrix",
            nt_nonzero_matrix = "M_matrix",
            control_nonzero_totals = "n_nonzero_tot_vector",
            guide_assignments = "initial_grna_assignment_list",
            source_groups = "grna_assignments_raw",
            analysis_groups = "grna_assignments",
            guide_counts_per_cell = "grnas_per_cell",
            assignment_exclusions = "cells_w_zero_or_twoplus_grnas",
            analysis_cells = "cells_in_use")
        for (field in names(aliases)) {
            expect_identical(slot(migrated@preparation, field), physical$preparation[[aliases[[field]]]], info = field)
        }
        for (owner in c("data", "spec", "preparation")) {
            component <- slot(migrated, owner)
            fields <- setdiff(intersect(slotNames(component), names(physical[[owner]])),
                c("identity", "backend", "stages", "discovery_pairs", "positive_control_pairs", "grna_target_data_frame", "response_matrix", "grna_matrix"))
            for (field in fields) expect_identical(slot(component, field), physical[[owner]][[field]], info = paste(version, owner, field))
        }
        for (field in c("discovery_pairs", "positive_control_pairs")) {
            expanded <- physical$spec[[field]]
            expect_identical(slot(migrated@preparation, field), expanded)
            original <- physical$spec$requested[[field]]
            expected <- if (is.data.frame(original)) original else expanded[, setdiff(names(expanded), "grna_group"), drop = FALSE]
            expect_identical(slot(migrated@spec, field), expected)
        }
        expect_false("grna_group" %in% names(migrated@data@grna_target_data_frame))
        expect_identical(migrated@preparation@guide_groups, physical$preparation$grna_target_data_frame$grna_group)
        expect_length(migrated@fits@response_fits, 0L)
        expect_length(migrated@fits@response_identity, 0L)
        expect_length(migrated@fits@perturbation_fits, 0L)
        expect_length(migrated@fits@perturbation_identity, 0L)
        expect_identical(migrated@results@import_metadata$migrated_from_schema, version)
        expect_identical(migrated@results@import_metadata$provenance, "unknown")
        expect_false(any(c("B1", "B2", "B3", "requested", "negative_control_pairs") %in% slotNames(migrated@spec)))
        expect_false(any(c("covariate_names", "integer_id", "mitochondrial_gene") %in% slotNames(migrated@data)))
        expect_setequal(names(migrated@preparation@stages), c("design", "assignment", "cell_qc", "pair_qc"))
        for (stage in names(migrated@preparation@stages)) {
            expect_true(sceptre:::stage_complete(migrated, stage))
            expect_false("complete" %in% names(migrated@preparation@stages[[stage]]))
        }
        fields <- c(run_calibration_check = "calibration_result", run_power_check = "power_result", run_discovery_analysis = "discovery_result")
        for (analysis in names(fields)) {
            expected <- physical$results[[fields[[analysis]]]]
            record <- migrated@results@records[[analysis]]
            if (nrow(expected) || !is.null(physical$results$records[[analysis]])) {
                expect_identical(record$table, expected)
                expect_identical(record$invalidated_by, "legacy provenance unknown")
                expect_identical(record$provenance$certainty, "unknown")
                expect_false(any(c("complete", "status") %in% names(record)))
                expect_false(sceptre:::result_record_is_current(migrated, analysis))
            }
        }
        expect_identical(serialize(old, NULL), before)
        expect_identical(sceptre:::migrate_sceptre_object(migrated), migrated)
        path <- tempfile(fileext = ".rds")
        saveRDS(migrated, path)
        expect_identical(readRDS(path), migrated)
        unlink(path)
    }
})

test_that("legacy workflow flags cannot bless malformed preparation", {
    old <- schema_fixture(1L)
    design <- attr(old, "covariate_matrix")
    attr(old, "covariate_matrix") <- design[-1L, , drop = FALSE]
    expect_error(sceptre:::migrate_sceptre_object(old), "Design")
    old <- schema_fixture(1L)
    attr(old, "discovery_pairs_with_info") <- data.frame()
    migrated <- sceptre:::migrate_sceptre_object(old)
    expect_true(sceptre:::stage_complete(migrated, "cell_qc"))
    expect_false(sceptre:::stage_complete(migrated, "pair_qc"))
})

test_that("identity depends on actual count content and aligned metadata without RNG", {
    old <- readRDS(test_path("fixtures", "saved-objects", "sceptre-0.99.0-21f9ea0.rds"))
    object <- sceptre:::migrate_sceptre_object(old)
    set.seed(92)
    before <- .Random.seed
    same <- sceptre:::initialize_screen_identity(object)
    expect_identical(same@data@identity, object@data@identity)
    expect_identical(.Random.seed, before)
    changed <- object
    changed@data@response_matrix@x[1] <- changed@data@response_matrix@x[1] + 1
    changed <- sceptre:::initialize_screen_identity(changed)
    expect_false(identical(changed@data@identity, object@data@identity))
    changed <- object
    changed@data@covariate_data_frame <- changed@data@covariate_data_frame[rev(seq_len(nrow(changed@data@covariate_data_frame))), , drop = FALSE]
    changed <- sceptre:::initialize_screen_identity(changed)
    expect_false(identical(changed@data@identity, object@data@identity))
})

test_that("backed descriptors reopen relocated identical stores and reject changed descriptors", {
    skip_if_not_installed("ondisc")
    inputs <- readRDS(file.path(refactor_fixture_dir(), "backed-inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$high_perm
    directory <- tempfile("schema-backed-")
    relocated <- tempfile("schema-relocated-")
    on.exit(unlink(c(directory, relocated), recursive = TRUE), add = TRUE)
    object <- refactor_capture(refactor_import(inputs, case, TRUE, directory))$value
    sceptre::write_ondisc_backed_sceptre_object(object, directory)
    saved <- readRDS(file.path(directory, "sceptre_object.rds"))
    expect_s4_class(saved@data@response_matrix, "BackedCountStore")
    expect_s4_class(saved@data@grna_matrix, "BackedCountStore")
    expect_null(saved@data@response_matrix@handle)
    expect_null(saved@data@grna_matrix@handle)
    expect_match(saved@data@response_matrix@descriptor$content_md5, "^[a-f0-9]{32}$")
    expect_identical(dim(saved@data@response_matrix), dim(object@data@response_matrix))
    expect_identical(dimnames(saved@data@response_matrix), dimnames(object@data@response_matrix))
    expect_error(sceptre:::count_store_value(saved@data@response_matrix), "closed")
    dir.create(relocated)
    file.copy(file.path(directory, c("response.odm", "grna.odm")), relocated)
    opened <- sceptre::read_ondisc_backed_sceptre_object(
        file.path(directory, "sceptre_object.rds"),
        file.path(relocated, "response.odm"), file.path(relocated, "grna.odm"))
    expect_identical(opened@data@identity, object@data@identity)
    expect_identical(opened@fits, object@fits)
    expect_equal(sceptre::get_response_matrix(opened)[1, ], sceptre::get_response_matrix(object)[1, ])
    saved@data@response_matrix@descriptor$content_md5 <- strrep("0", 32)
    saveRDS(saved, file.path(directory, "tampered.rds"))
    expect_error(sceptre::read_ondisc_backed_sceptre_object(
        file.path(directory, "tampered.rds"),
        file.path(relocated, "response.odm"), file.path(relocated, "grna.odm")), "content does not match")
})

test_that("legacy result tables remain inspectable without being blessed as current", {
    old <- readRDS(test_path("fixtures", "saved-objects", "sceptre-0.99.0-21f9ea0.rds"))
    table <- data.frame(response_id = "saved_gene", p_value = 0.03)
    attr(old, "discovery_result") <- table
    flags <- attr(old, "functs_called")
    flags["run_discovery_analysis"] <- TRUE
    attr(old, "functs_called") <- flags
    migrated <- sceptre:::migrate_sceptre_object(old)
    record <- migrated@results@records$run_discovery_analysis
    expect_identical(record$invalidated_by, "legacy provenance unknown")
    expect_identical(record$provenance$certainty, "unknown")
    expect_false(sceptre:::result_record_is_current(migrated, "run_discovery_analysis"))
    expect_warning(saved <- sceptre::get_result(migrated, "run_discovery_analysis", historical = TRUE), "historical")
    expect_identical(attr(saved, "sceptre_result_record"), record)
    attr(saved, "sceptre_result_record") <- NULL
    expect_identical(saved, table)
    expect_error(sceptre::get_result(migrated, "run_discovery_analysis"), "current")
})


test_that("schema two result history survives without duplicate live tables or flags", {
    old <- schema_fixture(2L)
    results <- attr(old, "results")
    records <- attr(results, "records")
    previous <- data.frame(response_id = "previous_gene", p_value = 0.18)
    records$run_discovery_analysis$history <- list(list(result = previous,
        complete = TRUE, status = "current", scope = "complete_family",
        provenance = list(effective = list(B3 = 250L))))
    attr(results, "records") <- records
    attr(old, "results") <- results
    migrated <- sceptre:::migrate_sceptre_object(old)
    history <- migrated@results@records$run_discovery_analysis$history[[1L]]
    expect_identical(history$table, previous)
    expect_identical(history$provenance$saved_provenance$effective$B3, 250L)
    expect_identical(history$invalidated_by, "legacy provenance unknown")
    expect_false(any(c("complete", "status", "result") %in% names(history)))
    expect_setequal(slotNames(migrated@results), c("records", "import_metadata"))
    expect_setequal(names(attributes(migrated)), c("class", slotNames(migrated)))
})


test_that("legacy canonical recovery preserves repeated original requests", {
    one <- data.frame(grna_target = "A", response_id = "r1")
    original <- one[c(1L, 1L), , drop = FALSE]
    metadata <- data.frame(grna_target = c("A", "A"), grna_group = c("a1", "a2"))
    expanded <- data.frame(grna_target = rep("A", 4L), response_id = rep("r1", 4L),
        grna_group = rep(c("a1", "a2"), 2L))
    for (strategy in c("singleton", "bonferroni")) {
        recovered <- sceptre:::recover_legacy_requested_pairs(expanded, NULL, metadata, strategy)
        expect_identical(recovered$certainty, "reconstructed")
        expect_equal(unname(as.matrix(recovered$pairs)), unname(as.matrix(original)))
        expect_equal(nrow(recovered$pairs), 2L)
    }
    expect_identical(sceptre:::recover_legacy_requested_pairs(expanded, original, metadata,
        "singleton")$pairs, original)
    union <- original
    union$grna_group <- union$grna_target
    expect_identical(sceptre:::recover_legacy_requested_pairs(union, NULL, metadata, "union")$pairs, original)
    malformed <- expanded[c(1L, 3L, 2L, 4L), ]
    recovered <- sceptre:::recover_legacy_requested_pairs(malformed, NULL, metadata, "singleton")
    expect_identical(recovered$certainty, "unknown")
    expect_identical(recovered$pairs, malformed[, c("grna_target", "response_id")])
})


test_that("migration never infers original mitochondrial labels from row identifiers", {
    old <- schema_fixture(1L)
    covariates <- attr(old, "covariate_data_frame")
    covariates$response_p_mito <- rep(0.25, nrow(covariates))
    attr(old, "covariate_data_frame") <- covariates
    migrated <- sceptre:::migrate_sceptre_object(old)
    expect_length(migrated@data@response_names, 0L)
    expect_identical(migrated@data@covariate_data_frame$response_p_mito, covariates$response_p_mito)
    replacement <- sceptre::get_response_matrix(migrated)
    replacement@x[1L] <- replacement@x[1L] + 1
    expect_error(sceptre:::set_response_matrix(migrated, replacement),
        "Original mitochondrial feature labels are unknown; use import_data")

    old <- schema_fixture(2L)
    data <- attr(old, "data")
    known <- rep("ordinary", nrow(attr(data, "response_matrix")[[1L]]))
    known[1L] <- "MT-synthetic"
    attr(data, "response_names") <- known
    attr(old, "data") <- data
    migrated <- sceptre:::migrate_sceptre_object(old)
    expect_identical(migrated@data@response_names, known)
    expect_false(identical(known, rownames(migrated@data@response_matrix)))
})


test_that("count-store adapters are typed and validate descriptors without backend access", {
    descriptor <- list(path = "/nonexistent/closed.odm", content_md5 = strrep("a", 32),
        dimensions = c(2L, 3L), dimnames = list(c("one", "two"), NULL), integer_id = 42L)
    adapter <- sceptre:::as_count_store(NULL, descriptor)
    expect_s4_class(adapter, "BackedCountStore")
    expect_identical(sceptre:::count_store_descriptor(adapter), descriptor)
    expect_identical(dim(adapter), c(2L, 3L))
    expect_identical(dimnames(adapter), descriptor$dimnames)
    expect_true(validObject(adapter))
    expect_identical(sceptre:::close_count_store(adapter), adapter)
    expect_error(sceptre:::count_store_value(adapter), "closed")
    descriptor$dimensions <- c(3L, 3L)
    expect_error(sceptre:::as_count_store(NULL, descriptor), "labels must align")
    expect_error(sceptre:::as_count_store(list()), "Count stores must be")
    expect_identical(getSlots("ScreenData")[["response_matrix"]], "CountStore")
    expect_false("backend" %in% slotNames(new("ScreenData")))
})
