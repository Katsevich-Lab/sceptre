test_that("physical legacy slots migrate with data and unknown history preserved", {
    old <- readRDS(test_path("fixtures", "saved-objects", "sceptre-0.99.0-21f9ea0.rds"))
    physical <- attributes(old)
    expect_true("response_matrix" %in% names(physical))
    expect_false("data" %in% names(physical))
    before <- serialize(old, NULL)
    migrated <- sceptre:::migrate_sceptre_object(old)
    expect_identical(migrated@schema_version, 2L)
    owners <- sceptre:::component_slot_owners()
    for (field in setdiff(intersect(names(owners), names(physical)), "response_precomputations")) {
        expect_identical(slot(slot(migrated, owners[[field]]), field), physical[[field]], info = field)
    }
    expect_identical(migrated@spec@requested$source, "legacy")
    expect_identical(migrated@spec@requested$certainty, "unknown")
    expect_identical(migrated@spec@requested$effective_snapshot$B3, physical$B3)
    expect_false("B3" %in% names(migrated@spec@requested))
    expect_identical(migrated@results@records$legacy$provenance, "unknown")
    expect_identical(migrated@results@records$legacy$unverified_fits, physical$response_precomputations)
    expect_length(migrated@fits@response_precomputations, 0L)
    expect_false(any(migrated@functs_called[c("run_calibration_check", "run_power_check", "run_discovery_analysis")]))
    expect_identical(serialize(old, NULL), before)
    expect_identical(sceptre:::migrate_sceptre_object(migrated), migrated)
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    saveRDS(migrated, path)
    expect_identical(readRDS(path), migrated)
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
    changed@data@response_matrix[[1]]@x[1] <- changed@data@response_matrix[[1]]@x[1] + 1
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
    expect_length(saved@data@response_matrix, 0L)
    expect_length(saved@data@grna_matrix, 0L)
    expect_match(saved@data@backend$response$content_md5, "^[a-f0-9]{32}$")
    dir.create(relocated)
    file.copy(file.path(directory, c("response.odm", "grna.odm")), relocated)
    opened <- sceptre::read_ondisc_backed_sceptre_object(
        file.path(directory, "sceptre_object.rds"),
        file.path(relocated, "response.odm"), file.path(relocated, "grna.odm"))
    expect_identical(opened@data@identity, object@data@identity)
    expect_identical(opened@fits, object@fits)
    expect_equal(sceptre::get_response_matrix(opened)[1, ], sceptre::get_response_matrix(object)[1, ])
    saved@data@backend$response$content_md5 <- strrep("0", 32)
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
    expect_identical(record$status, "historical")
    expect_identical(record$provenance$certainty, "unknown")
    expect_false(sceptre:::result_record_is_current(migrated, "run_discovery_analysis"))
    expect_warning(saved <- sceptre::get_result(migrated, "run_discovery_analysis", historical = TRUE), "historical")
    expect_identical(attr(saved, "sceptre_result_record"), record)
    attr(saved, "sceptre_result_record") <- NULL
    expect_identical(saved, table)
    expect_error(sceptre::get_result(migrated, "run_discovery_analysis"), "current")
})
