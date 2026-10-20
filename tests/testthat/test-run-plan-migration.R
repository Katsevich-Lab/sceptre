# This fixture was produced before RunPlan existed, using the frozen installed
# schema-three package and two public discovery runs. The capture script asserted
# currentness before and after serialization and retained the first real run in
# history. The outer schema stays three because values and input hashes agree.
list_plan_fixture <- function() readRDS(test_path("fixtures", "saved-objects",
    "sceptre-schema3-list-plans-2a0938e.rds"))

test_that("schema three list plans normalize without changing scientific identity", {
    old <- list_plan_fixture()
    before <- serialize(old, NULL)
    old_record <- old@results@records$run_discovery_analysis
    expect_true(is.list(old_record$provenance$plan))
    expect_false(isS4(old_record$provenance$plan))
    expect_length(old_record$history, 1L)
    migrated <- sceptre:::migrate_sceptre_object(old)
    record <- migrated@results@records$run_discovery_analysis
    expect_identical(migrated@schema_version, 3L)
    expect_s4_class(record$provenance$plan, "RunPlan")
    expect_s4_class(record$history[[1L]]$provenance$plan, "RunPlan")
    expect_true(sceptre:::result_record_is_current(migrated, "run_discovery_analysis"))
    expect_identical(record$identity, old_record$identity)
    expect_identical(record$table, old_record$table)
    expect_identical(record$history[[1L]]$identity, old_record$history[[1L]]$identity)
    expect_identical(record$history[[1L]]$table, old_record$history[[1L]]$table)
    expect_identical(migrated@data, old@data)
    expect_identical(migrated@fits@response_fits, attr(attr(old, "fits"), "entries"))
    expect_identical(migrated@fits@response_identity, attr(attr(old, "fits"), "identity"))
    expect_length(migrated@fits@perturbation_fits, 0L)
    expect_length(migrated@fits@perturbation_identity, 0L)
    expect_identical(migrated@spec@perturbation_fit_method, character())
    expect_identical(get_perturbation_fit_method(migrated), "glm.fit")
    physical_spec <- migrated@spec
    attr(physical_spec, "perturbation_fit_method") <- NULL
    expect_identical(physical_spec, old@spec)
    expect_identical(record$provenance$spec, old_record$provenance$spec)
    expect_identical(migrated@preparation, old@preparation)
    expect_identical(serialize(old, NULL), before)
    expect_identical(sceptre:::migrate_sceptre_object(migrated), migrated)
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    saveRDS(migrated, path)
    restored <- sceptre:::migrate_sceptre_object(readRDS(path))
    expect_identical(restored, migrated)
    expect_true(sceptre:::result_record_is_current(restored, "run_discovery_analysis"))
})

test_that("normalizing list plans preserves prior invalidation and historical tables", {
    old <- list_plan_fixture()
    old@results@records$run_discovery_analysis$invalidated_by <- "original upstream change"
    original <- old@results@records$run_discovery_analysis
    migrated <- sceptre:::migrate_sceptre_object(old)
    record <- migrated@results@records$run_discovery_analysis
    expect_s4_class(record$provenance$plan, "RunPlan")
    expect_identical(record$invalidated_by, original$invalidated_by)
    expect_identical(record$identity, original$identity)
    expect_false(sceptre:::result_record_is_current(migrated, "run_discovery_analysis"))
    expect_warning(table <- sceptre::get_result(migrated, "run_discovery_analysis", historical = TRUE), "historical")
    attr(table, "sceptre_result_record") <- NULL
    expect_identical(table, original$table)
})

test_that("malformed saved plans stay readable and cannot become current", {
    old <- list_plan_fixture()
    bad <- old@results@records$run_discovery_analysis$provenance$plan
    bad$B1 <- NULL
    old@results@records$run_discovery_analysis$provenance$plan <- bad
    original <- old@results@records$run_discovery_analysis
    migrated <- sceptre:::migrate_sceptre_object(old)
    record <- migrated@results@records$run_discovery_analysis
    expect_null(record$provenance$plan)
    expect_identical(record$provenance$unrecognized_plan$payload, bad)
    expect_identical(record$invalidated_by, "unrecognized saved run plan")
    expect_identical(record$identity, original$identity)
    expect_identical(record$table, original$table)
    expect_false(sceptre:::result_record_is_current(migrated, "run_discovery_analysis"))
    expect_error(sceptre::get_result(migrated, "run_discovery_analysis"), "current")
    expect_warning(table <- sceptre::get_result(migrated, "run_discovery_analysis", historical = TRUE), "historical")
    attr(table, "sceptre_result_record") <- NULL
    expect_identical(table, original$table)
    expect_s4_class(record$history[[1L]]$provenance$plan, "RunPlan")
})

test_that("a rejected historical plan does not invalidate a genuine current run", {
    old <- list_plan_fixture()
    old@results@records$run_discovery_analysis$history[[1L]]$provenance$plan$B2 <- -1L
    old@results@records$run_discovery_analysis$history[[1L]]$invalidated_by <- "previous historical reason"
    previous <- old@results@records$run_discovery_analysis$history[[1L]]
    migrated <- sceptre:::migrate_sceptre_object(old)
    record <- migrated@results@records$run_discovery_analysis
    expect_true(sceptre:::result_record_is_current(migrated, "run_discovery_analysis"))
    expect_null(record$invalidated_by)
    historical <- record$history[[1L]]
    expect_null(historical$provenance$plan)
    expect_true("previous historical reason" %in% historical$invalidated_by)
    expect_true("unrecognized saved run plan" %in% historical$invalidated_by)
    expect_identical(historical$identity, previous$identity)
    expect_identical(historical$table, previous$table)
    expect_identical(sceptre:::migrate_sceptre_object(migrated), migrated)
})


test_that("a changed old list budget cannot reuse its original execution certificate", {
    old <- list_plan_fixture()
    record <- old@results@records$run_discovery_analysis
    original_identity <- record$identity
    original_certificate <- record$provenance$plan$identity
    tampered <- record$provenance$plan
    tampered$B3 <- tampered$B3 + 1L
    # This still has all expected fields and nonnegative integer budgets. Only
    # the preserved original certificate exposes the changed execution scope.
    expect_identical(names(tampered), names(record$provenance$plan))
    expect_type(tampered$B3, "integer")
    old@results@records$run_discovery_analysis$provenance$plan <- tampered
    migrated <- sceptre:::migrate_sceptre_object(old)
    rejected <- migrated@results@records$run_discovery_analysis
    expect_null(rejected$provenance$plan)
    expect_identical(rejected$provenance$unrecognized_plan$payload, tampered)
    expect_identical(rejected$provenance$unrecognized_plan$payload$identity, original_certificate)
    expect_identical(rejected$identity, original_identity)
    expect_identical(rejected$table, record$table)
    expect_match(rejected$provenance$unrecognized_plan$error, "frozen execution payload")
    expect_identical(rejected$invalidated_by, "unrecognized saved run plan")
    expect_false(sceptre:::result_record_is_current(migrated, "run_discovery_analysis"))
    expect_warning(table <- sceptre::get_result(migrated, "run_discovery_analysis", historical = TRUE), "historical")
    attr(table, "sceptre_result_record") <- NULL
    expect_identical(table, record$table)
})
