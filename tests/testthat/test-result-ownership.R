ownership_object <- function() {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_complement_perm
    object <- refactor_capture(refactor_prepare(inputs, case))$value
    refactor_capture(run_discovery_analysis(object, print_progress = FALSE))$value
}

# Run an actual := write in a data.table-aware caller environment.
ownership_add_column <- eval(quote(function(table) {
    table[, caller_added_column := 7L]
}), envir = asNamespace("data.table"))

test_that("public current results cannot mutate current or later historical records", {
    first <- ownership_object()
    later <- refactor_capture(run_discovery_analysis(first, print_progress = FALSE))$value
    original <- data.table::copy(first@results@records$run_discovery_analysis$table)
    public <- get_result(first, "run_discovery_analysis")
    expect_identical(class(public), class(original))
    expect_identical(names(attributes(public)),
        names(attributes(first@results@records$run_discovery_analysis$table)))
    expect_equal(public, original)
    data.table::set(public, i = 1L, j = "p_value", value = 0.123456)
    expect_no_warning(ownership_add_column(public))
    data.table::setattr(public, "caller_note", "edited")
    expect_equal(get_result(first, "run_discovery_analysis"), original)
    expect_equal(later@results@records$run_discovery_analysis$history[[1L]]$table, original)
    expect_true(refactor_private("result_record_is_current")(first, "run_discovery_analysis"))
    expect_true(refactor_private("result_record_is_current")(later, "run_discovery_analysis"))
})

test_that("admission owns caller tables, RunPlan pairs and nested attributes", {
    object <- ownership_object()
    proposal <- get_result(object, "run_discovery_analysis")
    pairs <- data.table::as.data.table(object@preparation@discovery_pairs_with_info)
    data.table::setattr(pairs, "details", list(table = data.table::data.table(value = 1L)))
    object@preparation@discovery_pairs_with_info <- pairs
    object <- refactor_private("commit_preparation_stage")(object, "pair_qc")
    plan <- refactor_private("make_run_plan")(object, "run_discovery_analysis")
    data.table::setattr(proposal, "details", list(table = data.table::data.table(value = 2L)))
    accepted <- refactor_private("commit_result_record")(object, "run_discovery_analysis",
        proposal, plan@pairs, 1L, plan)
    stored <- accepted@results@records$run_discovery_analysis
    certificate <- stored$identity
    expect_identical(refactor_private("as_run_plan")(stored$provenance$plan)@identity, certificate)
    data.table::set(proposal, j = "p_value", value = 0.321)
    data.table::set(attr(proposal, "details")$table, j = "value", value = 9L)
    data.table::set(plan@pairs, i = 1L, j = "response_id", value = "caller-edit")
    data.table::set(attr(plan@pairs, "details")$table, j = "value", value = 8L)
    expect_false(all(stored$table$p_value == 0.321, na.rm = TRUE))
    expect_identical(attr(stored$table, "details")$table$value, 2L)
    expect_false("caller-edit" %in% stored$provenance$plan@pairs$response_id)
    expect_identical(attr(stored$provenance$plan@pairs, "details")$table$value, 1L)
    expect_identical(refactor_private("as_run_plan")(stored$provenance$plan)@identity, certificate)
})

test_that("historical extraction owns every exposed record and preserves certificates", {
    first <- ownership_object()
    later <- refactor_capture(run_discovery_analysis(first, print_progress = FALSE))$value
    before <- later@results@records$run_discovery_analysis
    public <- get_result(later, "run_discovery_analysis", historical = TRUE)
    graph <- attr(public, "sceptre_result_record")
    expect_identical(refactor_private("as_run_plan")(graph$provenance$plan)@identity, before$identity)
    expect_identical(serialize(refactor_private("run_plan_payload")(graph$provenance$plan), NULL),
        serialize(refactor_private("run_plan_payload")(before$provenance$plan), NULL))
    data.table::set(graph$table, j = "p_value", value = 0.2)
    data.table::set(graph$history[[1L]]$table, j = "p_value", value = 0.3)
    data.table::set(graph$provenance$plan@pairs, j = "pass_qc", value = FALSE)
    data.table::setattr(graph$history[[1L]]$provenance$plan@pairs, "caller_attribute", 1L)
    data.table::setattr(graph$table, "nested", data.table::data.table(x = 1L))
    expect_equal(later@results@records$run_discovery_analysis, before)
    expect_equal(first@results@records$run_discovery_analysis$table, before$history[[1L]]$table)
    expect_true(refactor_private("result_record_is_current")(later, "run_discovery_analysis"))
})

test_that("owned snapshots retain plain-frame attributes and automatic row names", {
    source <- structure(data.frame(x = 1:3), detail = list(NULL,
        data.table::data.table(z = 1:2)))
    copied <- refactor_private("owned_result_value")(source)
    expect_identical(serialize(source, NULL), serialize(copied, NULL))
    data.table::set(attr(copied, "detail")[[2L]], j = "z", value = 3L)
    expect_identical(attr(source, "detail")[[2L]]$z, 1:2)
})

test_that("old saved list plans keep their original certificates through extraction", {
    old <- readRDS(test_path("fixtures", "saved-objects",
        "sceptre-schema3-list-plans-2a0938e.rds"))
    saved <- old@results@records$run_discovery_analysis$provenance$plan
    public <- get_result(old, "run_discovery_analysis", historical = TRUE)
    exposed <- attr(public, "sceptre_result_record")
    expect_identical(refactor_private("as_run_plan")(exposed$provenance$plan)@identity,
        saved$identity)
    expect_identical(serialize(refactor_private("run_plan_payload")(exposed$provenance$plan), NULL),
        serialize(refactor_private("run_plan_payload")(saved), NULL))
    data.table::set(exposed$provenance$plan@pairs, j = "pass_qc", value = FALSE)
    data.table::set(exposed$table, j = "p_value", value = 0.456)
    expect_identical(old@results@records$run_discovery_analysis$provenance$plan, saved)
    expect_no_error(get_result(old, "run_discovery_analysis"))
})
