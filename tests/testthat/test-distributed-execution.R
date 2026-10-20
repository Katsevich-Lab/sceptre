distributed_fixture <- function(name, mechanism = NULL) {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))[[name]]
    if (!is.null(mechanism)) case$mechanism <- mechanism
    refactor_capture(refactor_prepare(inputs, case))$value
}

distributed_proposals <- function(object, manifest) lapply(manifest$units,
    function(unit) execute_work_unit(object, manifest, unit$id))

distributed_resign <- function(value) {
    value$identity <- refactor_private("execution_digest")(value[setdiff(names(value), "identity")])
    value
}

distributed_observe <- function(object, analysis) list(
    result = refactor_plain(get_result(object, analysis)), fits = object@fits,
    identity = object@results@records[[analysis]]$identity)

test_that("whole-key permutation and CRT units reproduce direct complete-family results", {
    cases <- list(c("low_complement_perm", NA), c("low_nt_perm", NA),
        c("high_complement_crt", NA), c("low_nt_perm", "crt"),
        c("low_bonferroni", "crt"), c("low_singleton", NA))
    for (case in cases) {
        object <- distributed_fixture(case[[1L]], if (is.na(case[[2L]])) NULL else case[[2L]])
        before <- serialize(object, NULL)
        set.seed(1938)
        rng <- .Random.seed
        manifest <- prepare_execution(object, "run_discovery_analysis", 3L, 3L)
        expect_identical(.Random.seed, rng)
        proposals <- distributed_proposals(object, manifest)
        expect_identical(.Random.seed, rng)
        actual <- collect_execution(object, manifest, rev(proposals))
        expect_identical(.Random.seed, rng)
        expected <- refactor_capture(run_discovery_analysis(object,
            print_progress = FALSE, output_amount = 3L))$value
        expect_identical(distributed_observe(actual, "run_discovery_analysis"),
            distributed_observe(expected, "run_discovery_analysis"), info = paste(case, collapse = "/"))
        expect_identical(serialize(object, NULL), before)
        expect_identical(sum(vapply(manifest$units, function(unit) length(unit$row_ids), integer(1))),
            sum(manifest$plan@pairs$pass_qc))
    }
})

test_that("calibration is sampled once globally and power finalizes once", {
    for (mechanism in c("permutations", "crt")) {
        object <- distributed_fixture("low_nt_perm", mechanism)
        set.seed(2738)
        manifest <- refactor_capture(prepare_execution(object, "run_calibration_check", 3L, 2L,
            n_calibration_pairs = 4L, calibration_group_size = 1L))$value
        after_preparation <- .Random.seed
        proposals <- distributed_proposals(object, manifest)
        actual <- collect_execution(object, manifest, proposals)
        expect_identical(.Random.seed, after_preparation)
        set.seed(2738)
        expected <- refactor_capture(run_calibration_check(object,
            n_calibration_pairs = 4L, calibration_group_size = 1L,
            print_progress = FALSE, output_amount = 2L))$value
        expect_identical(.Random.seed, after_preparation)
        expect_identical(distributed_observe(actual, "run_calibration_check"),
            distributed_observe(expected, "run_calibration_check"))
        manifest <- prepare_execution(object, "run_power_check", 3L, 2L)
        actual <- collect_execution(object, manifest, distributed_proposals(object, manifest))
        expected <- refactor_capture(run_power_check(object, print_progress = FALSE, output_amount = 2L))$value
        expect_identical(distributed_observe(actual, "run_power_check"),
            distributed_observe(expected, "run_power_check"))
    }
})

test_that("manifests and proposals roundtrip as values and isolate caller mutation", {
    object <- distributed_fixture("high_complement_crt")
    manifest <- prepare_execution(object, "run_discovery_analysis", 2L, 2L)
    proposals <- distributed_proposals(object, manifest)
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    saveRDS(list(manifest = manifest, proposals = rev(proposals)), path)
    restored <- readRDS(path)
    actual <- collect_execution(object, restored$manifest, restored$proposals)
    expected <- collect_execution(object, manifest, proposals)
    expect_identical(distributed_observe(actual, "run_discovery_analysis"),
        distributed_observe(expected, "run_discovery_analysis"))
    before <- serialize(object, NULL)
    data.table::set(proposals[[1L]]$result, j = "p_value", value = 0.5)
    expect_identical(serialize(object, NULL), before)
    expect_error(collect_execution(object, manifest, proposals), "certificate")
    expect_false(any(vapply(manifest, function(value) is.function(value) ||
        methods::is(value, "externalptr") || methods::is(value, "sceptre_object"), logical(1))))
})

test_that("collection rejects incomplete, duplicate, corrupt and foreign work atomically", {
    object <- distributed_fixture("high_complement_crt")
    manifest <- prepare_execution(object, "run_discovery_analysis", 2L, 2L)
    proposals <- distributed_proposals(object, manifest)
    before <- serialize(object, NULL)
    expect_error(collect_execution(object, manifest, proposals[-1L]), "every expected")
    expect_error(collect_execution(object, manifest, rep(proposals[1L], 2L)), "duplicate")
    bad <- proposals
    bad[[1L]]$result[[refactor_private("execution_row_column")()]] <-
        rep(bad[[1L]]$result[[refactor_private("execution_row_column")()]][[1L]], nrow(bad[[1L]]$result))
    bad[[1L]] <- distributed_resign(bad[[1L]])
    expect_error(collect_execution(object, manifest, bad), "global hypothesis rows")
    bad <- proposals
    bad[[1L]]$result$pass_qc[[1L]] <- FALSE
    bad[[1L]] <- distributed_resign(bad[[1L]])
    expect_error(collect_execution(object, manifest, bad), "QC metadata")
    bad <- proposals
    bad[[1L]]$runtime$R <- "foreign runtime"
    bad[[1L]] <- distributed_resign(bad[[1L]])
    expect_error(collect_execution(object, manifest, bad), "certificate")
    foreign <- refactor_capture(run_qc(object, response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1), additional_cells_to_remove = 1L))$value
    expect_error(collect_execution(foreign, manifest, proposals), "prepared screen context")
    expect_identical(serialize(object, NULL), before)
})

test_that("actual fit scopes and cross-unit response conflicts are checked", {
    object <- distributed_fixture("low_bonferroni", "crt")
    manifest <- prepare_execution(object, "run_discovery_analysis", 3L, 2L)
    proposals <- distributed_proposals(object, manifest)
    bad <- proposals
    group <- names(bad[[1L]]$perturbation_fits)[[1L]]
    bad[[1L]]$perturbation_fits[[group]]$identity$treatment <- "wrong membership"
    bad[[1L]] <- distributed_resign(bad[[1L]])
    expect_error(collect_execution(object, manifest, bad), "actual ordered inputs")
    object <- distributed_fixture("high_complement_crt")
    manifest <- prepare_execution(object, "run_discovery_analysis", 2L, 2L)
    proposals <- distributed_proposals(object, manifest)
    shared <- intersect(names(proposals[[1L]]$response_fits), names(proposals[[2L]]$response_fits))
    # Some families use disjoint responses. Add an eligible same-response request
    # on each guide to obtain a real duplicate fit when needed.
    if (!length(shared)) {
        inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
        pairs <- data.frame(grna_target = c("A", "B"), response_id = "r5")
        object <- refactor_capture(set_analysis_parameters(object, discovery_pairs = pairs,
            formula_object = ~ batch + z, resampling_mechanism = "crt"))$value
        object <- refactor_capture(assign_grnas(object, method = "thresholding", threshold = 40,
            print_progress = FALSE))$value
        object <- refactor_capture(run_qc(object, response_n_umis_range = c(0, 1),
            response_n_nonzero_range = c(0, 1)))$value
        manifest <- prepare_execution(object, "run_discovery_analysis", 2L, 2L)
        proposals <- distributed_proposals(object, manifest)
        shared <- intersect(names(proposals[[1L]]$response_fits), names(proposals[[2L]]$response_fits))
    }
    expect_true(length(shared) > 0L)
    bad <- proposals
    bad[[2L]]$response_fits[[shared[[1L]]]]$fitted_coefs[[1L]] <-
        bad[[2L]]$response_fits[[shared[[1L]]]]$fitted_coefs[[1L]] + 1
    bad[[2L]] <- distributed_resign(bad[[2L]])
    expect_error(collect_execution(object, manifest, bad), "Conflicting response fits")
})

test_that("worker failure returns no completed result or partial cache", {
    object <- distributed_fixture("low_complement_perm")
    local_mocked_bindings(fit_response = function(...) stop("distributed fit failed"))
    manifest <- prepare_execution(object, "run_discovery_analysis", 2L)
    before <- serialize(object, NULL)
    expect_error(execute_work_unit(object, manifest, manifest$units[[1L]]$id), "distributed fit failed")
    expect_identical(serialize(object, NULL), before)
    expect_false(refactor_private("result_record_is_current")(object, "run_discovery_analysis"))
})

test_that("a fresh R worker rebuilds full permutation resources from serializable inputs", {
    object <- distributed_fixture("low_nt_perm")
    manifest <- prepare_execution(object, "run_discovery_analysis", 2L, 3L)
    directory <- tempfile("distributed-process-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    paths <- file.path(directory, c("screen.rds", "execution.rds", "worker.R", "proposal.rds", "worker.log", "runtime.rds"))
    saveRDS(object, paths[[1L]])
    saveRDS(manifest, paths[[2L]])
    package_path <- getNamespaceInfo(asNamespace("sceptre"), "path")
    saveRDS(list(package_path = package_path, libraries = .libPaths(),
        source = file.exists(file.path(package_path, "R", "s4_class.R")),
        unit = manifest$units[[1L]]$id), paths[[6L]])
    writeLines(c(
        "args <- commandArgs(trailingOnly = TRUE)",
        "runtime <- readRDS(args[[4L]])",
        ".libPaths(runtime$libraries)",
        "if (runtime$source) pkgload::load_all(runtime$package_path, recompile = FALSE) else library(sceptre, lib.loc = dirname(runtime$package_path))",
        "screen <- readRDS(args[[1L]])",
        "execution <- readRDS(args[[2L]])",
        "proposal <- sceptre::execute_work_unit(screen, execution, runtime$unit)",
        "saveRDS(proposal, args[[3L]])"), paths[[3L]])
    status <- system2(file.path(R.home("bin"), "Rscript"),
        c("--vanilla", shQuote(paths[[3L]]), shQuote(paths[[1L]]), shQuote(paths[[2L]]),
            shQuote(paths[[4L]]), shQuote(paths[[6L]])), stdout = paths[[5L]], stderr = paths[[5L]])
    expect_identical(status, 0L, info = paste(readLines(paths[[5L]], warn = FALSE), collapse = "\n"))
    if (!identical(status, 0L)) return(invisible(NULL))
    remote <- readRDS(paths[[4L]])
    local <- execute_work_unit(object, manifest, manifest$units[[2L]]$id)
    actual <- collect_execution(object, manifest, list(local, remote))
    expected <- refactor_capture(run_discovery_analysis(object, print_progress = FALSE, output_amount = 3L))$value
    expect_identical(distributed_observe(actual, "run_discovery_analysis"),
        distributed_observe(expected, "run_discovery_analysis"))
})

test_that("collection owns raw tables and coefficient attributes retained by callers", {
    object <- distributed_fixture("low_complement_perm")
    manifest <- prepare_execution(object, "run_discovery_analysis", 2L, 2L)
    proposals <- distributed_proposals(object, manifest)
    id <- names(proposals[[1L]]$response_fits)[[1L]]
    coefficients <- proposals[[1L]]$response_fits[[id]]$fitted_coefs
    data.table::setattr(coefficients, "caller_details", list(table = data.table::data.table(value = 1L)))
    proposals[[1L]]$response_fits[[id]]$fitted_coefs <- coefficients
    proposals[[1L]] <- distributed_resign(proposals[[1L]])
    accepted <- collect_execution(object, manifest, proposals)
    original <- refactor_plain(get_result(accepted, "run_discovery_analysis"))
    data.table::set(attr(coefficients, "caller_details")$table, j = "value", value = 9L)
    data.table::setattr(coefficients, "caller_note", "changed")
    data.table::set(proposals[[1L]]$result, j = "p_value", value = 0.123)
    stored <- accepted@fits@response_fits[[id]]$fitted_coefs
    expect_identical(attr(stored, "caller_details")$table$value, 1L)
    expect_null(attr(stored, "caller_note"))
    expect_identical(refactor_plain(get_result(accepted, "run_discovery_analysis")), original)
})

test_that("row IDs require exact integers and structural field names are unique", {
    object <- distributed_fixture("low_complement_perm")
    manifest <- prepare_execution(object, "run_discovery_analysis", 2L, 2L)
    proposals <- distributed_proposals(object, manifest)
    bad <- proposals
    column <- refactor_private("execution_row_column")()
    bad[[1L]]$result[[column]] <- as.numeric(bad[[1L]]$result[[column]]) + 0.5
    bad[[1L]] <- distributed_resign(bad[[1L]])
    expect_error(collect_execution(object, manifest, bad), "global hypothesis rows")
    bad_manifest <- manifest
    bad_manifest <- c(bad_manifest, list(protocol = 1L))
    expect_error(execute_work_unit(object, bad_manifest, manifest$units[[1L]]$id), "manifest protocol")
    bad <- proposals
    bad[[1L]] <- c(bad[[1L]], list(protocol = 1L))
    expect_error(collect_execution(object, manifest, bad), "proposal certificate")
})

test_that("existing compatible caches are preserved across split runs", {
    object <- distributed_fixture("low_nt_perm")
    object <- refactor_capture(run_calibration_check(object, n_calibration_pairs = 4L,
        calibration_group_size = 1L, print_progress = FALSE))$value
    manifest <- prepare_execution(object, "run_discovery_analysis", 2L, 2L)
    actual <- collect_execution(object, manifest, distributed_proposals(object, manifest))
    expected <- refactor_capture(run_discovery_analysis(object, print_progress = FALSE, output_amount = 2L))$value
    expect_identical(distributed_observe(actual, "run_discovery_analysis"),
        distributed_observe(expected, "run_discovery_analysis"))
    expect_identical(actual@fits@response_fits, object@fits@response_fits)
    object <- distributed_fixture("high_complement_crt")
    object <- refactor_capture(run_discovery_analysis(object, print_progress = FALSE, output_amount = 2L))$value
    manifest <- prepare_execution(object, "run_discovery_analysis", 2L, 2L)
    actual <- collect_execution(object, manifest, distributed_proposals(object, manifest))
    expected <- refactor_capture(run_discovery_analysis(object, print_progress = FALSE, output_amount = 2L))$value
    expect_identical(distributed_observe(actual, "run_discovery_analysis"),
        distributed_observe(expected, "run_discovery_analysis"))
})

test_that("ambiguous numerical prepared metadata is rejected before worker execution", {
    object <- distributed_fixture("low_complement_perm")
    object@preparation@discovery_pairs_with_info$p_value <- 0.5
    object <- refactor_private("commit_preparation_stage")(object, "pair_qc")
    expect_error(prepare_execution(object, "run_discovery_analysis", 2L), "metadata collides")
})

test_that("backed permutation work reopens counts in a fresh process before collection", {
    skip_if_not_installed("ondisc")
    inputs <- readRDS(file.path(refactor_fixture_dir(), "backed-inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$high_perm
    directory <- tempfile("distributed-backed-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    object <- refactor_capture(refactor_prepare(inputs, case, TRUE, directory))$value
    manifest <- prepare_execution(object, "run_discovery_analysis", 2L, 2L)
    write_ondisc_backed_sceptre_object(object, directory)
    screen <- file.path(directory, "sceptre_object.rds")
    response <- file.path(directory, "response.odm")
    guide <- file.path(directory, "grna.odm")
    closed <- readRDS(screen)
    expect_null(closed@data@response_matrix@handle)
    reopened <- read_ondisc_backed_sceptre_object(screen, response, guide)
    execution_path <- file.path(directory, "execution.rds")
    proposal_path <- file.path(directory, "proposal.rds")
    runtime_path <- file.path(directory, "runtime.rds")
    worker_path <- file.path(directory, "worker.R")
    log_path <- file.path(directory, "worker.log")
    saveRDS(manifest, execution_path)
    package_path <- getNamespaceInfo(asNamespace("sceptre"), "path")
    saveRDS(list(package_path = package_path, libraries = .libPaths(),
        source = file.exists(file.path(package_path, "R", "s4_class.R")),
        unit = manifest$units[[1L]]$id, screen = screen, response = response, guide = guide), runtime_path)
    writeLines(c("args <- commandArgs(trailingOnly = TRUE)",
        "runtime <- readRDS(args[[1L]])", ".libPaths(runtime$libraries)",
        "if (runtime$source) pkgload::load_all(runtime$package_path, recompile = FALSE) else library(sceptre, lib.loc = dirname(runtime$package_path))",
        "screen <- sceptre::read_ondisc_backed_sceptre_object(runtime$screen, runtime$response, runtime$guide)",
        "proposal <- sceptre::execute_work_unit(screen, readRDS(args[[2L]]), runtime$unit)",
        "saveRDS(proposal, args[[3L]])"), worker_path)
    status <- system2(file.path(R.home("bin"), "Rscript"),
        c("--vanilla", shQuote(worker_path), shQuote(runtime_path), shQuote(execution_path), shQuote(proposal_path)),
        stdout = log_path, stderr = log_path)
    expect_identical(status, 0L, info = paste(readLines(log_path, warn = FALSE), collapse = "\n"))
    if (!identical(status, 0L)) return(invisible(NULL))
    remote <- readRDS(proposal_path)
    local <- execute_work_unit(reopened, manifest, manifest$units[[2L]]$id)
    actual <- collect_execution(reopened, manifest, list(local, remote))
    expected <- refactor_capture(run_discovery_analysis(reopened, print_progress = FALSE, output_amount = 2L))$value
    expect_identical(distributed_observe(actual, "run_discovery_analysis"),
        distributed_observe(expected, "run_discovery_analysis"))
    invalid <- reopened
    invalid@spec@run_permutations <- FALSE
    expect_error(prepare_execution(invalid, "run_discovery_analysis", 2L), "Backed CRT remains unsupported")
})


test_that("execution certificates bind actual R methods and loaded native bytes", {
    object <- distributed_fixture("low_complement_perm")
    certificate <- refactor_private("execution_build_certificate")()
    expect_identical(certificate, refactor_private("execution_build_certificate")())
    dll <- getLoadedDLLs()[["sceptre"]][["path"]]
    expect_identical(certificate$native, digest::digest(file = dll, algo = "sha256"))
    manifest <- prepare_execution(object, "run_discovery_analysis", 2L)
    version <- as.character(getNamespaceVersion(asNamespace("sceptre")))
    local_mocked_bindings(fit_response = function(...) stop("different implementation"))
    expect_identical(as.character(getNamespaceVersion(asNamespace("sceptre"))), version)
    expect_false(identical(certificate$R, refactor_private("execution_build_certificate")()$R))
    expect_error(execute_work_unit(object, manifest, manifest$units[[1L]]$id), "certificate is incompatible")
})

test_that("native build mismatch is rejected without re-signing an incoming manifest", {
    object <- distributed_fixture("low_complement_perm")
    manifest <- prepare_execution(object, "run_discovery_analysis", 2L)
    original <- serialize(manifest, NULL)
    local_mocked_bindings(execution_native_certificate = function() "different native bytes")
    expect_error(execute_work_unit(object, manifest, manifest$units[[1L]]$id), "certificate is incompatible")
    expect_identical(serialize(manifest, NULL), original)
})


test_that("redefining a dispatched package S4 method changes its build certificate", {
    object <- distributed_fixture("low_complement_perm")
    manifest <- prepare_execution(object, "run_discovery_analysis", 2L)
    original <- methods::selectMethod("dim", "BackedCountStore")
    on.exit(methods::setMethod("dim", "BackedCountStore", original), add = TRUE)
    methods::setMethod("dim", "BackedCountStore", function(x) c(1L, 2L))
    expect_false(identical(manifest$runtime$build, refactor_private("execution_build_certificate")()))
    expect_error(execute_work_unit(object, manifest, manifest$units[[1L]]$id), "certificate is incompatible")
    methods::setMethod("dim", "BackedCountStore", original)
    expect_identical(manifest$runtime$build, refactor_private("execution_build_certificate")())
})


test_that("valid warm perturbation fits cannot be changed or removed by proposals", {
    object <- distributed_fixture("high_complement_crt")
    object <- refactor_capture(run_discovery_analysis(object, print_progress = FALSE,
        output_amount = 2L))$value
    manifest <- prepare_execution(object, "run_discovery_analysis", 3L, 2L)
    proposals <- distributed_proposals(object, manifest)
    group <- names(proposals[[1L]]$perturbation_fits)[[1L]]
    expect_identical(proposals[[1L]]$perturbation_fits[[group]],
        object@fits@perturbation_fits[[group]])
    before <- serialize(object, NULL)
    for (change in c("coefficient", "diagnostics", "remove")) {
        bad <- proposals
        if (change == "coefficient") {
            bad[[1L]]$perturbation_fits[[group]]$fitted_coefs[[1L]] <-
                bad[[1L]]$perturbation_fits[[group]]$fitted_coefs[[1L]] + 1
        } else if (change == "diagnostics") {
            bad[[1L]]$perturbation_fits[[group]]$diagnostics$iterations <-
                bad[[1L]]$perturbation_fits[[group]]$diagnostics$iterations + 1L
        } else bad[[1L]]$perturbation_fits[group] <- list(NULL)
        bad[[1L]] <- distributed_resign(bad[[1L]])
        expect_error(collect_execution(object, manifest, bad), "Conflicting perturbation fit")
        expect_identical(serialize(object, NULL), before)
    }
    stale <- object
    stale@fits@perturbation_fits[[group]]$identity$treatment <- strrep("0", 64L)
    manifest <- prepare_execution(stale, "run_discovery_analysis", 3L, 2L)
    fresh <- collect_execution(stale, manifest, distributed_proposals(stale, manifest))
    expect_identical(fresh@fits@perturbation_fits[[group]], object@fits@perturbation_fits[[group]])
    direct <- refactor_capture(run_discovery_analysis(stale, print_progress = FALSE,
        output_amount = 2L))$value
    expect_identical(distributed_observe(fresh, "run_discovery_analysis"),
        distributed_observe(direct, "run_discovery_analysis"))
})

test_that("many CRT work units carry compact scope digests rather than cell vectors", {
    object <- distributed_fixture("high_complement_crt")
    manifest <- prepare_execution(object, "run_discovery_analysis", 100L, 2L)
    proposals <- distributed_proposals(object, manifest)
    expect_gt(length(proposals), 1L)
    context <- refactor_private("validate_execution_manifest")(object, manifest)$context
    expected <- list(response_identity = refactor_private("execution_digest")(context$fit_identity),
        perturbation_identity = refactor_private("execution_digest")(context$perturbation_identity))
    for (proposal in proposals) {
        headers <- proposal[names(expected)]
        expect_identical(headers, expected)
        expect_true(all(vapply(headers, function(value) is.character(value) &&
            length(value) == 1L && grepl("^[[:xdigit:]]{64}$", value), logical(1))))
        expect_lt(as.numeric(object.size(headers)), as.numeric(object.size(list(
            response_identity = context$fit_identity, perturbation_identity = context$perturbation_identity))))
    }
    collected <- collect_execution(object, manifest, proposals)
    expect_identical(collected@fits@response_identity, context$fit_identity)
    expect_identical(collected@fits@perturbation_identity, context$perturbation_identity)
})
