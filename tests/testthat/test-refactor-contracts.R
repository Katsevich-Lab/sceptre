test_that("reference provenance and every fixture digest are intact", {
    manifest <- readRDS(file.path(refactor_fixture_dir(), "manifest.rds"))
    expect_identical(manifest$source_commit, "21f9ea098b69c444a884a49bc254de711968e3b5")
    expect_true(manifest$reference_source_verified)
    expect_true(manifest$clean_native_rebuild)
    files <- file.path(refactor_fixture_dir(), names(manifest$artifacts_md5))
    expect_true(all(file.exists(files)))
    expect_identical(unname(tools::md5sum(files)), unname(manifest$artifacts_md5))
    # CI and differing numerical environments must supply a verified original
    # package run, never expectations regenerated from the candidate package.
    oracle <- refactor_runtime_oracle()
    if (!is.null(oracle)) expect_identical(oracle$manifest$source_commit, manifest$source_commit)
})

test_that("temporary runtime references reject unpinned stale and incomplete evidence", {
    oracle <- refactor_runtime_oracle()
    if (is.null(oracle)) {
        expect_true(refactor_frozen_environment())
    } else {
        reject <- function(manifest, message) {
            directory <- tempfile("refactor-invalid-oracle-")
            dir.create(directory)
            on.exit(unlink(directory, recursive = TRUE), add = TRUE)
            saveRDS(manifest, file.path(directory, "runtime-manifest.rds"))
            withr::with_envvar(c(SCEPTRE_REFACTOR_ORACLE = directory), {
                expect_error(refactor_runtime_oracle(), message)
            })
        }
        changed <- oracle$manifest
        changed$source_commit <- strrep("0", 40L)
        reject(changed, "unpinned")
        changed <- oracle$manifest
        changed$inputs_md5[[1L]] <- strrep("0", 32L)
        reject(changed, "immutable fixture inputs")
        changed <- oracle$manifest
        changed$harness_md5 <- strrep("0", 32L)
        reject(changed, "observation harness differs")
        changed <- oracle$manifest
        changed$source_md5[[1L]] <- strrep("0", 32L)
        reject(changed, "source differs")
        changed <- oracle$manifest
        changed$environment$toolchain$config[["CXXFLAGS"]] <- "changed-compiler-flags"
        reject(changed, "different numerical environment")
        changed <- oracle$manifest
        changed$environment$RNGkind <- "changed-rng-kind"
        reject(changed, "different numerical environment")
        reject(oracle$manifest, "output files are missing")
    }
})

test_that("cached runtime evidence rejects live RNG flags and debug-mode changes", {
    oracle <- refactor_runtime_oracle()
    if (is.null(oracle)) {
        expect_true(refactor_frozen_environment())
    } else {
        # Prime the real directory cache before each change. A directory-only
        # cache incorrectly accepts every one of these stale references.
        expect_identical(refactor_runtime_oracle()$manifest, oracle$manifest)
        original_kind <- RNGkind()
        kind <- if (RNGkind()[[1L]] == "L'Ecuyer-CMRG") "Mersenne-Twister" else "L'Ecuyer-CMRG"
        withr::with_preserve_seed({
            tryCatch(withr::with_seed(2718L, {
                expect_error(refactor_runtime_oracle(), "different numerical environment")
            }, .rng_kind = kind), finally = do.call(RNGkind, as.list(original_kind)))
        })
        expect_identical(refactor_runtime_oracle()$manifest, oracle$manifest)
        withr::with_envvar(c(PKG_CPPFLAGS = paste(Sys.getenv("PKG_CPPFLAGS"), "-DSCEPTRE_STALE_ORACLE=1")), {
            expect_error(refactor_runtime_oracle(), "different numerical environment")
        })
        expect_identical(refactor_runtime_oracle()$manifest, oracle$manifest)
        withr::with_envvar(c(SCEPTRE_REFACTOR_ORACLE_DEBUG = if (oracle$manifest$debug_compile) "false" else "true"), {
            expect_error(refactor_runtime_oracle(), "different numerical environment")
        })
        expect_identical(refactor_runtime_oracle()$manifest, oracle$manifest)
    }
})

test_that("toolchain caching observes edits to effective Makevars at the same path", {
    makevars <- tempfile("refactor-Makevars-")
    on.exit(unlink(makevars), add = TRUE)
    writeLines("CXXFLAGS += -DSCEPTRE_REFERENCE_FIRST=1", makevars)
    withr::with_envvar(c(R_MAKEVARS_USER = makevars), {
        first <- refactor_oracle_toolchain()
        writeLines("CXXFLAGS += -DSCEPTRE_REFERENCE_SECOND=1", makevars)
        second <- refactor_oracle_toolchain()
        expect_false(identical(first$makevars$user$md5, second$makevars$user$md5))
        expect_false(identical(first$config[["CXXFLAGS"]], second$config[["CXXFLAGS"]]))
    })
})

test_that("numerical environment diagnostics identify fields without exposing values", {
    expected <- list(toolchain = list(environment_flags = c(CXXFLAGS = "private-value")), RNGkind = "old")
    actual <- list(toolchain = list(environment_flags = c(CXXFLAGS = "different-value", CC = "new")), RNGkind = "new")
    expect_identical(refactor_oracle_environment_differences(expected, actual),
        c("environment.toolchain.environment_flags.CXXFLAGS", "environment.toolchain.environment_flags.CC", "environment.RNGkind"))
    expect_identical(refactor_oracle_environment_differences(expected, expected), character())
})

test_that("dynamic Makevars do not hide changes to referenced external flags", {
    makevars <- tempfile("refactor-dynamic-Makevars-")
    on.exit(unlink(makevars), add = TRUE)
    writeLines("CXXFLAGS += $(SCEPTRE_REFERENCE_EXTERNAL_FLAGS)", makevars)
    withr::with_envvar(c(R_MAKEVARS_USER = makevars, SCEPTRE_REFERENCE_EXTERNAL_FLAGS = "-DFIRST=1"), {
        first <- refactor_oracle_toolchain()
        withr::with_envvar(c(SCEPTRE_REFERENCE_EXTERNAL_FLAGS = "-DSECOND=1"), {
            second <- refactor_oracle_toolchain()
            expect_false(identical(first$config[["CXXFLAGS"]], second$config[["CXXFLAGS"]]))
        })
    })
})

test_that("developer selection invalidates compiler and runtime oracle caches", {
    # Model a stable clang dispatcher with two developer backends. This never
    # calls xcode-select or changes the host's selected developer directory.
    mocks <- new.env(parent = environment(refactor_oracle_environment))
    state <- new.env(parent = emptyenv())
    state$developer <- "first"
    state$probes <- 0L
    mocks$refactor_oracle_file <- function(path, hash = TRUE) {
        if (identical(path, "/var/db/xcode_select_link")) {
            return(list(path = paste0("/mock/developer/", state$developer), exists = TRUE,
                        size = 1, mtime = 1, ctime = 1))
        }
        list(path = path, exists = TRUE, size = 1, mtime = 1, ctime = 1)
    }
    mocks$refactor_oracle_makevars <- function() {
        list(user = list(path = "", exists = FALSE), site = list(path = "", exists = FALSE))
    }
    mocks$Sys.which <- function(x) stats::setNames(rep("/mock/clang-dispatcher", length(x)), x)
    mocks$system2 <- function(command, args, ...) {
        if (identical(args, c("CMD", "config", "--all"))) {
            return(c("CC = mockcc", "CXX = mockcxx", "CXX20 = mockcxx", "CPPFLAGS = -DFROZEN=1"))
        }
        stopifnot(identical(args, "--version"))
        state$probes <- state$probes + 1L
        paste("mock compiler backend", state$developer)
    }
    clone <- function(fun, environment) eval(call("function", formals(fun), body(fun)), environment)
    toolchain_cache <- new.env(parent = mocks)
    toolchain_cache$cached <- NULL
    toolchain_cache$binary_state <- clone(get("binary_state", environment(refactor_oracle_toolchain)), toolchain_cache)
    mocks$refactor_oracle_toolchain <- clone(refactor_oracle_toolchain, toolchain_cache)
    mocks$refactor_oracle_environment <- clone(refactor_oracle_environment, mocks)
    first <- mocks$refactor_oracle_environment()
    probes <- state$probes
    expect_gt(probes, 0L)
    expect_identical(mocks$refactor_oracle_environment(), first)
    expect_identical(state$probes, probes)
    withr::with_envvar(c(PATH = paste(Sys.getenv("PATH"), "unused-reference-path", sep = .Platform$path.sep)), {
        expect_identical(mocks$refactor_oracle_environment(), first)
        expect_gt(state$probes, probes)
    })

    # Seed an already verified cache entry, then put its old environment in the
    # on-disk manifest. A backend change must force readback and reject it.
    directory <- tempfile("refactor-developer-oracle-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    frozen <- readRDS(file.path(refactor_fixture_dir(), "manifest.rds"))
    manifest <- list(schema_version = 1L, type = "sceptre-same-environment-reference",
        source_commit = frozen$source_commit, source_verified_before_and_after = TRUE,
        clean_native_rebuild = TRUE, fresh_process = TRUE,
        debug_compile = identical(Sys.getenv("SCEPTRE_REFACTOR_ORACLE_DEBUG", "true"), "true"),
        source_files = frozen$source_files, environment = first,
        outputs_md5 = stats::setNames(rep("mock", length(refactor_oracle_output_files(!is.null(first$dependencies$ondisc)))),
                                    refactor_oracle_output_files(!is.null(first$dependencies$ondisc))))
    saveRDS(manifest, file.path(directory, "runtime-manifest.rds"))
    runtime_cache <- new.env(parent = mocks)
    runtime_cache$cached <- list(directory = normalizePath(directory, winslash = "/"), manifest = manifest)
    runtime <- clone(refactor_runtime_oracle, runtime_cache)
    withr::with_envvar(c(SCEPTRE_REFACTOR_ORACLE = directory), {
        expect_identical(runtime()$manifest, manifest)
        state$developer <- "second"
        second <- mocks$refactor_oracle_environment()
        expect_identical(first$toolchain$config, second$toolchain$config)
        expect_identical(first$toolchain$binaries, second$toolchain$binaries)
        expect_false(identical(first$toolchain$developer_selection, second$toolchain$developer_selection))
        expect_false(identical(first$toolchain$compilers, second$toolchain$compilers))
        expect_gt(state$probes, probes)
        expect_error(runtime(), "different numerical environment")
        probes <- state$probes
        expect_identical(mocks$refactor_oracle_environment(), second)
        expect_identical(state$probes, probes)
        # Environment overrides must also re-probe the otherwise unchanged
        # dispatcher, even when xcode-select's selection itself stays fixed.
        withr::with_envvar(c(DEVELOPER_DIR = "/mock/override", SDKROOT = "/mock/sdk"), {
            mocks$refactor_oracle_environment()
            expect_gt(state$probes, probes)
            expect_error(runtime(), "different numerical environment")
        })
    })
})

test_that("source byte verification identifies changed paths and MD5 values", {
    directory <- tempfile("refactor-source-bytes-")
    dir.create(file.path(directory, "R"), recursive = TRUE)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    file <- file.path(directory, "R", "test.R")
    base::writeLines("original <- TRUE", file, useBytes = TRUE)
    sources <- c("R/test.R" = "sha256-placeholder")
    expected <- refactor_oracle_hashes(stats::setNames(file, names(sources)))
    expect_true(refactor_oracle_verify_source(directory, sources, expected))
    base::writeLines("changed <- TRUE", file, useBytes = TRUE)
    actual <- unname(tools::md5sum(file))
    expect_error(refactor_oracle_verify_source(directory, sources, expected),
                 paste0("R/test.R: expected ", expected[[1L]], ", actual ", actual), fixed = TRUE)
})

test_that("stored-platform qualification checks its recorded compiler and RNG", {
    current <- refactor_oracle_environment()
    manifest <- list(R = current$R, platform = current$platform, external = current$external,
        RNGkind = current$RNGkind, dependencies = lapply(current$dependencies, function(version) list(version = version)))
    compiler <- unname(current$toolchain$config[["CXX"]])
    program <- scan(text = compiler, what = "", quiet = TRUE, comment.char = "")[[1L]]
    native <- list(compiler = compiler, compiler_version = current$toolchain$compilers[[program]])
    expect_true(refactor_frozen_environment(manifest, native))
    native$compiler_version <- "different recorded compiler"
    expect_false(refactor_frozen_environment(manifest, native))
    native$compiler_version <- current$toolchain$compilers[[program]]
    manifest$RNGkind <- "different recorded RNG"
    expect_false(refactor_frozen_environment(manifest, native))
})

test_that("refactor preserves the exact public API and existing S4 facade", {
    expected <- readRDS(file.path(refactor_fixture_dir(), "api.rds"))
    expect_identical(refactor_api(), expected)
})

test_that("supported scientific workflows reproduce the pinned baseline", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    cases <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))
    expected <- refactor_expected("workflows.rds")
    for (name in names(cases)) {
        actual <- refactor_run_case(inputs, cases[[name]])
        # Candidate and original run in the same numerical environment. The
        # declared tolerance and exact identities/decisions remain unchanged.
        expect_equal(actual, expected[[name]], tolerance = 1e-8, info = name)
        for (stage in c("prepared", "calibrated", "final")) {
            a <- actual$value[[stage]]
            b <- expected[[name]]$value[[stage]]
            for (field in c("cells_in_use", "grna_assignments", "B1", "B2", "B3", "n_ok_discovery_pairs")) {
                expect_identical(a[[field]], b[[field]], info = paste(name, stage, field))
            }
            for (table in c("discovery_pairs_with_info", "positive_control_pairs_with_info", "negative_control_pairs", "calibration_result", "power_result", "discovery_result")) {
                for (column in intersect(c("grna_target", "grna_group", "response_id", "pass_qc", "significant", "stage"), names(a[[table]]$columns))) {
                    expect_identical(a[[table]]$columns[[column]], b[[table]]$columns[[column]], info = paste(name, stage, table, column))
                }
            }
        }
        if (identical(name, "low_bonferroni")) {
            info <- actual$value$prepared$discovery_pairs_with_info$columns
            selected <- info$response_id == "r3" & info$grna_group %in% c("a1", "a2")
            expect_identical(sort(info$pass_qc[selected]), c(FALSE, TRUE))
        }
    }
})

test_that("ordered assignment remapping preserves source and fit views", {
    actual <- refactor_index_observations()
    expect_identical(actual, readRDS(file.path(refactor_fixture_dir(), "indexing.rds")))
    expect_identical(actual$treatment_source, c(4L, 5L))
    expect_identical(actual$nt_source, c(8L, 2L, 11L))
    expect_identical(actual$combined_source, c(4L, 5L, 8L, 2L, 11L))
    expect_true(actual$caller_unchanged)
})

test_that("repeated analyses and existing cache invalidation match baseline", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    cases <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))
    actual <- refactor_cache_observations(inputs, cases$low_complement_perm)
    expected <- refactor_expected("cache.rds")
    expect_equal(actual, expected, tolerance = 1e-8)
    expect_equal(actual$value$warm$discovery_result, actual$value$cold$discovery_result, tolerance = 1e-8)
    expect_equal(actual$value$cleared$discovery_result, actual$value$cold$discovery_result, tolerance = 1e-8)
    expect_identical(actual$value$retained_side_cache, actual$value$cold$response_precomputations)
    expect_length(actual$value$formula_cache, 0L)
    expect_length(actual$value$changed_qc_cache, 0L)
    expect_false(identical(actual$value$reassigned$grna_assignments, actual$value$cold$grna_assignments))
    expect_equal(actual$value$reassigned$discovery_result, actual$value$reassigned_cold$discovery_result, tolerance = 1e-8)
})

test_that("baseline saved objects remain analyzable", {
    object <- readRDS(file.path(refactor_fixture_dir(), "saved-object.rds"))
    result <- refactor_capture(sceptre::run_discovery_analysis(object, print_progress = FALSE, output_amount = 2L))
    actual <- list(value = refactor_observe(result$value), warnings = result$warnings)
    expect_equal(actual, refactor_expected("saved-object-result.rds"), tolerance = 1e-8)
})

test_that("rejected requests and degenerate baseline behavior remain explicit", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_complement_perm
    actual <- refactor_edge_observations(inputs, case)
    expect_equal(actual, refactor_expected("edges.rds"), tolerance = 1e-8)
    expect_identical(actual, readRDS(file.path(refactor_fixture_dir(), "edges.rds")))
    expect_match(actual$duplicate_requests$error, "duplicate rows", fixed = TRUE)
})

test_that("serial and two-worker execution preserve pair identity and results", {
    skip_on_os("windows")
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    cases <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))
    serial <- refactor_run_case(inputs, cases$high_complement_crt)
    parallel <- refactor_run_case(inputs, cases$high_complement_crt, parallel = TRUE)
    # Existing caches may have worker-specific insertion order; compare keyed fits.
    for (stage in c("calibrated", "final")) {
        for (out in c("serial", "parallel")) {
            value <- get(out)
            fits <- value$value[[stage]]$response_precomputations
            value$value[[stage]]$response_precomputations <- fits[sort(names(fits))]
            assign(out, value)
        }
    }
    expect_equal(parallel, serial, tolerance = 1e-8)
})

test_that("backed permutation workflows and independent reattachment match sparse", {
    if (identical(Sys.getenv("SCEPTRE_REQUIRE_ONDISC"), "true")) {
        expect_true(requireNamespace("ondisc", quietly = TRUE))
    }
    skip_if_not_installed("ondisc")
    inputs <- readRDS(file.path(refactor_fixture_dir(), "backed-inputs.rds"))
    cases <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))
    directory <- tempfile("refactor-ondisc-")
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    actual <- refactor_run_case(inputs, cases$high_perm, backed = TRUE, directory = directory)
    expected <- refactor_expected("backed-workflow.rds")
    if (!is.null(refactor_runtime_oracle())) {
        expect_equal(actual, refactor_expected("backed-observations.rds"), tolerance = 1e-8)
    }
    # Backed import generates backing-file IDs using R's RNG. Scientific parity
    # does not require that import to consume the same draws as sparse import.
    backend_metadata <- refactor_expected("backed-metadata.rds")
    if (identical(as.character(packageVersion("ondisc")), backend_metadata$ondisc_version)) {
        expect_identical(actual$warnings, backend_metadata$warnings)
        expect_identical(actual$value$rng_after, backend_metadata$rng_after)
    }
    actual$value$rng_after <- expected$value$rng_after <- NULL
    actual$warnings <- expected$warnings <- NULL
    expect_equal(actual, expected, tolerance = 1e-8)

    object <- refactor_capture(refactor_prepare(inputs, cases$high_perm, TRUE, directory))$value
    sceptre::write_ondisc_backed_sceptre_object(object, directory)
    rm(object)
    gc()
    reopened <- sceptre::read_ondisc_backed_sceptre_object(
        file.path(directory, "sceptre_object.rds"), file.path(directory, "response.odm"), file.path(directory, "grna.odm"))
    reopened_counts <- do.call(rbind, lapply(rownames(inputs$response_matrix), function(id) {
        refactor_private("load_row")(sceptre::get_response_matrix(reopened), id)
    }))
    expect_equal(unname(reopened_counts), unname(inputs$response_matrix))
    continued <- refactor_capture(sceptre::run_discovery_analysis(reopened, print_progress = FALSE, output_amount = 3L))
    discovery_only <- refactor_expected("backed-discovery.rds")
    expect_equal(refactor_observe(continued$value)$discovery_result, discovery_only$value$discovery_result, tolerance = 1e-8)

    other <- tempfile("refactor-ondisc-other-")
    on.exit(unlink(other, recursive = TRUE), add = TRUE)
    invisible(refactor_capture(refactor_import(inputs, cases$high_perm, TRUE, other)))
    expect_error(sceptre::read_ondisc_backed_sceptre_object(
        file.path(directory, "sceptre_object.rds"), file.path(other, "response.odm"), file.path(directory, "grna.odm")), "distinct IDs")
    expect_error(sceptre::read_ondisc_backed_sceptre_object(
        file.path(directory, "sceptre_object.rds"), file.path(directory, "response.odm"), file.path(other, "grna.odm")), "distinct IDs")
})
