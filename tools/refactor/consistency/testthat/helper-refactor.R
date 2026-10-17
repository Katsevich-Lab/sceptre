# Characterization helpers only; these are not package API or replacement code.
refactor_fixture_dir <- function() {
    getOption("sceptre.refactor.fixture_dir", testthat::test_path("fixtures", "refactor"))
}

refactor_private <- function(name) getFromNamespace(name, "sceptre")

# Runtime references never supply inputs or API/mapping expectations. They run
# the immutable inputs through the original package in this numerical setting.
refactor_oracle_file <- function(path, hash = TRUE) {
    exists <- length(path) == 1L && nzchar(path) && file.exists(path)
    if (!exists) return(list(path = path, exists = FALSE))
    path <- normalizePath(path.expand(path), winslash = "/", mustWork = TRUE)
    if (hash) return(list(path = path, exists = TRUE, md5 = unname(tools::md5sum(path))))
    info <- file.info(path)
    list(path = path, exists = TRUE, size = info$size,
         mtime = as.numeric(info$mtime), ctime = as.numeric(info$ctime))
}

refactor_oracle_makevars <- function() {
    personal <- path.expand("~/.R")
    fallback <- if (.Platform$OS.type == "windows") {
        file.path(personal, c("Makevars.ucrt", "Makevars.win64", "Makevars.win", "Makevars"))
    } else file.path(personal, c(paste0("Makevars-", R.version$platform), "Makevars"))
    candidates <- c(Sys.getenv("R_MAKEVARS_USER"), fallback)
    candidates <- candidates[nzchar(candidates) & file.exists(path.expand(candidates))]
    user <- if (length(candidates)) candidates[[1L]] else ""
    site <- Sys.getenv("R_MAKEVARS_SITE", unset = NA_character_)
    if (is.na(site)) site <- file.path(R.home("etc"), "Makevars.site")
    list(user = refactor_oracle_file(user), site = refactor_oracle_file(site))
}

refactor_oracle_toolchain <- local({
    cached <- NULL
    binary_state <- function(programs) lapply(programs, function(program) {
        found <- unname(Sys.which(program))
        refactor_oracle_file(if (nzchar(found)) found else program, hash = FALSE)
    })
    function() {
        variables <- c("CC", "CXX", "CXX11", "CXX14", "CXX17", "CXX20", "CXX23", "CXX26",
            "CFLAGS", "CXXFLAGS", "CXX11FLAGS", "CXX14FLAGS", "CXX17FLAGS", "CXX20FLAGS", "CXX23FLAGS", "CXX26FLAGS",
            "CPPFLAGS", "LDFLAGS", "PKG_CPPFLAGS", "PKG_CXXFLAGS", "PKG_LIBS", "SHLIB_CXXFLAGS", "SHLIB_CXXLDFLAGS",
            "MAKE", "MAKEFLAGS", "MFLAGS", "R_MAKEVARS_USER", "R_MAKEVARS_SITE", "R_ARCH", "R_PLATFORM",
            "PATH", "BINPREF", "RTOOLS40_HOME", "RTOOLS42_HOME", "R_CUSTOM_TOOLS_PATH", "R_CUSTOM_TOOLS_SOFT",
            "CPATH", "CPLUS_INCLUDE_PATH", "C_INCLUDE_PATH", "COMPILER_PATH", "GCC_EXEC_PREFIX", "LIBRARY_PATH",
            "SDKROOT", "DEVELOPER_DIR", "MACOSX_DEPLOYMENT_TARGET")
        present <- names(Sys.getenv())
        variables <- sort(unique(c(variables, present[grepl(
            "^(CXX|CC$|CPP|SHLIB|DYLIB|LTO|PKG_.*FLAGS|PKG_LIBS|CFLAGS|LDFLAGS|RTOOLS[0-9]+_HOME)", present
        )])), method = "radix")
        configuration_files <- c(file.path(R.home("etc"), c("Makeconf", "Renviron")),
            file.path(R.home("bin"), c("R", "config")),
            file.path(R.home("share"), "make", c("vars.mk", "shlib.mk")))
        key <- list(environment_flags = Sys.getenv(variables, unset = NA_character_),
            makevars = refactor_oracle_makevars(),
            configuration_files = lapply(configuration_files, refactor_oracle_file),
            developer_selection = refactor_oracle_file("/var/db/xcode_select_link", hash = FALSE))
        dynamic_makevars <- any(vapply(key$makevars, function(file) {
            file$exists && any(grepl("[$]|^[[:space:]]*(-?include|sinclude)[[:space:]]",
                                     readLines(file$path, warn = FALSE)))
        }, logical(1)))
        if (!dynamic_makevars && !is.null(cached) && identical(cached$key, key) &&
            identical(cached$binaries, binary_state(cached$programs))) return(cached$value)

        # One batched query includes effective user/site Makevars. Cache only
        # while its live input hashes, flags, PATH and compiler binaries agree.
        output <- system2(file.path(R.home("bin"), "R"), c("CMD", "config", "--all"),
                          stdout = TRUE, stderr = TRUE)
        if (!is.null(attr(output, "status"))) stop("Cannot fingerprint the active R compiler configuration.")
        output <- output[grepl("^[A-Z0-9_]+[[:space:]]*=", output)]
        keys <- sub("[[:space:]]*=.*$", "", output)
        values <- sub("^[^=]+=[[:space:]]*", "", output)
        config <- stats::setNames(values, keys)
        config <- config[grepl("^(CC|CXX|CPP|CFLAGS|CPIC|SHLIB|DYLIB|LDFLAGS|BLAS_LIBS|LAPACK_LIBS|LTO|MAKE|AR$|NM$|RANLIB$)", names(config))]
        commands <- config[grepl("^(CC|CXX[0-9]*|SHLIB_CXX[0-9]*LD)$", names(config))]
        commands <- commands[nzchar(commands)]
        command_tokens <- lapply(commands, function(command) {
            scan(text = command, what = "", quiet = TRUE, comment.char = "")
        })
        tokens <- unique(unlist(command_tokens))
        programs <- unique(c(vapply(command_tokens, `[[`, character(1), 1L),
                             tokens[!startsWith(tokens, "-") & nzchar(Sys.which(tokens))]))
        binaries <- binary_state(programs)
        versions <- if (!is.null(cached) && identical(cached$key, key) &&
                        identical(cached$programs, programs) &&
                        identical(cached$binaries, binaries)) cached$value$compilers else {
            stats::setNames(lapply(programs, function(program) {
                version <- system2(program, "--version", stdout = TRUE, stderr = TRUE)
                if (!is.null(attr(version, "status"))) stop("Cannot fingerprint compiler version: ", program)
                paste(version, collapse = " ")
            }), programs)
        }
        # PATH selects the queried binaries and invalidates this query cache.
        # R CMD check may prepend its own directories without changing any
        # selected compiler. Compare resolved binaries/configuration, not that
        # incidental transport path, when qualifying numerical equivalence.
        numerical_flags <- key$environment_flags[names(key$environment_flags) != "PATH"]
        value <- list(config = config, compilers = versions, binaries = binaries,
                      makevars = key$makevars, environment_flags = numerical_flags,
                      configuration_files = key$configuration_files,
                      developer_selection = key$developer_selection)
        cached <<- list(key = key, programs = programs, binaries = binaries, value = value)
        value
    }
})

refactor_oracle_environment <- function() {
    packages <- c("BH", "data.table", "dplyr", "Matrix", "Rcpp", "withr", "ondisc")
    list(
        R = R.version.string, platform = R.version$platform,
        toolchain = refactor_oracle_toolchain(), RNGkind = RNGkind(),
        external = extSoftVersion(), lapack = La_version(),
        dependencies = stats::setNames(lapply(packages, function(package) {
            if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else NULL
        }), packages),
        threads = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))
    )
}

# Report field names only; environment values can contain private local paths.
refactor_oracle_environment_differences <- function(expected, actual, path = "environment") {
    if (identical(expected, actual)) return(character())
    named <- function(x) is.list(x) || (!is.null(names(x)) && is.atomic(x))
    if (named(expected) && named(actual) &&
        !is.null(names(expected)) && !is.null(names(actual))) {
        keys <- union(names(expected), names(actual))
        field <- function(x, key) if (key %in% names(x)) x[[key]] else NULL
        return(unlist(lapply(keys, function(key) {
            refactor_oracle_environment_differences(field(expected, key), field(actual, key), paste0(path, ".", key))
        }), use.names = FALSE))
    }
    path
}

refactor_frozen_environment <- function(
    manifest = readRDS(file.path(refactor_fixture_dir(), "manifest.rds")),
    native = readRDS(file.path(refactor_fixture_dir(), "native-expected.rds"))$provenance$platform
) {
    current <- refactor_oracle_environment()
    command <- unname(current$toolchain$config[["CXX"]])
    program <- scan(text = command, what = "", quiet = TRUE, comment.char = "")[[1L]]
    compiler_matches <- identical(command, native$compiler) &&
        identical(current$toolchain$compilers[[program]], native$compiler_version)
    packages <- setdiff(names(current$dependencies), "ondisc")
    identical(current$R, manifest$R) && identical(current$platform, manifest$platform) &&
        identical(current$external, manifest$external) && compiler_matches &&
        identical(current$RNGkind, manifest$RNGkind) &&
        all(vapply(packages, function(package) {
            identical(current$dependencies[[package]], manifest$dependencies[[package]]$version)
        }, logical(1)))
}

refactor_oracle_input_files <- function() {
    c("inputs.rds", "cases.rds", "backed-inputs.rds", "saved-object.rds")
}

refactor_oracle_output_files <- function(backed) {
    c("workflows.rds", "cache.rds", "edges.rds", "saved-object-result.rds",
      "backed-workflow.rds", "backed-discovery.rds",
      if (backed) c("backed-metadata.rds", "backed-observations.rds"))
}

refactor_oracle_hashes <- function(files) {
    stats::setNames(unname(tools::md5sum(files)), names(files))
}

refactor_oracle_verify_source <- function(reference, source_files, source_md5) {
    files <- file.path(reference, names(source_files))
    if (!all(file.exists(files))) stop("Pinned runtime reference source files are missing.")
    actual <- refactor_oracle_hashes(stats::setNames(files, names(source_files)))
    if (!identical(names(source_md5), names(source_files)) || !identical(actual, source_md5)) {
        expected <- source_md5[names(actual)]
        mismatched <- names(actual)[is.na(expected) | actual != expected]
        detail <- vapply(mismatched, function(file) {
            paste0(file, ": expected ", expected[[file]], ", actual ", actual[[file]])
        }, character(1))
        stop("Runtime reference source differs from the immutable baseline source hashes (MD5):\n",
             paste(detail, collapse = "\n"),
             if (!identical(names(source_md5), names(source_files))) "\nSource manifest names/order differ.")
    }
    for (directory in c("R", "src", "data")) {
        actual <- list.files(file.path(reference, directory), recursive = TRUE,
                             all.files = TRUE, no.. = TRUE)
        actual <- actual[!grepl("\\.(o|so|dll|dylib|d)$", actual)]
        expected <- substring(names(source_files)[startsWith(names(source_files), paste0(directory, "/"))],
                              nchar(directory) + 2L)
        if (!setequal(actual, expected)) stop("Unexpected runtime reference source files under ", directory)
    }
    invisible(TRUE)
}

refactor_runtime_oracle <- local({
    cached <- NULL
    function() {
        directory <- Sys.getenv("SCEPTRE_REFACTOR_ORACLE")
        if (!nzchar(directory)) {
            if (identical(Sys.getenv("CI"), "true") || !refactor_frozen_environment()) {
                stop("A verified same-environment pinned reference is required. Run tools/refactor/generate-runtime-oracle.R and set SCEPTRE_REFACTOR_ORACLE.")
            }
            return(NULL)
        }
        directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
        current <- refactor_oracle_environment()
        debug <- identical(Sys.getenv("SCEPTRE_REFACTOR_ORACLE_DEBUG", "true"), "true")
        if (!is.null(cached) && identical(cached$directory, directory) &&
            identical(cached$manifest$environment, current) &&
            identical(cached$manifest$debug_compile, debug)) return(cached)
        manifest_path <- file.path(directory, "runtime-manifest.rds")
        if (!file.exists(manifest_path)) stop("Runtime oracle manifest is missing.")
        oracle <- readRDS(manifest_path)
        frozen <- readRDS(file.path(refactor_fixture_dir(), "manifest.rds"))
        required <- refactor_oracle_output_files(!is.null(current$dependencies$ondisc))
        if (!identical(oracle$schema_version, 1L) ||
            !identical(oracle$type, "sceptre-same-environment-reference") ||
            !identical(oracle$source_commit, "21f9ea098b69c444a884a49bc254de711968e3b5") ||
            !isTRUE(oracle$source_verified_before_and_after) || !isTRUE(oracle$clean_native_rebuild) ||
            !isTRUE(oracle$fresh_process) ||
            !identical(oracle$debug_compile, debug) ||
            !identical(oracle$source_files, frozen$source_files) ||
            !identical(oracle$environment, current) ||
            !identical(sort(names(oracle$outputs_md5)), sort(required))) {
            fields <- refactor_oracle_environment_differences(oracle$environment, current)
            stop("Runtime oracle is unpinned, incomplete, or from a different numerical environment.",
                 if (length(fields)) paste0(" Differing fields: ", paste(fields, collapse = ", ")) else "")
        }
        inputs <- refactor_oracle_input_files()
        expected_inputs <- frozen$artifacts_md5[inputs]
        actual_inputs <- refactor_oracle_hashes(stats::setNames(file.path(refactor_fixture_dir(), inputs), inputs))
        if (!identical(oracle$inputs_md5, expected_inputs) || !identical(actual_inputs, expected_inputs)) {
            stop("Runtime oracle inputs differ from immutable fixture inputs.")
        }
        helper <- file.path(refactor_fixture_dir(), "..", "..", "helper-refactor.R")
        if (!identical(oracle$harness_md5, unname(refactor_oracle_hashes(helper)))) {
            stop("Runtime oracle observation harness differs from the tested harness.")
        }
        refactor_oracle_verify_source(oracle$reference_source, oracle$source_files, oracle$source_md5)
        outputs <- stats::setNames(file.path(directory, required), required)
        if (!all(file.exists(outputs)) || !identical(refactor_oracle_hashes(outputs), oracle$outputs_md5[required])) {
            stop("Runtime oracle output files are missing or their hashes differ.")
        }
        cached <<- list(directory = directory, manifest = oracle)
        cached
    }
})

refactor_expected <- function(file) {
    oracle <- refactor_runtime_oracle()
    frozen <- file.path(refactor_fixture_dir(), file)
    if (is.null(oracle)) return(readRDS(frozen))
    if (!file %in% names(oracle$manifest$outputs_md5)) stop("File is not a runtime scientific expectation: ", file)
    expected <- readRDS(file.path(oracle$directory, file))
    # Retain the stored Mac scientific checks when its recorded numerical
    # setting and debug compile mode apply. Other platforms still require the
    # original source's exact decisions and the unchanged floating tolerance.
    backed_metadata_matches <- file != "backed-metadata.rds" || identical(
        refactor_oracle_environment()$dependencies$ondisc,
        readRDS(file.path(refactor_fixture_dir(), "manifest.rds"))$dependencies$ondisc$version
    )
    if (refactor_frozen_environment() && isTRUE(oracle$manifest$debug_compile) && backed_metadata_matches &&
        file != "backed-observations.rds") {
        testthat::expect_equal(expected, readRDS(frozen), tolerance = 1e-8,
                               info = paste("stored platform reference", file))
    }
    expected
}

refactor_api <- function() {
    ns <- asNamespace("sceptre")
    # R/Matrix synthesize .__T__ method-table bindings; these are not user API.
    exports <- getNamespaceExports(ns)
    exports <- sort(exports[!startsWith(exports, ".__")], method = "radix")
    describe_formals <- function(f) {
        if (!is.function(f)) return(NULL)
        lapply(formals(f), function(x) paste(deparse(x, width.cutoff = 500L), collapse = "\n"))
    }
    method_names <- c("plot", "print")
    exported_methods <- lapply(method_names, function(name) {
        found <- methods::findMethods(name, where = ns)
        lapply(found, describe_formals)
    })
    names(exported_methods) <- method_names
    list(
        exports = exports,
        formals = stats::setNames(lapply(exports, function(x) describe_formals(get(x, ns))), exports),
        methods = exported_methods,
        slots = methods::getSlots(methods::getClass("sceptre_object", where = ns))
    )
}

refactor_make_inputs <- function() {
    withr::with_seed(20261016, {
        n <- 160L
        guide_ids <- c("a1", "a2", "b1", "b2", "nt1", "nt2", "nt3", "nt4")
        targets <- c("A", "A", "B", "B", rep("non-targeting", 4L))
        guide <- matrix(0, length(guide_ids), n,
                        dimnames = list(guide_ids, paste0("cell_", seq_len(n))))
        assigned <- rep(seq_along(guide_ids), length.out = n)
        guide[cbind(assigned, seq_len(n))] <- 50
        guide[1L, which(assigned == 1L)[seq_len(4L)]] <- 43
        response <- matrix(stats::rpois(8L * n, 3), 8L, n,
                           dimnames = list(paste0("r", seq_len(8L)), colnames(guide)))
        response[1L, assigned %in% 1:2] <- response[1L, assigned %in% 1:2] + 20
        response[2L, assigned %in% 3:4] <- 0
        response[3L, ] <- 0
        response[3L, assigned == 1L] <- 5
        response[3L, assigned >= 5L] <- 3
        response[4L, ] <- 0
        high <- guide
        high[1L, seq(10L, n, 16L)] <- 50
        high[3L, seq(5L, n, 16L)] <- 50
        high[6L, seq(3L, n, 16L)] <- 50
        list(
            response_matrix = response,
            low_grna_matrix = guide,
            high_grna_matrix = high,
            grna_target_data_frame = data.frame(grna_id = guide_ids, grna_target = targets),
            extra_covariates = data.frame(
                batch = factor(rep(c("b", "a", "b", "b", "a"), length.out = n)),
                z = sin(seq_len(n)), row.names = colnames(guide)
            ),
            discovery_pairs = data.frame(
                grna_target = c("B", "A", "A", "B", "A"),
                response_id = c("r5", "r3", "r1", "r2", "r4")
            ),
            positive_control_pairs = data.frame(grna_target = c("A", "B"), response_id = c("r1", "r2")),
            remove_cells = c(2L, 17L, 49L, 103L, 158L)
        )
    })
}

refactor_case_definitions <- function() {
    list(
        low_nt_perm = list(moi = "low", control = "nt_cells", mechanism = "permutations", integration = "union", side = "both", approximation = "skew_normal", output = 3L),
        low_complement_perm = list(moi = "low", control = "complement", mechanism = "permutations", integration = "union", side = "left", approximation = "no_approximation", output = 2L),
        high_complement_crt = list(moi = "high", control = "complement", mechanism = "crt", integration = "union", side = "right", approximation = "skew_normal", output = 3L),
        low_complement_crt = list(moi = "low", control = "complement", mechanism = "crt", integration = "union", side = "both", approximation = "no_approximation", output = 1L),
        high_perm = list(moi = "high", control = "complement", mechanism = "permutations", integration = "union", side = "both", approximation = "no_approximation", output = 3L),
        low_bonferroni = list(moi = "low", control = "nt_cells", mechanism = "permutations", integration = "bonferroni", side = "left", approximation = "skew_normal", output = 2L),
        low_singleton = list(moi = "low", control = "complement", mechanism = "permutations", integration = "singleton", side = "right", approximation = "no_approximation", output = 1L)
    )
}

refactor_expand_inputs <- function(inputs, copies = 8L) {
    n <- ncol(inputs$response_matrix)
    ids <- paste0("cell_", seq_len(n * copies))
    for (name in c("response_matrix", "low_grna_matrix", "high_grna_matrix")) {
        inputs[[name]] <- inputs[[name]][, rep(seq_len(n), copies), drop = FALSE]
        colnames(inputs[[name]]) <- ids
    }
    inputs$extra_covariates <- inputs$extra_covariates[rep(seq_len(n), copies), , drop = FALSE]
    rownames(inputs$extra_covariates) <- ids
    inputs$remove_cells <- as.integer(unlist(lapply(seq_len(copies) - 1L, function(i) inputs$remove_cells + i * n)))
    inputs
}

refactor_import <- function(inputs, case, backed = FALSE, directory = NULL) {
    sceptre::import_data(
        response_matrix = inputs$response_matrix,
        grna_matrix = inputs[[paste0(case$moi, "_grna_matrix")]],
        grna_target_data_frame = inputs$grna_target_data_frame,
        extra_covariates = inputs$extra_covariates, moi = case$moi,
        use_ondisc = backed, directory_to_write = directory
    )
}

refactor_prepare <- function(inputs, case, backed = FALSE, directory = NULL, assignment_threshold = 40) {
    object <- refactor_import(inputs, case, backed, directory)
    object <- sceptre::set_analysis_parameters(
        object, discovery_pairs = inputs$discovery_pairs,
        positive_control_pairs = inputs$positive_control_pairs,
        formula_object = stats::as.formula("~ batch + z", env = baseenv()),
        side = case$side, grna_integration_strategy = case$integration,
        control_group = case$control, resampling_mechanism = case$mechanism,
        resampling_approximation = case$approximation
    )
    object <- sceptre::assign_grnas(object, method = "thresholding", threshold = assignment_threshold, print_progress = FALSE)
    sceptre::run_qc(object,
        response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
        n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L,
        additional_cells_to_remove = inputs$remove_cells
    )
}

# Preserve table types/order without serializing data.table's external selfref.
refactor_plain <- function(x) {
    if (is.data.frame(x)) {
        return(list(class = class(x), names = names(x), row_names = rownames(x),
                    columns = lapply(x, refactor_plain)))
    }
    if (is.list(x)) return(lapply(x, refactor_plain))
    x
}

refactor_observe <- function(object) {
    selected <- c(
        "cells_in_use", "grna_assignments", "cell_removal_metrics", "M_matrix",
        "n_nonzero_tot_vector", "discovery_pairs_with_info", "positive_control_pairs_with_info",
        "negative_control_pairs", "n_ok_discovery_pairs", "n_ok_positive_control_pairs",
        "B1", "B2", "B3", "side_code", "functs_called", "last_function_called",
        "response_precomputations", "calibration_result", "power_result", "discovery_result"
    )
    stats::setNames(lapply(selected, function(name) refactor_plain(methods::slot(object, name))), selected)
}

refactor_capture <- function(code) {
    warnings <- character()
    value <- withCallingHandlers(suppressMessages(code), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    list(value = value, warnings = warnings)
}

refactor_run_case <- function(inputs, case, parallel = FALSE, backed = FALSE, directory = NULL) {
    withr::with_seed(271828, refactor_capture({
        object <- refactor_prepare(inputs, case, backed, directory)
        prepared <- refactor_observe(object)
        object <- sceptre::run_calibration_check(object, n_calibration_pairs = 4L,
            calibration_group_size = 2L, print_progress = FALSE,
            parallel = parallel, n_processors = 2L, output_amount = case$output)
        calibrated <- refactor_observe(object)
        object <- sceptre::run_power_check(object, print_progress = FALSE,
            parallel = parallel, n_processors = 2L, output_amount = case$output)
        object <- sceptre::run_discovery_analysis(object, print_progress = FALSE,
            parallel = parallel, n_processors = 2L, output_amount = case$output)
        list(prepared = prepared, calibrated = calibrated, final = refactor_observe(object),
             rng_after = get(".Random.seed", envir = .GlobalEnv))
    }))
}

refactor_cache_observations <- function(inputs, case) {
    refactor_capture({
        object <- refactor_prepare(inputs, case)
        cold <- sceptre::run_discovery_analysis(object, print_progress = FALSE, output_amount = 2L)
        warm <- sceptre::run_discovery_analysis(cold, print_progress = FALSE, output_amount = 2L)
        cleared <- cold
        cleared@response_precomputations <- list()
        cleared <- sceptre::run_discovery_analysis(cleared, print_progress = FALSE, output_amount = 2L)
        changed <- sceptre::set_analysis_parameters(cold,
            discovery_pairs = inputs$discovery_pairs[rev(seq_len(nrow(inputs$discovery_pairs))), ],
            positive_control_pairs = inputs$positive_control_pairs,
            formula_object = stats::as.formula("~ batch + z", env = baseenv()),
            control_group = case$control, resampling_mechanism = case$mechanism,
            resampling_approximation = case$approximation, side = "right")
        side_cache <- changed@response_precomputations
        changed <- sceptre::run_discovery_analysis(changed, print_progress = FALSE, output_amount = 2L)
        formula_changed <- sceptre::set_analysis_parameters(cold,
            discovery_pairs = inputs$discovery_pairs,
            formula_object = stats::as.formula("~ batch", env = baseenv()),
            control_group = case$control, resampling_mechanism = case$mechanism)
        reassigned <- sceptre::assign_grnas(cold, method = "thresholding", threshold = 45, print_progress = FALSE)
        reassigned_cache <- reassigned@response_precomputations
        reassigned <- sceptre::run_qc(reassigned,
            response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
            n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L,
            additional_cells_to_remove = inputs$remove_cells)
        reassigned <- sceptre::run_discovery_analysis(reassigned, print_progress = FALSE, output_amount = 2L)
        reassigned_cold <- refactor_prepare(inputs, case, assignment_threshold = 45)
        reassigned_cold <- sceptre::run_discovery_analysis(reassigned_cold, print_progress = FALSE, output_amount = 2L)
        qc_changed <- sceptre::run_qc(cold,
            response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
            n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L,
            additional_cells_to_remove = c(inputs$remove_cells, 25L))
        list(cold = refactor_observe(cold), warm = refactor_observe(warm),
             cleared = refactor_observe(cleared), changed_side = refactor_observe(changed),
             retained_side_cache = side_cache,
             formula_cache = formula_changed@response_precomputations,
             reassigned_cache = reassigned_cache,
             reassigned = refactor_observe(reassigned), reassigned_cold = refactor_observe(reassigned_cold),
             changed_qc_cache = qc_changed@response_precomputations)
    })
}

refactor_index_observations <- function() {
    object <- methods::new("sceptre_object")
    object <- refactor_private("set_response_matrix")(object, matrix(101:112, nrow = 1L))
    object@low_moi <- TRUE
    object@control_group_complement <- FALSE
    object@cells_in_use <- c(2L, 4L, 5L, 8L, 11L)
    raw <- list(grna_group_idxs = list(A = c(5L, 4L)),
                indiv_nt_grna_idxs = list(ntA = 8L, ntB = c(2L, 11L)))
    # Use the actual baseline assignment schema (all_nt_idxs is derived below).
    object@grna_assignments_raw <- raw
    before <- serialize(raw, NULL)
    out <- refactor_private("update_grna_assignments_given_qc")(object)
    again <- refactor_private("update_grna_assignments_given_qc")(object)
    groups <- out@grna_assignments
    combined <- c(groups$grna_group_idxs$A, groups$all_nt_idxs)
    list(assignments = groups, repeated = again@grna_assignments,
         caller_unchanged = identical(serialize(raw, NULL), before),
         treatment_source = object@cells_in_use[groups$grna_group_idxs$A],
         nt_source = object@cells_in_use[groups$all_nt_idxs],
         combined_source = object@cells_in_use[combined],
         extracted_response = (101:112)[object@cells_in_use][combined],
         extracted_covariates = cbind(intercept = 1, cell = 101:112)[object@cells_in_use, , drop = FALSE][combined, , drop = FALSE])
}

refactor_edge_observations <- function(inputs, case) {
    outcome <- function(code) {
        tryCatch(refactor_capture(code), error = function(e) list(error_class = class(e), error = conditionMessage(e)))
    }
    duplicate <- inputs
    duplicate$discovery_pairs <- rbind(inputs$discovery_pairs, inputs$discovery_pairs[1L, ])
    no_pairs <- inputs
    no_pairs$discovery_pairs <- data.frame(grna_target = "A", response_id = "r4")
    no_pairs$positive_control_pairs <- inputs$positive_control_pairs[FALSE, ]
    list(
        duplicate_requests = outcome(refactor_prepare(duplicate, case)),
        no_passing_pairs = outcome({
            object <- refactor_prepare(no_pairs, case)
            refactor_observe(sceptre::run_discovery_analysis(object, print_progress = FALSE))
        }),
        intercept_only_nt = outcome({
            object <- refactor_prepare(inputs, readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_nt_perm)
            object <- sceptre::set_analysis_parameters(object, discovery_pairs = inputs$discovery_pairs,
                formula_object = stats::as.formula("~ 1", env = baseenv()), control_group = "nt_cells",
                resampling_mechanism = "permutations")
            refactor_observe(sceptre::run_discovery_analysis(object, print_progress = FALSE))
        })
    )
}
