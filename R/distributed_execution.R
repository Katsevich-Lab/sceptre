# Protocol certificates protect accidental mismatch/corruption, not adversarial
# forgery. The implementation fingerprint binds loaded R/schema/method code and
# native bytes; deployments must also pin one exact package artifact per
# process.
# Fingerprint loaded implementations, rather than a declared version or source
# revision. Deparsing excludes source references, bytecode and closure
# environments
# (which contain process paths); the function expressions themselves are
# retained.
execution_function_code <- function(fun) {
    if (is.primitive(fun)) {
        return(list(primitive = deparse(fun)))
    }
    list(
        formals = deparse(formals(fun), width.cutoff = 500L),
        body = deparse(body(fun), width.cutoff = 500L)
    )
}

execution_native_certificate <- function() {
    dll <- getLoadedDLLs()[["sceptre"]]
    if (is.null(dll) || !file.exists(dll[["path"]])) {
        stop(
            "Cannot certify execution: the loaded sceptre native library is ",
            "unavailable."
        )
    }
    digest::digest(file = dll[["path"]], algo = "sha256")
}

execution_class_certificate <- function(name, namespace) {
    definition <- methods::getClassDef(name, where = namespace)
    if (!identical(definition@package, "sceptre")) {
        return(NULL)
    }
    list(
        name = name,
        slots = methods::getSlots(definition),
        prototype = definition@prototype,
        contains = sort(
            as.character(names(definition@contains)),
            method = "radix"
        ),
        virtual = definition@virtual,
        validity = if (is.null(definition@validity)) {
            NULL
        } else {
            execution_function_code(definition@validity)
        }
    )
}

execution_generic_certificate <- function(generic, namespace) {
    implementations <- methods::findMethods(generic, where = namespace)
    implementations <- Filter(
        function(method) identical(environment(method), namespace),
        as.list(implementations)
    )
    ordered <- sort(
        as.character(names(implementations)),
        method = "radix"
    )
    stats::setNames(
        lapply(ordered, function(signature) {
            method <- implementations[[signature]]
            # Include the effective implementation for each package signature
            # too:
            # a redefinition outside this namespace can change actual dispatch.
            selected <- methods::selectMethod(
                generic,
                as.character(method@target)
            )
            list(
                target = as.character(method@target),
                defined = as.character(method@defined),
                code = execution_function_code(method),
                selected = execution_function_code(selected)
            )
        }),
        ordered
    )
}

execution_build_certificate <- function() {
    namespace <- asNamespace("sceptre")
    bindings <- sort(ls(namespace, all.names = TRUE), method = "radix")
    functions <- lapply(bindings, function(name) {
        get(name, envir = namespace, inherits = FALSE)
    })
    names(functions) <- bindings
    functions <- lapply(Filter(is.function, functions), execution_function_code)
    classes <- lapply(
        sort(methods::getClasses(where = namespace), method = "radix"),
        execution_class_certificate,
        namespace = namespace
    )
    classes <- Filter(Negate(is.null), classes)
    generics <- sort(
        unique(as.character(methods::getGenerics(where = namespace))),
        method = "radix"
    )
    registered <- stats::setNames(
        lapply(generics, execution_generic_certificate, namespace = namespace),
        generics
    )
    list(
        R = execution_digest(list(
            functions = functions,
            classes = classes,
            methods = registered
        )),
        native = execution_native_certificate()
    )
}

execution_runtime_certificate <- function(object) {
    list(
        R = as.character(getRversion()),
        platform = R.version$platform,
        package = as.character(getNamespaceVersion(asNamespace("sceptre"))),
        build = execution_build_certificate(),
        stats = as.character(utils::packageVersion("stats")),
        Matrix = as.character(utils::packageVersion("Matrix")),
        Rcpp = as.character(utils::packageVersion("Rcpp")),
        ondisc = if (
            any(vapply(
                list(object@data@response_matrix, object@data@grna_matrix),
                function(store) methods::is(store, "BackedCountStore"),
                logical(1)
            ))
        ) {
            if ("ondisc" %in% loadedNamespaces()) {
                as.character(getNamespaceVersion(asNamespace("ondisc")))
            } else {
                "unloaded"
            }
        } else {
            NULL
        },
        BLAS = unname(extSoftVersion()["BLAS"]),
        matprod = getOption("matprod"),
        RNGkind = RNGkind()
    )
}

execution_method_certificate <- function() {
    list(
        response = response_fit_method_id(),
        perturbation = perturbation_fit_method_id(),
        resampling = "global-native-draw-sharing-v1"
    )
}

execution_digest <- function(value) digest::digest(value, algo = "sha256")
execution_row_column <- function() ".sceptre_global_pair_row"

require_execution_backend <- function(object) {
    if (
        !isTRUE(object@spec@run_permutations) &&
            any(vapply(
                list(object@data@response_matrix, object@data@grna_matrix),
                function(store) methods::is(store, "BackedCountStore"),
                logical(1)
            ))
    ) {
        stop(
            "Backed CRT remains unsupported; use permutation resampling or ",
            "in-memory counts."
        )
    }
    invisible(NULL)
}

execution_context_certificate <- function(object) {
    list(
        data = object@data@identity,
        stages = lapply(object@preparation@stages, function(stage) {
            stage$identity
        })
    )
}

execution_units <- function(context, plan, n_work_units, permutations) {
    field <- if (permutations) "response_idx" else "group_idx"
    keys <- unique(context$index$pair_rows[[field]])
    if (!length(keys)) {
        stop("Distributed execution requires eligible hypothesis rows.")
    }
    n <- min(as.integer(n_work_units), length(keys))
    chunks <- split(
        keys,
        as.integer(floor((seq_along(keys) - 1L) * n / length(keys))) + 1L
    )
    eligible <- which(plan@pairs$pass_qc)
    lapply(seq_along(chunks), function(i) {
        selected <- as.integer(chunks[[i]])
        rows <- unlist(
            lapply(selected, function(key) {
                eligible[which(context$index$pair_rows[[field]] == key)]
            }),
            use.names = FALSE
        )
        list(
            id = sprintf("unit-%06d", i),
            keys = selected,
            row_ids = as.integer(rows)
        )
    })
}

execution_full_context <- function(object, plan) {
    context <- prepare_analysis_context(
        object,
        plan@analysis,
        dplyr::filter(plan@pairs, pass_qc),
        plan = plan
    )
    column <- execution_row_column()
    if (column %in% names(plan@pairs)) {
        stop("The hypothesis family uses a reserved execution row column.")
    }
    numerical <- c(
        "p_value",
        "fold_change",
        "se_fold_change",
        "log_2_fold_change",
        "stage",
        "z_orig",
        "xi",
        "omega",
        "alpha"
    )
    if (
        any(names(plan@pairs) %in% numerical) ||
            any(grepl("^z_null_", names(plan@pairs)))
    ) {
        stop(
            "Prepared hypothesis metadata collides with generated numerical ",
            "output columns."
        )
    }
    # Append after the canonical plan/context comparison; never subset the
    # index,
    # assignments, design, fitting view or global permutation dimensions.
    context$pairs[[column]] <- which(plan@pairs$pass_qc)
    context
}

validate_execution_envelope <- function(object, execution) {
    fields <- c(
        "protocol",
        "runtime",
        "methods",
        "context",
        "plan",
        "output_amount",
        "units",
        "identity"
    )
    if (
        !is.list(execution) ||
            anyDuplicated(names(execution)) ||
            !setequal(names(execution), fields) ||
            !identical(execution$protocol, 1L) ||
            !identical(
                execution$runtime,
                execution_runtime_certificate(object)
            ) ||
            !identical(execution$methods, execution_method_certificate())
    ) {
        stop(
            "Execution manifest protocol, runtime or method certificate is ",
            "incompatible."
        )
    }
    invisible(NULL)
}

validate_execution_manifest <- function(object, execution) {
    require_stage(object, "run_qc")
    require_execution_backend(object)
    validate_execution_envelope(object, execution)
    if (
        !identical(
            execution$identity,
            execution_digest(execution[setdiff(names(execution), "identity")])
        )
    ) {
        stop("Execution manifest certificate does not match its payload.")
    }
    if (!identical(execution$context, execution_context_certificate(object))) {
        stop(
            "Execution manifest no longer matches the full prepared screen ",
            "context."
        )
    }
    plan <- validate_run_plan(object, execution$plan)
    if (
        !plan@analysis %in% "run_calibration_check" &&
            !identical(plan@pairs, analysis_pairs(object, plan@analysis))
    ) {
        stop("Execution plan does not cover the complete prepared family.")
    }
    if (
        !identical(
            execution$output_amount,
            as.integer(execution$output_amount)
        ) ||
            length(execution$output_amount) != 1L ||
            !execution$output_amount %in% seq_len(3L)
    ) {
        stop("Invalid execution output level.")
    }
    context <- execution_full_context(object, plan)
    units <- execution_units(
        context,
        plan,
        length(execution$units),
        object@spec@run_permutations
    )
    if (!identical(execution$units, units)) {
        stop("Execution work units do not match canonical whole-key traversal.")
    }
    list(plan = plan, context = context)
}

prepare_execution_plan <- function(
    object,
    analysis,
    n_calibration_pairs,
    calibration_group_size
) {
    if (identical(analysis, "run_calibration_check")) {
        run_calibration_check_pt_1(
            object,
            n_calibration_pairs,
            calibration_group_size,
            parallel = FALSE,
            n_processors = 2L
        )$plan
    } else {
        pairs <- analysis_pairs(object, analysis)
        if (!nrow(pairs) || !any(pairs$pass_qc)) {
            stop("The analysis has no eligible hypotheses.")
        }
        make_run_plan(object, analysis, pairs)
    }
}

#' Prepare serializable work units for one complete analysis
#'
#' Prepare once, transport the manifest with a separately saved/reopened screen,
#' execute every unit, and collect the proposals once. Hypothesis adjustment is
#' performed across the complete family at collection. The runtime and package
#' artifact must match in all processes; backing files remain immutable.
#' @param sceptre_object A screen with successful assignment and both QC stages.
#' @param analysis Canonical `run_calibration_check`, `run_power_check`, or
#'   `run_discovery_analysis`.
#' @param n_work_units Positive requested number of contiguous whole-key units.
#'   At most one unit per eligible response (permutation) or guide group (CRT).
#' @param output_amount Output level 1, 2 or 3, as in the public analyses.
#' @param n_calibration_pairs,calibration_group_size Optional calibration
#' choices.
#'   Calibration groups are sampled only during preparation.
#' @return An owned serializable manifest containing a frozen RunPlan, ordered
#'   work units and protocol/runtime certificates. Prepared metadata must not
#'   reuse generated numerical output column names. It contains no counts or
#'   live
#'   handles. This is a value artifact, not a saved completed analysis.
#' @examples
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#' pairs <- construct_cis_pairs(object)
#' object <- object |>
#'   set_analysis_parameters(
#'     side = "left", resampling_mechanism = "permutations",
#'     discovery_pairs = pairs
#'   ) |>
#'   assign_grnas(method = "thresholding") |>
#'   run_qc()
#' execution <- prepare_execution(object, "run_discovery_analysis", 2L)
#' # Each worker uses the same full prepared screen and one stable unit ID.
#' proposals <- lapply(execution$units, function(unit) {
#'   execute_work_unit(object, execution, unit$id)
#' })
#' completed <- collect_execution(object, execution, rev(proposals))
#' head(get_result(completed, "run_discovery_analysis"))
#' @export
prepare_execution <- function(
    sceptre_object,
    analysis,
    n_work_units = 1L,
    output_amount = 1L,
    n_calibration_pairs = NULL,
    calibration_group_size = NULL
) {
    object <- migrate_sceptre_object(sceptre_object)
    if (has_current_wilcoxon_result(object)) {
        stop(
            "Experimental Wilcoxon is unsupported in distributed execution. ",
            "Use run_wilcoxon() with direct execution or a separate ",
            "prepared SCEPTRE object."
        )
    }
    require_stage(object, "run_qc")
    require_execution_backend(object)
    analysis <- analysis_record_name(analysis)
    check_native_scalar_count(n_work_units, "n_work_units")
    if (
        n_work_units < 1L ||
            length(output_amount) != 1L ||
            is.na(output_amount) ||
            !output_amount %in% seq_len(3L)
    ) {
        stop("Use positive work-unit count and output_amount 1, 2 or 3.")
    }
    plan <- prepare_execution_plan(
        object,
        analysis,
        n_calibration_pairs,
        calibration_group_size
    )
    context <- execution_full_context(object, plan)
    manifest <- list(
        protocol = 1L,
        runtime = execution_runtime_certificate(object),
        methods = execution_method_certificate(),
        context = execution_context_certificate(object),
        plan = plan,
        output_amount = as.integer(output_amount),
        units = execution_units(
            context,
            plan,
            n_work_units,
            object@spec@run_permutations
        )
    )
    manifest$identity <- execution_digest(manifest)
    owned_result_value(manifest)
}

execute_unit_analysis <- function(object, prepared, controls) {
    context <- prepared$context
    if (object@spec@run_permutations) {
        draws <- get_synthetic_permutation_idxs(
            object@preparation@analysis_groups,
            B = context$test_options$total_resamples,
            calibration_check = context$index$calibration_check,
            control_group_complement = context$index$control_group_complement,
            calibration_group_size = if (
                length(prepared$plan@calibration_group_size)
            ) {
                prepared$plan@calibration_group_size
            } else {
                1L
            },
            n_cells = length(context$index$analysis_to_source)
        )
        run_permutation_analysis(context, controls, draws)
    } else {
        run_crt_analysis(context, controls)
    }
}

execution_unit_fits <- function(object, prepared, unit, output) {
    context <- prepared$context
    response_ids <- unique(as.character(prepared$plan@pairs$response_id[
        unit$row_ids
    ]))
    response_fits <- if (context$shared_response_fit) {
        output$response_precomputations[
            intersect(names(output$response_precomputations), response_ids)
        ]
    } else {
        list()
    }
    groups <- if (object@spec@run_permutations) {
        character()
    } else {
        as.character(context$index$group_ids[unit$keys])
    }
    perturbation_fits <- if (is.null(output$perturbation_fits)) {
        NULL
    } else {
        stats::setNames(
            lapply(groups, function(id) output$perturbation_fits[[id]]),
            groups
        )
    }
    list(response = response_fits, perturbation = perturbation_fits)
}

execution_unit_by_id <- function(execution, unit_id) {
    ids <- vapply(execution$units, function(unit) unit$id, character(1))
    if (
        !is.character(unit_id) ||
            length(unit_id) != 1L ||
            is.na(unit_id) ||
            !unit_id %in% ids
    ) {
        stop("Unknown execution unit ID.")
    }
    execution$units[[match(unit_id, ids)]]
}

#' Execute one certified work unit without committing an analysis
#' @param sceptre_object The full matching prepared screen, with live count
#' stores.
#' @param execution Manifest returned by `prepare_execution()`.
#' @param unit_id One exact stable unit ID from `execution$units`.
#' @return An owned serializable proposal containing global row IDs, raw tests,
#'   compact fits and digest certificates of their full context scopes. Ordered
#'   cell vectors are retained in the manifest, not repeated in proposals.
#'   No guide integration or family correction is
#'   performed. Workers must return all proposals to `collect_execution()`.
#' @inherit prepare_execution examples
#' @export
execute_work_unit <- function(sceptre_object, execution, unit_id) {
    object <- migrate_sceptre_object(sceptre_object)
    prepared <- validate_execution_manifest(object, execution)
    unit <- execution_unit_by_id(execution, unit_id)
    context <- prepared$context
    controls <- list(
        output_amount = execution$output_amount,
        print_progress = FALSE,
        parallel = FALSE,
        n_processors = 1L,
        log_dir = tempdir(),
        work_keys = unit$keys
    )
    output <- execute_unit_analysis(object, prepared, controls)
    fits <- execution_unit_fits(object, prepared, unit, output)
    proposal <- list(
        protocol = execution$protocol,
        runtime = execution$runtime,
        methods = execution$methods,
        execution = execution$identity,
        unit_id = unit$id,
        result = output$result,
        response_fits = fits$response,
        response_identity = if (context$shared_response_fit) {
            execution_digest(context$fit_identity)
        } else {
            NULL
        },
        perturbation_fits = fits$perturbation,
        perturbation_identity = if (is.null(fits$perturbation)) {
            NULL
        } else {
            execution_digest(context$perturbation_identity)
        }
    )
    proposal$identity <- execution_digest(proposal)
    owned_result_value(proposal)
}

validate_execution_proposal_envelope <- function(proposal, unit, execution) {
    fields <- c(
        "protocol",
        "runtime",
        "methods",
        "execution",
        "unit_id",
        "result",
        "response_fits",
        "response_identity",
        "perturbation_fits",
        "perturbation_identity",
        "identity"
    )
    if (
        !is.list(proposal) ||
            anyDuplicated(names(proposal)) ||
            !setequal(names(proposal), fields) ||
            !identical(proposal$protocol, execution$protocol) ||
            !identical(proposal$runtime, execution$runtime) ||
            !identical(proposal$methods, execution$methods) ||
            !identical(proposal$execution, execution$identity) ||
            !identical(proposal$unit_id, unit$id) ||
            !identical(
                proposal$identity,
                execution_digest(proposal[setdiff(names(proposal), "identity")])
            )
    ) {
        stop(
            "Work proposal certificate does not match the execution or payload."
        )
    }
    invisible(NULL)
}

validate_execution_proposal <- function(proposal, unit, execution, context) {
    validate_execution_proposal_envelope(proposal, unit, execution)
    result <- proposal$result
    column <- execution_row_column()
    if (
        !is.data.frame(result) ||
            !all(c(column, names(execution$plan@pairs)) %in% names(result)) ||
            !identical(result[[column]], unit$row_ids) ||
            anyDuplicated(result[[column]])
    ) {
        stop(
            "Work proposal has missing, duplicated or reordered global ",
            "hypothesis rows."
        )
    }
    expected <- execution$plan@pairs[unit$row_ids, , drop = FALSE]
    if (
        !identical(
            lapply(as.list(result)[names(expected)], identity),
            lapply(expected, identity)
        )
    ) {
        stop("Work proposal changed hypothesis labels or QC metadata.")
    }
    response_ids <- unique(as.character(expected$response_id))
    if (
        !is.list(proposal$response_fits) ||
            anyDuplicated(names(proposal$response_fits)) ||
            (context$shared_response_fit &&
                (!setequal(names(proposal$response_fits), response_ids) ||
                    !identical(
                        proposal$response_identity,
                        execution_digest(context$fit_identity)
                    ))) ||
            (!context$shared_response_fit &&
                (length(proposal$response_fits) ||
                    !is.null(proposal$response_identity)))
    ) {
        stop("Work proposal response fit scope is incompatible.")
    }
    invisible(NULL)
}

execution_perturbation_identity <- function(context, group) {
    if (context$shared_response_fit) {
        trt <- make_treatment_index_lookup(context$index)(group)$trt_idxs
        perturbation_fit_identity(
            context$common_view@source_cells,
            context$design,
            trt
        )
    } else {
        view <- make_fit_view(context$index, group)
        perturbation_fit_identity(
            view@source_cells,
            context$design[view@fit_to_analysis, , drop = FALSE],
            view@trt_fit_pos
        )
    }
}

order_execution_proposals <- function(execution, proposals) {
    if (!is.list(proposals) || length(proposals) != length(execution$units)) {
        stop("Collection requires every expected work unit exactly once.")
    }
    # Own caller-retained fits and raw tables before any cache admission.
    proposals <- owned_result_value(proposals)
    expected_ids <- vapply(
        execution$units,
        function(unit) unit$id,
        character(1)
    )
    ids <- vapply(
        proposals,
        function(proposal) {
            if (
                !is.list(proposal) ||
                    !is.character(proposal$unit_id) ||
                    length(proposal$unit_id) != 1L
            ) {
                stop("Malformed work proposal.")
            }
            proposal$unit_id
        },
        character(1)
    )
    if (anyNA(ids) || anyDuplicated(ids) || !setequal(ids, expected_ids)) {
        stop("Collection has duplicate, missing or unexpected unit IDs.")
    }
    proposals <- proposals[match(expected_ids, ids)]
    proposals
}

validate_execution_perturbation_fits <- function(
    context,
    groups,
    perturbations,
    proposed
) {
    for (group in groups) {
        fit <- proposed[[group]]
        actual_identity <- execution_perturbation_identity(
            context,
            group
        )
        existing <- perturbations[[group]]
        if (
            !is.null(existing) &&
                identical(existing$identity, actual_identity) &&
                !identical(existing, fit)
        ) {
            stop(
                "Conflicting perturbation fit for an existing valid warm ",
                "cache entry."
            )
        }
        if (
            !is.null(fit) &&
                !identical(
                    fit$identity,
                    actual_identity
                )
        ) {
            stop(
                "Work proposal perturbation fit does not match its actual ",
                "ordered inputs."
            )
        }
    }
    invisible(NULL)
}

collect_execution_perturbations <- function(
    object,
    context,
    unit,
    proposal,
    perturbations
) {
    groups <- if (
        object@spec@run_permutations || context$index$calibration_check
    ) {
        character()
    } else {
        as.character(context$index$group_ids[unit$keys])
    }
    if (length(groups)) {
        if (
            !is.list(proposal$perturbation_fits) ||
                !setequal(names(proposal$perturbation_fits), groups) ||
                anyDuplicated(names(proposal$perturbation_fits)) ||
                !identical(
                    proposal$perturbation_identity,
                    execution_digest(context$perturbation_identity)
                )
        ) {
            stop("Work proposal perturbation scope is incompatible.")
        }
        validate_execution_perturbation_fits(
            context,
            groups,
            perturbations,
            proposal$perturbation_fits
        )
        perturbations <- merge_perturbation_fits(
            perturbations,
            proposal$perturbation_fits
        )
    } else if (!is.null(proposal$perturbation_fits)) {
        stop("Unexpected perturbation cache in this execution mode.")
    }
    perturbations
}

collect_execution_fits <- function(object, execution, proposals, context) {
    responses <- if (context$shared_response_fit) {
        context$response_fits
    } else {
        object@fits@response_fits
    }
    perturbations <- context$perturbation_fits
    for (i in seq_along(proposals)) {
        proposal <- proposals[[i]]
        unit <- execution$units[[i]]
        validate_execution_proposal(proposal, unit, execution, context)
        for (id in names(proposal$response_fits)) {
            if (
                !is.null(responses[[id]]) &&
                    !identical(responses[[id]], proposal$response_fits[[id]])
            ) {
                stop("Conflicting response fits across work units.")
            }
            responses[id] <- proposal$response_fits[id]
        }
        perturbations <- collect_execution_perturbations(
            object,
            context,
            unit,
            proposal,
            perturbations
        )
    }
    list(responses = responses, perturbations = perturbations)
}

order_execution_response_fits <- function(object, context, responses) {
    # Cache insertion order follows the direct global response traversal, even
    # when one response was fitted independently by several CRT guide units.
    if (context$shared_response_fit) {
        existing <- names(context$response_fits)
        global <- as.character(context$index$response_ids[unique(
            context$index$pair_rows$response_idx
        )])
        responses <- responses[c(existing, setdiff(global, existing))]
    }
    attr(responses, "fit_identity") <- if (context$shared_response_fit) {
        context$fit_identity
    } else {
        object@fits@response_identity
    }
    responses
}

assemble_execution_result <- function(object, plan, proposals) {
    raw <- bind_analysis_results(lapply(proposals, function(proposal) {
        proposal$result
    }))
    raw[[execution_row_column()]] <- NULL
    assembled <- assemble_pair_results(raw, plan@pairs)
    switch(
        plan@analysis,
        run_calibration_check = process_calibration_result(assembled, object),
        run_power_check = apply_grouping_to_result(assembled, object),
        run_discovery_analysis = process_discovery_result(assembled, object)
    )
}

#' Collect every work proposal and commit one complete-family result
#' @param sceptre_object The full matching prepared screen.
#' @param execution Manifest returned by `prepare_execution()`.
#' @param proposals List containing exactly one proposal per expected unit, in
#'   any arrival order. Missing, duplicated or incompatible proposals are
#'   errors.
#' @return Updated screen with one completed result and compatible compact fits.
#'   Errors leave the caller unchanged. Family correction is performed once,
#'   after restoring canonical traversal and adding failed-QC rows.
#' @inherit prepare_execution examples
#' @export
collect_execution <- function(sceptre_object, execution, proposals) {
    object <- migrate_sceptre_object(sceptre_object)
    prepared <- validate_execution_manifest(object, execution)
    proposals <- order_execution_proposals(execution, proposals)
    context <- prepared$context
    fits <- collect_execution_fits(object, execution, proposals, context)
    responses <- fits$responses
    perturbations <- fits$perturbations
    responses <- order_execution_response_fits(object, context, responses)
    result <- assemble_execution_result(object, prepared$plan, proposals)
    output <- list(
        response_precomputations = responses,
        perturbation_fits = if (
            object@spec@run_permutations || context$index$calibration_check
        ) {
            NULL
        } else {
            perturbations
        }
    )
    if (!is.null(output$perturbation_fits)) {
        attr(
            output$perturbation_fits,
            "fit_identity"
        ) <- context$perturbation_identity
    }
    object <- commit_result_record(
        object,
        prepared$plan@analysis,
        result,
        prepared$plan@pairs,
        execution$output_amount,
        plan = prepared$plan
    )
    commit_analysis_fits(object, output)
}
