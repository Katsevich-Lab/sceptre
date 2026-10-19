#' Internal per-run scientific plan
#'
#' Typed family, budgets and frozen input identity, used during execution and
#' retained as result provenance. Calibration-only integer slots are empty for
#' other analyses. Structural validity is separate from screen-context matching.
#' @name RunPlan-class
#' @keywords internal
NULL

# One canonical value representation serves hashing and the schema-three list
# boundary. In particular, keep the old field order and explicitly named NULLs.
run_plan_payload <- function(plan, include_identity = FALSE) {
    fields <- c(
        "analysis",
        "pairs",
        "B1",
        "B2",
        "B3",
        "n_calibration_pairs",
        "calibration_group_size"
    )
    if (methods::is(plan, "RunPlan")) {
        payload <- stats::setNames(
            lapply(fields, function(name) slot(plan, name)),
            fields
        )
        for (name in c("n_calibration_pairs", "calibration_group_size")) {
            if (!length(payload[[name]])) payload[name] <- list(NULL)
        }
    } else {
        if (
            !is.list(plan) ||
                is.null(names(plan)) ||
                anyDuplicated(names(plan)) ||
                !setequal(
                    names(plan),
                    c(fields, if ("identity" %in% names(plan)) "identity")
                )
        ) {
            stop("Unrecognized run-plan payload fields.")
        }
        payload <- plan[fields]
    }
    if (include_identity) {
        payload["identity"] <- list(
            if (methods::is(plan, "RunPlan")) plan@identity else plan$identity
        )
    }
    payload
}

# Ordinary data frames already have R copy-on-modify behavior. Snapshot only
# data.tables, whose by-reference writes would otherwise change a chosen family.
# Preserve attribute order because schema-three hashes include its representation.
copy_run_plan_pairs <- function(pairs) {
    if (!data.table::is.data.table(pairs)) {
        return(pairs)
    }
    copy <- data.table::copy(pairs)
    attributes(copy) <- attributes(copy)[names(attributes(pairs))]
    copy
}

# Lossless boundary conversion, not a constructor that re-signs old provenance.
# Old list plans must carry their original identity and exactly the known fields.
as_run_plan <- function(plan) {
    if (methods::is(plan, "RunPlan")) {
        validObject(plan)
    } else {
        payload <- run_plan_payload(plan, include_identity = TRUE)
        if (!"identity" %in% names(plan) || !is.list(payload$identity)) {
            stop("A saved run plan requires its original input identity.")
        }
        for (name in c("n_calibration_pairs", "calibration_group_size")) {
            if (is.null(payload[[name]])) payload[name] <- list(integer())
        }
        pairs <- payload$pairs
        plan <- do.call(methods::new, c(list(Class = "RunPlan"), payload))
        # Initialization can coerce a data.table subclass to data.frame. Keep
        # the original validated payload so conversion never changes its hash.
        methods::slot(plan, "pairs", check = FALSE) <- pairs
        validObject(plan)
    }
    # Verify the original payload certificate before any contextual substitution
    # (such as checking a saved result against the current prepared family).
    if (
        !identical(
            plan@identity$execution,
            digest::digest(run_plan_payload(plan), algo = "sha256")
        )
    ) {
        stop("The run plan no longer matches its frozen execution payload.")
    }
    plan
}

setValidity("RunPlan", function(object) {
    analyses <- analysis_names()
    if (
        length(object@analysis) != 1L ||
            is.na(object@analysis) ||
            !object@analysis %in% analyses
    ) {
        return("RunPlan requires one canonical analysis mode")
    }
    pairs <- object@pairs
    if (
        !all(c("response_id", "grna_group", "pass_qc") %in% names(pairs)) ||
            anyDuplicated(names(pairs)) ||
            !is.logical(pairs$pass_qc) ||
            anyNA(pairs$pass_qc)
    ) {
        return(
            "RunPlan requires prepared hypothesis columns and complete logical QC status"
        )
    }
    for (name in c("response_id", "grna_group")) {
        ids <- pairs[[name]]
        if (!is.character(ids) || anyNA(ids) || any(!nzchar(ids))) {
            return(
                "RunPlan hypothesis identifiers must be nonmissing character values"
            )
        }
    }
    budget_error <- tryCatch(
        {
            checked_resampling_total(object@B1, object@B2, object@B3)
            NULL
        },
        error = function(error) conditionMessage(error)
    )
    if (!is.null(budget_error)) {
        return(budget_error)
    }
    calibration <- identical(object@analysis, "run_calibration_check")
    sizes <- list(object@n_calibration_pairs, object@calibration_group_size)
    if (calibration) {
        if (
            any(vapply(
                sizes,
                function(x) length(x) != 1L || anyNA(x) || x < 1L,
                logical(1)
            ))
        ) {
            return("Calibration RunPlan sizes must be positive scalar integers")
        }
        if (
            !nrow(pairs) ||
                nrow(pairs) > object@n_calibration_pairs ||
                any(
                    lengths(strsplit(pairs$grna_group, "&", fixed = TRUE)) !=
                        object@calibration_group_size
                )
        ) {
            return(
                "Selected calibration pairs do not match the requested run scope"
            )
        }
    } else if (any(lengths(sizes))) {
        return(
            "Calibration-only RunPlan sizes must be empty for other analyses"
        )
    }
    identity <- object@identity
    fields <- c(
        "data",
        "cells",
        "design",
        "assignments",
        "choices",
        "execution"
    )
    if (
        !setequal(names(identity), fields) ||
            anyDuplicated(names(identity)) ||
            !is.integer(identity$cells) ||
            anyNA(identity$cells) ||
            any(identity$cells < 1L) ||
            anyDuplicated(identity$cells) ||
            !is.list(identity$choices) ||
            any(vapply(
                identity[c("data", "design", "assignments", "execution")],
                function(x) {
                    !is.character(x) ||
                        length(x) != 1L ||
                        anyNA(x) ||
                        !nzchar(x)
                },
                logical(1)
            ))
    ) {
        return("RunPlan requires its frozen input identity")
    }
    TRUE
})
