# The saved object records which workflow steps have finished. Define their
# order and prerequisites here. Calibration, power and discovery are separate
# final analyses: finishing one does not mark either of the others as finished.
#' List the workflow steps and the order in which they may run
#' @return Named integer vector for the seven saved workflow flags. Import,
#'   parameters, assignment and QC have increasing step numbers. Calibration,
#'   power and discovery all have number 5: none requires another final analysis.
#' @noRd
workflow_steps <- function() {
    c(
        import_data = 1L,
        set_analysis_parameters = 2L,
        assign_grnas = 3L,
        run_qc = 4L,
        run_calibration_check = 5L,
        run_power_check = 5L,
        run_discovery_analysis = 5L
    )
}

#' Mark import as complete and all later workflow steps as incomplete
#' @return Named logical vector in workflow order, ready to save in the object.
#' @noRd
workflow_initial_state <- function() {
    steps <- workflow_steps()
    stats::setNames(steps == 1L, names(steps))
}

#' Check that prerequisites finished and mark the requested step as complete
#' @param functs_called Named logical vector containing all seven workflow flags.
#' @param curr_funct Name of the workflow function being called.
#' @return Updated flags: the current step becomes TRUE and steps later in the
#'   workflow become FALSE. Other final analyses at the same step number, and
#'   any extra unrecognized flags, retain their values.
#' @details Require all steps at the immediately preceding number. For example,
#'   QC requires guide assignment. Repeating discovery does not clear calibration
#'   or power, since all three share the final step number. This updates flags
#'   only; the caller runs the analysis and handles saved fits separately.
#' @noRd
workflow_transition <- function(functs_called, curr_funct) {
    steps <- workflow_steps()
    curr_rank <- steps[[curr_funct]]
    prerequisites <- names(steps)[steps == curr_rank - 1L]
    if (!all(functs_called[prerequisites])) {
        stop(paste0(
            "The function",
            if (length(prerequisites) >= 2L) "s " else " ",
            paste0("`", prerequisites, "()`", collapse = ", "),
            " must be called before `",
            curr_funct,
            "()` is called."
        ))
    }
    functs_called[curr_funct] <- TRUE
    functs_called[names(steps)[steps > curr_rank]] <- FALSE
    functs_called
}

#' Put saved workflow flags in the order used for display
#' @param functs_called Named logical vector of saved workflow flags.
#' @return Known flags in workflow order, followed by unrecognized flags in their
#'   original order. Values are unchanged; displaying a flag does not check
#'   prerequisites or decide whether a step is complete.
#' @noRd
workflow_status <- function(functs_called) {
    # Keep extra flags from saved objects in the display, after known steps,
    # so showing status does not silently discard that information.
    order <- c(
        intersect(names(workflow_steps()), names(functs_called)),
        setdiff(names(functs_called), names(workflow_steps()))
    )
    functs_called[order]
}

#' Update the object's saved workflow status for the requested function
#' @param sceptre_object Sceptre object; the caller must keep the returned object.
#' @param curr_funct Name of the workflow function being called.
#' @return The object with `functs_called` updated by `workflow_transition()` and
#'   `last_function_called` set to this function. Other fields are unchanged.
#'   This helper records status; it does not run the numerical analysis.
#' @noRd
perform_status_check_and_update <- function(sceptre_object, curr_funct) {
    sceptre_object@functs_called <- workflow_transition(
        sceptre_object@functs_called,
        curr_funct
    )
    sceptre_object@last_function_called <- curr_funct
    sceptre_object
}
