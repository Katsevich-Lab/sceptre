# The S4 flags remain the saved/public representation. This table is the single
# owner of prerequisites, invalidation order, and display/initialization order.
# The three analyses are siblings: completing one does not complete another.
#' Define the existing workflow ranks and saved-flag order
#'
#' @return Named integer vector of seven ranks. Calibration, power and discovery
#'   share rank five: they are independent siblings, not sequential prerequisites.
#' @noRd
workflow_steps <- function() {
    c(import_data = 1L, set_analysis_parameters = 2L, assign_grnas = 3L,
      run_qc = 4L, run_calibration_check = 5L, run_power_check = 5L,
      run_discovery_analysis = 5L)
}

#' Initialize the saved workflow flags after a successful import
#'
#' @return Named logical vector in workflow order: import true, all other flags
#'   false. Produces no object mutation, RNG draws or file side effects.
#' @noRd
workflow_initial_state <- function() {
    steps <- workflow_steps()
    stats::setNames(steps == 1L, names(steps))
}

#' Check prerequisites and update the existing workflow flags
#'
#' @param functs_called Named logical vector with all seven supported saved flags.
#' @param curr_funct One supported workflow function name.
#' @return Updated flag vector: current flag true, all strictly later ranks false;
#'   flags at the same rank and unknown extra flags remain unchanged.
#' @details Requires all flags at the immediately preceding rank, preserving the
#'   existing prerequisite error text. Repeating a final analysis preserves its
#'   sibling flags. Callers supply valid names/flags; the input is not mutated and
#'   this operation neither runs analyses nor changes caches, RNG or files.
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
            " must be called before `", curr_funct, "()` is called."
        ))
    }
    functs_called[curr_funct] <- TRUE
    functs_called[names(steps)[steps > curr_rank]] <- FALSE
    functs_called
}

#' Order workflow flags for display without changing their meaning
#'
#' @param functs_called Named saved logical workflow flags.
#' @return Present known flags in workflow order, followed by unknown flags in
#'   their original order for saved-object compatibility. Does not infer completion
#'   or enforce prerequisites; the input is unchanged.
#' @noRd
workflow_status <- function(functs_called) {
    # Old supported objects use this same order. Preserve any unknown flags
    # rather than silently discarding information in a display operation.
    order <- c(intersect(names(workflow_steps()), names(functs_called)),
               setdiff(names(functs_called), names(workflow_steps())))
    functs_called[order]
}

#' Apply a workflow transition to the existing S4 facade
#'
#' @param sceptre_object Sceptre S4 facade; caller must retain the returned object.
#' @param curr_funct One supported workflow function name.
#' @return Object with transitioned `functs_called` and updated
#'   `last_function_called`; other slots are unchanged. Preserves the saved flag
#'   schema without running numerical work, invalidating fits, using RNG or files.
#' @noRd
perform_status_check_and_update <- function(sceptre_object, curr_funct) {
    sceptre_object@functs_called <- workflow_transition(
        sceptre_object@functs_called, curr_funct
    )
    sceptre_object@last_function_called <- curr_funct
    sceptre_object
}
