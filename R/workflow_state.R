# The S4 flags remain the saved/public representation. This table is the single
# owner of prerequisites, invalidation order, and display/initialization order.
# The three analyses are siblings: completing one does not complete another.
workflow_steps <- function() {
    c(import_data = 1L, set_analysis_parameters = 2L, assign_grnas = 3L,
      run_qc = 4L, run_calibration_check = 5L, run_power_check = 5L,
      run_discovery_analysis = 5L)
}

workflow_initial_state <- function() {
    steps <- workflow_steps()
    stats::setNames(steps == 1L, names(steps))
}

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

workflow_status <- function(functs_called) {
    # Old supported objects use this same order. Preserve any unknown flags
    # rather than silently discarding information in a display operation.
    order <- c(intersect(names(workflow_steps()), names(functs_called)),
               setdiff(names(functs_called), names(workflow_steps())))
    functs_called[order]
}

perform_status_check_and_update <- function(sceptre_object, curr_funct) {
    sceptre_object@functs_called <- workflow_transition(
        sceptre_object@functs_called, curr_funct
    )
    sceptre_object@last_function_called <- curr_funct
    sceptre_object
}
