# These tags are persisted in cache identities and result provenance. Changing
# a method requires a deliberate new tag and qualification of its consumers.
response_fit_method_id <- function() "poisson-glm.fit-then-dispersion-v1"
perturbation_fit_method_id <- function(method = "glm.fit") {
    if (identical(method, "glm.fit")) {
        "binomial-glm.fit-linkinv-v1"
    } else {
        "binomial-ziang-newton-linkinv-fallback-glm.fit-v1"
    }
}
analysis_names <- function() {
    c("run_calibration_check", "run_power_check", "run_discovery_analysis")
}
