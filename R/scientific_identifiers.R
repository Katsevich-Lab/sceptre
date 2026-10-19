# These tags are persisted in cache identities and result provenance. Changing
# a method requires a deliberate new tag and qualification of its consumers.
response_fit_method_id <- function() "poisson-glm.fit-then-dispersion-v1"
perturbation_fit_method_id <- function() "binomial-glm.fit-linkinv-v1"
analysis_names <- function() c("run_calibration_check", "run_power_check",
    "run_discovery_analysis")
