# An opt-in direct-R example on simulated counts with a known down effect.
# Run in an isolated library after installation; candidate is experimental.
library(sceptre)
set.seed(71019)
n <- 600L
z <- rnorm(n)
x <- rbinom(n, 1L, plogis(-1.4 + .7*z))
Y <- rbind(null=rpois(n,exp(1+.5*z)),
    down=rpois(n,exp(1+.5*z+log(.5)*x)),
    background=rpois(n,10))
G <- rbind(gA=50*x, nt1=50*rbinom(n,1L,.3), nt2=50*rbinom(n,1L,.3))
colnames(Y) <- colnames(G) <- paste0("c",seq_len(n))
o <- import_data(Y,G,data.frame(grna_id=rownames(G),
    grna_target=c("A","non-targeting","non-targeting")),moi="high",
    extra_covariates=data.frame(z=z,row.names=colnames(Y)))
o <- set_analysis_parameters(o,discovery_pairs=data.frame(grna_target="A",
    response_id=c("null","down")),formula_object=~z,
    resampling_mechanism="crt",resampling_approximation="no_approximation")
o <- assign_grnas(o,method="thresholding",threshold=5L,print_progress=FALSE)
o <- run_qc(o,response_n_umis_range=c(0,1),response_n_nonzero_range=c(0,1),
    n_nonzero_trt_thresh=7L,n_nonzero_cntrl_thresh=7L)
o <- set_perturbation_fit_method(o,"fast_logistic")
fit <- withr::with_seed(71020,run_discovery_analysis(o,parallel=FALSE,print_progress=FALSE))
print(get_result(fit, "run_discovery_analysis"))
print(lapply(fit@fits@perturbation_fits,function(x)x$diagnostics))
# A cache hit retains method/fallback diagnostics and reconstructs exactly.
warm <- withr::with_seed(71020,run_discovery_analysis(fit,parallel=FALSE,print_progress=FALSE))
stopifnot(identical(get_result(fit, "run_discovery_analysis"),get_result(warm, "run_discovery_analysis")))
# Restore the reference selector explicitly for distributed/pipeline execution.
reference <- set_perturbation_fit_method(fit,"glm.fit")
