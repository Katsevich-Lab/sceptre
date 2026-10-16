args <- commandArgs(TRUE)
root <- args[[1L]]
.libPaths(c("/private/tmp/sceptre-resource-reference-release-library", "/private/tmp/sceptre-rf0-ondisc-library", .libPaths()))
library(sceptre)
Rcpp::sourceCpp(file.path(root,"probe.cpp"), rebuild=TRUE, showOutput=FALSE)
cases <- list(
 ordinary=list(y=c(0,1,5),mu=c(1,1,1)),
 exact_means=list(y=c(1,2,3),mu=c(1,2,3)),
 all_zero_response=list(y=c(0,0,0),mu=c(1,1,1)),
 zero_means=list(y=c(1,1,1),mu=c(0,0,0)),
 tiny_means=list(y=c(1,2,3),mu=rep(1e-200,3)),
 huge_response=list(y=c(1,2,1e300),mu=c(1,1,1)),
 tiny_response=list(y=c(1e-300,2e-300,3e-300),mu=c(1,1,1)))
original <- getFromNamespace("estimate_theta","sceptre")
observed <- lapply(cases,function(z) {
 stopifnot(length(z$y)==length(z$mu))
 result <- probe_estimate(z$y,z$mu,2)
 baseline <- original(z$y,z$mu,2,50,.Machine$double.eps^(1/4))
 stopifnot(identical(result$result,baseline))
 c(list(input=z),result,list(original_exact=TRUE))
})
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(271828)
for(i in seq_len(100)) {
 y <- sample.int(20,8,replace=TRUE)-1
 mu <- 10^runif(8,-200,200)
 got <- probe_estimate(y,mu,6)
 baseline <- original(y,mu,6,50,.Machine$double.eps^(1/4))
 stopifnot(identical(got$result,baseline))
 observed[[paste0("bounded_grid_",i)]] <- c(list(input=list(y=y,mu=mu,dfr=6)),got,list(original_exact=TRUE))
}
xs <- c(0,.Machine$double.xmin,.Machine$double.xmin*.Machine$double.eps,1e-310,1,Inf,NaN)
special <- lapply(xs,function(x) list(x=x,digamma=probe_special(x,FALSE),trigamma=probe_special(x,TRUE)))
report <- list(source_commit="21f9ea098b69c444a884a49bc254de711968e3b5", probe_note="Entry renamed/export suppressed and catch instrumented only; arithmetic/control flow unchanged. Special function probes are diagnostics, not new production fallbacks.", RNGkind=RNGkind(), R=R.version.string, platform=R.version$platform, BH=as.character(packageVersion("BH")), Rcpp=as.character(packageVersion("Rcpp")), original_library=find.package("sceptre"), cases=observed, special=special)
saveRDS(report,file.path(root,"report.rds"))
jsonlite::write_json(report,file.path(root,"report.json"),auto_unbox=TRUE,pretty=TRUE,na="null",digits=17)
cat("case_count",length(observed),"all_original_exact",all(vapply(observed,`[[`,logical(1),"original_exact")),"\n")
print(table(vapply(observed,`[[`,character(1),"type")))
print(lapply(special,function(z) c(x=z$x,digamma=z$digamma$type,trigamma=z$trigamma$type)))
