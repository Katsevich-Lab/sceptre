# Run with R_LIBS pointing at an isolated release installation. Acquire the
# shared heavy lease before invoking. No user library or dataset is modified.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
out <- args[[1L]]
harness <- normalizePath(args[[2L]])
if (dir.exists(out)) stop("Use a new evidence directory")
dir.create(out, recursive = TRUE)
library(sceptre)
data.table::setDTthreads(1L)
f <- function(name) getFromNamespace(name, "sceptre")
source(harness)
files <- c("DESCRIPTION", "NAMESPACE", list.files("R", "[.]R$", full.names=TRUE),
    list.files("src", "[.](cpp|h)$", full.names=TRUE))
manifest <- list(started=as.character(Sys.time()), source_revision=system2("git",c("rev-parse","HEAD"),stdout=TRUE),
    source_sha256=setNames(vapply(files,digest::digest,character(1),file=TRUE,algo="sha256"),files),
    harness=harness,harness_sha256=digest::digest(harness,file=TRUE,algo="sha256"),
    installed=find.package("sceptre"),RcppEigen=as.character(packageVersion("RcppEigen")),
    session=sessionInfo(),interpretation="Experimental pilot; independent cell datasets are replicates; no pooled-endpoint binomial precision")
inst <- list.files(find.package("sceptre"),recursive=TRUE,full.names=TRUE)
manifest$installed_sha256 <- setNames(vapply(inst,digest::digest,character(1),file=TRUE,algo="sha256"),inst)
saveRDS(manifest,file.path(out,"manifest.rds"))

# Paired actual-count public-workflow tests. W03 owns the generator and its truth.
# These are real tests of regenerated counts, never synthetic Uniform p-values.
cases <- validation_cases()
cases <- cases[cases$id %in% c("high_crt_empirical","high_crt_skew","confounded_crt"),]
rows <- list()
for (i in seq_len(nrow(cases))) {
    case <- cases[i,,drop=FALSE]
    j <- match(case$id,validation_cases()$id)
    for (r in seq_len(40L)) {
        for (method in c("glm.fit","fast_logistic")) {
            options(sceptre.validation.configure=local({m<-method;function(object,case) set_perturbation_fit_method(object,m)}))
            x <- validation_replicate(case,r,20261000L+1000L*j+2L*r,600L)
            x$method <- method
            rows[[length(rows)+1L]] <- x
        }
    }
    cat(case$id,"complete\n")
}
options(sceptre.validation.configure=NULL)
rows <- do.call(rbind,rows)
write.csv(rows,file.path(out,"paired-science.csv"),row.names=FALSE)
summary <- do.call(rbind,lapply(split(rows,rows$method),function(x){
    a<-validation_summary(x);a$method<-x$method[1];a
}))
write.csv(summary,file.path(out,"science-summary.csv"),row.names=FALSE)

# Runtime includes preparation unless explicitly labelled prepared or warm.
# Rprofmem measures allocation volume, not peak RSS or resident memory.
measure <- function(fun,reps=20L) {
    elapsed <- replicate(reps,unname(system.time(fun())[["elapsed"]]))
    log <- tempfile();Rprofmem(log);value<-fun();Rprofmem(NULL)
    allocation <- suppressWarnings(as.numeric(sub(" .*","",readLines(log,warn=FALSE))))
    unlink(log)
    list(median_seconds=median(elapsed),min_seconds=min(elapsed),
        allocated_bytes=sum(allocation,na.rm=TRUE),value=value)
}
bench <- list()
set.seed(71019)
for (n in c(5000L,20000L)) {
    D<-cbind('(Intercept)'=1,matrix(rnorm(n*5),n));colnames(D)[-1]<-paste0('z',1:5)
    idx<-which(rbinom(n,1,plogis(D%*%c(-1.2,-.3,-.1,.1,.2,.3)))==1L)
    prepared<-f("prepare_grna_logistic_design")(D)
    ref<-f("fit_perturbation_model")(idx,D)
    fast<-f("fit_perturbation_model")(idx,D,"fast_logistic",prepared)
    tight<-glm.fit(D,as.integer(seq_len(n)%in%idx),family=binomial(),
        control=glm.control(epsilon=1e-12,maxit=100L))
    lanes<-list(reference_cold=function()f("fit_perturbation_model")(idx,D),
        candidate_cold=function()f("fit_perturbation_model")(idx,D,"fast_logistic"),
        candidate_prepared=function()f("fit_perturbation_model")(idx,D,"fast_logistic",prepared),
        candidate_warm=function()f("reconstruct_perturbation_probabilities")(fast$fit,D))
    for (label in names(lanes)) {
        z<-measure(lanes[[label]])
        bench[[length(bench)+1L]]<-data.frame(n=n,p=ncol(D),lane=label,
            median_seconds=z$median_seconds,min_seconds=z$min_seconds,
            allocated_bytes=z$allocated_bytes,
            max_probability_error=max(abs(fast$probabilities-tight$fitted.values)),
            max_coefficient_error=max(abs(fast$fit$fitted_coefs-tight$coefficients)),
            method_used=fast$fit$diagnostics$method_used,
            candidate_iterations=fast$fit$diagnostics$iterations,
            candidate_cache_bytes=as.numeric(object.size(fast$fit)),
            reference_cache_bytes=as.numeric(object.size(ref$fit)))
    }
}
write.csv(do.call(rbind,bench),file.path(out,"fit-benchmark.csv"),row.names=FALSE)

# Real historical schema-three physical layouts, preserved read-only.
old_path<-"/private/tmp/sceptre-cleaned-voice-ci-20261019/direct/high_crt/discovery_without_calibration/sceptre_object.rds"
if(file.exists(old_path)) {
    old<-readRDS(old_path);before<-serialize(old,NULL)
    migrated<-f("migrate_sceptre_object")(old)
    saveRDS(list(path=old_path,sha256=digest::digest(old_path,file=TRUE,algo="sha256"),
        original_has_selector=!is.null(attr(old@spec,"perturbation_fit_method",exact=TRUE)),
        default=get_perturbation_fit_method(migrated),valid=validObject(migrated),
        result_current=f("result_record_is_current")(migrated,"run_discovery_analysis"),
        original_unchanged=identical(serialize(old,NULL),before),
        response_fit_count=length(migrated@fits@response_fits),
        perturbation_fit_count=length(migrated@fits@perturbation_fits)),
        file.path(out,"historical-migration.rds"))
}
# Actual DC-TAP counts, read-only, evenly spaced capped cells. These tests compare
# numerical outputs only; the real tested pairs have no asserted null/effect truth.
real_root <- "/Users/ekatsevi/data/projects/sceptre3/perturbseq-survey/dctap_k562_highmoi_GSE303901/sceptre"
real_paths <- file.path(real_root,c("response_matrix.rds","grna_matrix_aligned.rds","grna_target_data_frame.csv"))
if(all(file.exists(real_paths))) {
    Y <- readRDS(real_paths[1]); G <- readRDS(real_paths[2])
    meta <- read.csv(real_paths[3])
    stopifnot(identical(colnames(Y),colnames(G)),identical(rownames(G),meta$grna_id))
    # Entire observed targeted panel, not a selected-gene depth or a claim of
    # whole-transcriptome depth. No gene subset is taken before these totals.
    response_total <- Matrix::colSums(Y)
    guide_total <- Matrix::colSums(G)
    cells <- unique(as.integer(round(seq(1,ncol(Y),length.out=min(3000L,ncol(Y))))))
    D <- cbind('(Intercept)'=1,log_panel_umis=log1p(response_total[cells]),
        log_guide_umis=log1p(guide_total[cells]))
    chosen <- head(which(meta$grna_target!="non-targeting"),6L)
    genes <- head(order(Matrix::rowSums(Y>0),decreasing=TRUE),2L)
    comparisons <- list()
    for(g in chosen) {
        idx <- which(as.numeric(G[g,cells])>=5)
        a <- f("fit_perturbation_model")(idx,D)
        b <- f("fit_perturbation_model")(idx,D,"fast_logistic")
        for(j in genes) {
            y <- as.numeric(Y[j,cells])
            fit <- f("fit_response")(y,D)
            score <- f("prepare_score_components")(y,D,fit)
            values <- lapply(list(a$probabilities,b$probabilities),function(prob) withr::with_seed(71025,{
                bank<-f("sample_crt_indices")(prob,499L)
                f("run_association_test")(y,score$mu,score$a,score$w,score$D,idx,length(idx),TRUE,
                    bank,499L,0L,0L,FALSE,FALSE,0L)
            }))
            comparisons[[length(comparisons)+1L]]<-data.frame(guide=meta$grna_id[g],response=rownames(Y)[j],
                n=length(cells),n_treated=length(idx),
                max_probability_difference=max(abs(a$probabilities-b$probabilities)),
                max_coefficient_difference=if(is.null(a$fit)||is.null(b$fit)) NA_real_ else
                    max(abs(a$fit$fitted_coefs-b$fit$fitted_coefs)),
                method_used=if(is.null(b$fit))"uncacheable_reference_fallback" else b$fit$diagnostics$method_used,
                reference_p=values[[1]]$p,candidate_p=values[[2]]$p,
                reference_fc=values[[1]]$fc,candidate_fc=values[[2]]$fc)
        }
    }
    write.csv(do.call(rbind,comparisons),file.path(out,"real-dctap-comparison.csv"),row.names=FALSE)
    saveRDS(list(paths=real_paths,sha256=setNames(vapply(real_paths,digest::digest,character(1),file=TRUE,algo="sha256"),real_paths),
        input_dimensions=list(response=dim(Y),guide=dim(G)),cells=cells,design_columns=colnames(D),
        guide_threshold=5L,interpretation="Illustrative assigned groups and numerical comparisons; no asserted biological discovery or real null/power truth"),
        file.path(out,"real-dctap-manifest.rds"))
}

manifest$finished<-as.character(Sys.time())
saveRDS(manifest,file.path(out,"manifest.rds"))
capture.output(str(manifest),file=file.path(out,"manifest.txt"))
cat("W07 bounded pilot complete; see CSVs and immutable source/runtime hashes.\n")
