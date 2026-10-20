set.seed(20261019)
glm_probe <- stats::glm.fit
old <- body(glm_probe)
body(glm_probe) <- substitute({
    answer <- BODY
    answer$pre_mask_start <- start
    answer
}, list(BODY=old))
byte_equal <- function(a,b) identical(writeBin(as.double(a),raw()),writeBin(as.double(b),raw()))
probe <- function(name,x,y,control=list()) {
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(glm_probe(x=x,y=y,family=stats::binomial(),control=control),
        warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
    if(inherits(fit,'error')) return(list(name=name,error=conditionMessage(fit)))
    beta <- fit$coefficients
    beta[is.na(beta)] <- 0
    eta <- as.vector(x %*% beta) + rep.int(0,nrow(x))
    mu <- stats::binomial()$linkinv(eta)
    plain <- stats::plogis(eta)
    list(name=name,n=nrow(x),p=ncol(x),rank=fit$rank,iter=fit$iter,converged=fit$converged,
        boundary=fit$boundary,missing=sum(is.na(fit$coefficients)),warnings=warnings,
        coefficient_bytes=byte_equal(beta,fit$pre_mask_start),
        eta_bytes=byte_equal(eta,fit$linear.predictors),mu_bytes=byte_equal(mu,fit$fitted.values),
        plogis_bytes=byte_equal(plain,fit$fitted.values),
        max_eta_difference=max(abs(eta-fit$linear.predictors)),
        max_mu_difference=max(abs(mu-fit$fitted.values)),
        plogis_difference=max(abs(plain-fit$fitted.values)),
        mu_range=range(fit$fitted.values),coef_size=as.numeric(object.size(beta)),
        probability_size=as.numeric(object.size(fit$fitted.values)))
}
n<-1000L
z<-stats::rnorm(n);q<-stats::rnorm(n);x<-cbind('(Intercept)'=1,z=z,q=q)
y<-stats::rbinom(n,1,stats::plogis(-.2+.5*z-.3*q))
cases<-list(
    ordinary=probe('ordinary',x,y),
    intercept=probe('intercept',x[,1,drop=FALSE],y),
    allzero=probe('allzero',x,rep(0L,n)),
    allone=probe('allone',x,rep(1L,n)),
    separated=probe('separated',x,as.integer(z>0)),
    duplicate=probe('duplicate',cbind(x,z_again=z),y),
    collinear=probe('collinear',cbind(x,zq=2*z+3*q),y),
    zerocol=probe('zerocol',cbind(x,empty=0),y),
    rankone=probe('rankone',cbind(a=1,b=1,c=0),rep(c(0,1),length.out=1)),
    nonconverged=probe('nonconverged',x,y,control=list(maxit=1L)),
    empty=probe('empty',matrix(numeric(),n,0),y))
# Use the actual default-control regime, with column scaling and nearly aliased
# designs. This does not call any package fitter or CRT sampler.
for(i in seq_len(160L)) {
    scale<-10^sample(c(-150,-40,-10,0,10,40,150),1)
    a<-stats::rnorm(80);b<-stats::rnorm(80)
    xx<-cbind(1,a*scale,b, a*scale + b*10^sample(c(-16,-12,-8,0),1))
    yy<-if(i%%3L==0L) as.integer(a>0) else stats::rbinom(80,1,0.3)
    cases[[paste0('stress',i)]]<-probe(paste0('stress',i),xx,yy)
}
saveRDS(cases,'/private/tmp/sceptre-crt-compact-probe.rds')
jsonlite::write_json(cases,'/private/tmp/sceptre-crt-compact-probe.json',pretty=TRUE,auto_unbox=TRUE,digits=17)
base<-cases[1:11]
for(a in base)cat(a$name,if(!is.null(a$error))a$error else paste('rank',a$rank,'converged',a$converged,'boundary',a$boundary,'coef exact',a$coefficient_bytes,'eta exact',a$eta_bytes,'mu exact',a$mu_bytes,'plogis exact',a$plogis_bytes,'warnings',paste(a$warnings,collapse=' | ')), '\n')
success<-Filter(function(x)is.null(x$error),cases)
cat('total',length(cases),'successful',length(success),'boundary',sum(vapply(success,`[[`,logical(1),'boundary')),'rankdef',sum(vapply(success,function(x)x$rank<x$p,logical(1))),'nonconverged',sum(!vapply(success,`[[`,logical(1),'converged')),'mu-byte-failures',sum(!vapply(success,`[[`,logical(1),'mu_bytes')),'plogis-byte-failures',sum(!vapply(success,`[[`,logical(1),'plogis_bytes')),'\n')
