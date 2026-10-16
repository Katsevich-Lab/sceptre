args <- commandArgs(TRUE)
.libPaths(c(args[[1L]], '/private/tmp/sceptre-rf0-ondisc-library', .libPaths()))
suppressPackageStartupMessages(library(sceptre))
convert <- get('set_matrix_accessibility', asNamespace('sceptre'), inherits=TRUE)
matrices <- list(symmetric=matrix(c(1,2,2,3),2L,2L),
    near_symmetric=matrix(c(1,2+4*.Machine$double.eps,2,3),2L,2L))
results <- lapply(matrices,function(input){
 dimnames(input) <- list(c('c1','c2'),c('c1','c2'))
 before <- serialize(input,NULL)
 direct <- lapply(c(TRUE,FALSE),function(row_access){
  actual <- convert(input,row_access)
  stopifnot(methods::validObject(actual))
  list(row_access=row_access,actual=as.matrix(actual),
       exact=identical(as.matrix(actual),input),class=class(actual))
 })
 guides <- matrix(c(8,9),1L,2L,dimnames=list('g',c('c1','c2')))
 object <- suppressMessages(sceptre::import_data(input,guides,
   data.frame(grna_id='g',grna_target='non-targeting'),'low'))
 public <- as.matrix(sceptre::get_response_matrix(object))
 stopifnot(identical(serialize(input,NULL),before))
 list(expected=input,Matrix_detects_symmetric=isSymmetric(input),
   constructor_class=class(methods::as(input,'TsparseMatrix')),
   direct=direct,public=public,public_exact=identical(public,input),
   expected_column_totals=colSums(input),actual_import_totals=object@covariate_data_frame$response_n_umis)
})
ns <- find.package('sceptre')
report <- list(source_commit=args[[2L]],package=ns,
 dll_md5=tools::md5sum(file.path(ns,'libs',paste0('sceptre',.Platform$dynlib.ext))),
 Matrix=as.character(utils::packageVersion('Matrix')),results=results,session=utils::sessionInfo())
saveRDS(report,args[[3L]])
print(lapply(results,function(x)x[c('expected','constructor_class','public','public_exact','expected_column_totals','actual_import_totals')]))
