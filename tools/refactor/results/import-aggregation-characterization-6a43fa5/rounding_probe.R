.libPaths(c('/private/tmp/sceptre-resource-reference-release-library','/private/tmp/sceptre-rf0-ondisc-library',.libPaths()))
suppressPackageStartupMessages(library(sceptre)); data.table::setDTthreads(1)
proto <- function(m,d) {
  ids <- unique(d$vector_id); n <- ncol(m)
  t(vapply(ids,function(id) Matrix::colSums(m[d$grna_id[d$vector_id==id],,drop=FALSE]),numeric(n)))
}
meta <- data.frame(grna_id=c('a','b','c'),grna_target='t',vector_id='v')
for(vals in list(c(2^53,1,1),c(.1,.2,.3),c(1e16,1,-1e16),c(1e308,1e308,0))) {
 m <- matrix(vals,3,1,dimnames=list(c('a','b','c'),'cell'))
 sparse <- methods::as(Matrix::Matrix(m,sparse=TRUE),'CsparseMatrix')
 ref <- sceptre:::collapse_grna_matrix(m,meta); alt <- proto(sparse,meta)
 cat('values:',paste(sprintf('%a',vals),collapse=','),'\n')
 cat('original:',sprintf('%a',ref[1,1]),'sparse:',sprintf('%a',alt[1,1]),'exact:',identical(ref,alt),'all.equal:',paste(all.equal(ref,alt,tolerance=1e-8),collapse='; '),'\n')
}
