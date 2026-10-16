.libPaths(c('/private/tmp/sceptre-resource-reference-release-library', '/private/tmp/sceptre-rf0-ondisc-library', .libPaths()))
suppressPackageStartupMessages(library(sceptre))
data.table::setDTthreads(1)
make_raw <- function(i,p,x,dims) {
  out <- Matrix::sparseMatrix(i=1,p=c(0,1),x=1,repr='C')
  out@i <- as.integer(i);out@p <- as.integer(p);out@x <- as.numeric(x);out@Dim <- as.integer(dims)
  out
}
for(label in c('valid_unique','sorted_duplicate','unsorted_unique')) {
  x <- switch(label, valid_unique=make_raw(c(0,1),c(0,2,2),c(2,3),c(2,2)), sorted_duplicate=make_raw(c(0,0),c(0,2,2),c(2,3),c(2,2)), unsorted_unique=make_raw(c(1,0),c(0,2,2),c(3,2),c(2,2)))
  cat(label, ': ', paste(methods::validObject(x,test=TRUE),collapse='; '),'\n',sep='')
}
t <- methods::new('dgTMatrix',i=c(0L,0L,1L),j=c(0L,0L,0L),x=c(2,3,7),Dim=c(2L,2L),Dimnames=list(c('g1','g2'),c('c2','c1')))
cat('duplicate_triplet_valid:',methods::validObject(t,test=TRUE),'\n')
print(as.matrix(t))
metadata <- data.frame(grna_id=c('g2','g1','g1'),grna_target=c('target','target','target'),vector_id=c('v2','v1','v1'))
print(sceptre:::collapse_grna_matrix(t,metadata))
