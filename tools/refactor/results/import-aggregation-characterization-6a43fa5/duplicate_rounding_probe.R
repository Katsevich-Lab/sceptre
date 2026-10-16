source('/private/tmp/sceptre-import-aggregation-characterization-6a43fa5/rounding_probe.R')
m <- methods::new('dgTMatrix',i=c(0L,1L,0L),j=c(0L,0L,0L),x=c(1,2^53,1),Dim=c(2L,1L),Dimnames=list(c('a','b'),'cell'))
d <- data.frame(grna_id=c('a','b'),grna_target='t',vector_id='v')
for(input in list(m,methods::as(m,'CsparseMatrix'),methods::as(m,'RsparseMatrix'))) {
 ref <- sceptre:::collapse_grna_matrix(input,d); alt <- proto(input,d)
 cat('class',class(input),'dense input values',paste(sprintf('%a',as.matrix(input)),collapse=','),'ref',sprintf('%a',ref[1,1]),'proto',sprintf('%a',alt[1,1]),'exact',identical(ref,alt),'\n')
}
