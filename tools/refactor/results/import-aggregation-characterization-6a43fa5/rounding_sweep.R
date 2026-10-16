.libPaths(c('/private/tmp/sceptre-resource-reference-release-library','/private/tmp/sceptre-rf0-ondisc-library',.libPaths()))
suppressPackageStartupMessages(library(sceptre));data.table::setDTthreads(1);set.seed(94201)
probe<-function(m,d)t(vapply(unique(d$vector_id),function(id)Matrix::colSums(m[d$grna_id[d$vector_id==id],,drop=FALSE]),numeric(ncol(m))))
mismatches<-list();outside<-list()
for(k in seq_len(200L)) {
 i<-sample.int(4L,30L,TRUE)-1L;j<-sample.int(3L,30L,TRUE)-1L
 x<-sample(c(0,.1,.2,.3,1,2,2^53,1e-20),30L,TRUE)
 m<-methods::new('dgTMatrix',i=i,j=j,x=x,Dim=c(4L,3L),Dimnames=list(paste0('g',1:4),c('z','a','m')))
 d<-data.frame(grna_id=paste0('g',sample.int(4L,8L,TRUE)),grna_target='t',vector_id=sample(c('vZ','vA'),8L,TRUE))
 ref_input<-unserialize(serialize(list(m=m,d=d),NULL));alt_input<-unserialize(serialize(list(m=m,d=d),NULL));before_ref<-serialize(ref_input,NULL);before_alt<-serialize(alt_input,NULL)
 ref<-sceptre:::collapse_grna_matrix(ref_input$m,ref_input$d);alt<-probe(alt_input$m,alt_input$d)
 stopifnot(identical(before_ref,serialize(ref_input,NULL)),identical(before_alt,serialize(alt_input,NULL)))
 if(!identical(ref,alt)) {
  z<-list(case=k,max_abs_diff=max(abs(ref-alt)),all_equal=all.equal(ref,alt,tolerance=1e-8),input=m,metadata=d,original=ref,probe=alt)
  mismatches[[length(mismatches)+1L]]<-z
  if(!isTRUE(z$all_equal))outside[[length(outside)+1L]]<-z
 }
}
saveRDS(list(seed=94201L,cases=200L,mismatches=mismatches,outside=outside),'/private/tmp/sceptre-import-aggregation-characterization-6a43fa5/rounding-sweep.rds')
cat('200 valid small duplicate-triplet cases; exact mismatches',length(mismatches),'outside 1e-8',length(outside),'\n')
if(length(mismatches)) {print(mismatches[[1]][c('case','max_abs_diff','all_equal','metadata','original','probe')])}
