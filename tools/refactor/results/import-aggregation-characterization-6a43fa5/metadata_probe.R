.libPaths(c('/private/tmp/sceptre-resource-reference-release-library','/private/tmp/sceptre-rf0-ondisc-library',.libPaths()));suppressPackageStartupMessages(library(sceptre));data.table::setDTthreads(1)
f<-readRDS('/private/tmp/sceptre-import-aggregation-characterization-6a43fa5/fixtures.rds')$numeric_C
probe<-function(m,d)t(vapply(unique(d$vector_id),function(id){ members<-d |> dplyr::filter(vector_id==id) |> dplyr::pull(grna_id);Matrix::colSums(m[members,,drop=FALSE])},numeric(ncol(m))))
for(label in c('grna_factor','vector_na')) {
 d<-f$d
 if(label=='grna_factor')d$grna_id<-factor(d$grna_id)
 if(label=='vector_na')d$vector_id[c(1,4)]<-NA_character_
 old<-sceptre:::collapse_grna_matrix(unserialize(serialize(f$m,NULL)),d)
 new<-probe(unserialize(serialize(f$m,NULL)),d)
 cat(label,'exact',identical(old,new),'\n');print(old);print(new)
}
