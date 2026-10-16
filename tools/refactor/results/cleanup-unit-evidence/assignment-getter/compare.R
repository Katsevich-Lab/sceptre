.libPaths(c('/private/tmp/sceptre-cl3-library-627e7e4','/private/tmp/sceptre-rf0-ondisc-library',.libPaths()))
suppressPackageStartupMessages(library(sceptre))
old<-new.env(parent=asNamespace('sceptre'));new<-new.env(parent=asNamespace('sceptre'))
sys.source('/private/tmp/sceptre-assignment-getter-probe/original.R',old)
sys.source('/Users/ekatsevi/code/research/sceptre/R/getters_and_setters.R',new)
count<-0L
for(assignments in list(list(single=c(2L,5L)),list(first=c(2L,5L),second=c(1L,4L)),list(first=integer(),second=c(1L,4L)))) {
  object<-methods::new('sceptre_object',grna_matrix=list(matrix(0,length(assignments),5L)),initial_grna_assignment_list=assignments,functs_called=c(assign_grnas=TRUE,run_qc=TRUE),cells_in_use=c(2L,4L,5L))
  before<-serialize(object,NULL)
  for(qc in c(FALSE,TRUE)) {
    a<-old$get_grna_assignments(object,qc);b<-new$get_grna_assignments(object,qc)
    stopifnot(identical(a,b),identical(serialize(object,NULL),before))
    count<-count+1L
  }
}
cat(count,'exact original/candidate getter cases; caller objects unchanged.\n')
