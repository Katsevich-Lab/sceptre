.libPaths(c('/private/tmp/sceptre-cl3-library-627e7e4','/private/tmp/sceptre-rf0-ondisc-library',.libPaths()))
suppressPackageStartupMessages(library(sceptre))
old<-new.env(parent=asNamespace('sceptre'));new<-new.env(parent=asNamespace('sceptre'))
sys.source('/private/tmp/sceptre-b2-pair-probe/original.R',old)
sys.source('/private/tmp/sceptre-b2-pair-probe/source/R/pair_constructor_functs.R',new)
response<-matrix(1,3,2,dimnames=list(c('left','right','other'),NULL))
object<-methods::new('sceptre_object',response_matrix=list(response),grna_target_data_frame=data.frame(grna_target=c('B','non-targeting','A','B')))
pc_frame<-data.frame(grna_target=c('B','A','B'),response_id=c('left','other','left'))
count<-0L
for(as_table in c(FALSE,TRUE)) for(mode in c('none','pc_pairs','pairs_containing_pc_targets')) {
  pc<-if(as_table)data.table::as.data.table(pc_frame) else pc_frame
  before<-serialize(pc,NULL)
  # Independent copies prevent the original's known reference mutation from
  # changing candidate inputs or later comparisons.
  a<-old$construct_trans_pairs(object,data.table::copy(pc),mode)
  b<-new$construct_trans_pairs(object,pc,mode)
  stopifnot(identical(a,b),identical(serialize(pc,NULL),before))
  count<-count+1L
}
cat(count,'exact complete original/candidate trans outputs including expand.grid metadata.\n')
