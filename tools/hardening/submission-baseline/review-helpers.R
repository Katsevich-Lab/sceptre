options(sceptre.refactor.fixture_dir='/Users/ekatsevi/.codex/worktrees/ab6a/sceptre/tests/testthat/fixtures/refactor')
library(sceptre)
source('/Users/ekatsevi/.codex/worktrees/ab6a/sceptre/tests/testthat/helper-refactor.R')
ns<-asNamespace('sceptre')
load_functions<-function(path){
 e<-new.env(parent=ns)
 for(p in sort(list.files(file.path(path,'R'),pattern='\\.[Rr]$',full.names=TRUE))) {
  for(x in parse(p,keep.source=FALSE)) if(is.call(x)&&identical(x[[1]],as.name('<-'))&&length(x)==3L&&is.symbol(x[[2]])&&is.call(x[[3]])&&identical(x[[3]][[1]],as.name('function'))) eval(x,e)
 }
 e
}
old<-load_functions('/private/tmp/sceptre-submission-matched-check/qualified/sceptre')
new<-load_functions('/Users/ekatsevi/.codex/worktrees/ab6a/sceptre')
assert<-function(value,label){stopifnot(isTRUE(value));cat('PASS',label,'\n')}
assert(identical(old$ondisc_required_api(),new$ondisc_required_api()),'ondisc ABI argument names and order')
inputs<-readRDS(file.path(refactor_fixture_dir(),'inputs.rds'))
cases<-readRDS(file.path(refactor_fixture_dir(),'cases.rds'))
for(case in c('low_complement_perm','low_nt_perm','high_complement_crt')) {
 object<-refactor_capture(refactor_prepare(inputs,cases[[case]]))$value
 before<-serialize(object,NULL)
 for(stage in c('design','assignment','cell_qc','pair_qc')) assert(identical(old$stage_input_identity(object,stage),new$stage_input_identity(object,stage)),paste(case,stage,'raw stage identity'))
 old_plan<-old$make_run_plan(object,'run_discovery_analysis')
 new_plan<-new$make_run_plan(object,'run_discovery_analysis')
 assert(identical(old_plan,new_plan),paste(case,'realized run plan'))
 old_context<-old$prepare_analysis_context(object,'run_discovery_analysis',dplyr::filter(old_plan@pairs,pass_qc),old_plan)
 new_context<-new$prepare_analysis_context(object,'run_discovery_analysis',dplyr::filter(new_plan@pairs,pass_qc),new_plan)
 assert(identical(old_context,new_context),paste(case,'ordered context and numerical settings'))
 set.seed(3121);rng<-.Random.seed
 m1<-old$prepare_execution(object,'run_discovery_analysis',3L,2L)
 set.seed(3121);m2<-new$prepare_execution(object,'run_discovery_analysis',3L,2L)
 assert(identical(m1,m2),paste(case,'same-runtime manifest payload'))
 p1<-lapply(m1$units,function(unit) old$execute_work_unit(object,m1,unit$id))
 p2<-lapply(m2$units,function(unit) new$execute_work_unit(object,m2,unit$id))
 assert(identical(p1,p2),paste(case,'unit raw tables fits and certificates'))
 a<-old$collect_execution(object,m1,rev(p1));b<-new$collect_execution(object,m2,rev(p2))
 assert(identical(a,b),paste(case,'complete collected value and family correction'))
 assert(identical(.Random.seed,rng),paste(case,'R RNG unchanged by execution'))
 assert(identical(serialize(object,NULL),before),paste(case,'caller value unchanged'))
 err<-function(expr)tryCatch({force(expr);NULL},error=function(e)list(class=class(e),message=conditionMessage(e)))
 assert(identical(err(old$collect_execution(object,m1,p1[-1L])),err(new$collect_execution(object,m2,p2[-1L]))),paste(case,'missing proposal error'))
 assert(identical(err(old$execute_work_unit(object,m1,'missing')),err(new$execute_work_unit(object,m2,'missing'))),paste(case,'unknown unit error'))
}
cat('Review scope: isolated source-function environments over the same installed pre-cleanup class/runtime; this does not replace final installed candidate checks.\n')
