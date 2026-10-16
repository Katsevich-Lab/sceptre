# Bounded, read-only source characterization. All writes are under this directory.
root <- '/private/tmp/sceptre-import-aggregation-characterization-6a43fa5'
.libPaths(c('/private/tmp/sceptre-resource-reference-release-library',
            '/private/tmp/sceptre-rf0-ondisc-library', .libPaths()))
suppressPackageStartupMessages(library(sceptre))
stopifnot(normalizePath(find.package('sceptre')) == '/private/tmp/sceptre-resource-reference-release-library/sceptre')
stopifnot(requireNamespace('ondisc',quietly=TRUE), packageVersion('ondisc')=='1.3.5')
stopifnot(requireNamespace('digest',quietly=TRUE),requireNamespace('jsonlite',quietly=TRUE))
data.table::setDTthreads(1)
set.seed(1847)
sha <- function(x) digest::digest(x,algo='sha256',serialize=TRUE)
file_sha <- function(p) digest::digest(file=p,algo='sha256')
clone <- function(x) unserialize(serialize(x,NULL))
ns <- asNamespace('sceptre')
original <- function(name) get(name,envir=ns,inherits=FALSE)
number <- function(x) list(decimal=format(x,digits=17,scientific=FALSE,trim=TRUE),hex=sprintf('%a',as.numeric(x)))
snapshot <- function(x,dense=TRUE) {
 out <- list(class=class(x),typeof=typeof(x),dim=dim(x),dimnames=dimnames(x),sha256=sha(x))
 if(isS4(x)) {
  out$slots <- lapply(intersect(c('i','j','p','x'),slotNames(x)),function(n) {v<-slot(x,n);if(is.double(v))number(v) else v})
  names(out$slots)<-intersect(c('i','j','p','x'),slotNames(x))
  out$validity<-methods::validObject(x,test=TRUE)
 }
 if(dense) {m<-as.matrix(x);out$column_major_values<-number(as.numeric(m));out$storage_mode<-storage.mode(m)}
 out
}
capture <- function(expr) {
 warnings<-messages<-character();error<-NULL
 value<-tryCatch(withCallingHandlers(expr,warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')},message=function(m){messages<<-c(messages,conditionMessage(m));invokeRestart('muffleMessage')}),error=function(e){error<<-conditionMessage(e);NULL})
 list(value=value,warnings=warnings,messages=messages,error=error)
}
# A scratch characterization prototype, NOT a candidate package implementation.
# Matrix's sparse row subset/colSums keeps metadata row repetitions and their order;
# t(vapply()) deliberately retains the ORIGINAL dense output simplification rules.
sparse_probe <- function(grna_matrix,grna_target_data_frame) {
 ids<-unique(grna_target_data_frame$vector_id)
 t(vapply(ids,function(id) {
  members<-grna_target_data_frame |> dplyr::filter(vector_id==id) |> dplyr::pull(grna_id)
  Matrix::colSums(grna_matrix[members,,drop=FALSE])
 },numeric(ncol(grna_matrix))))
}
source_env<-new.env(parent=baseenv());sys.source(file.path(root,'original-import_functs.R'),source_env,keep.source=FALSE)
source_functions<-c('create_csc_matrix','collapse_grna_matrix','collapse_grna_target_data_frame','import_data_memory','import_data_from_cellranger_memory')
source_match<-setNames(lapply(source_functions,function(name) {
 a<-original(name);b<-get(name,source_env)
 # Deparse/parse strips source-reference attributes, not executable operations.
 list(formals_identical=identical(formals(a),formals(b)),body_identical=identical(parse(text=deparse(body(a)),keep.source=FALSE),parse(text=deparse(body(b)),keep.source=FALSE)),installed_function_sha256=sha(list(formals(a),deparse(body(a)))))
}),source_functions)
stopifnot(all(vapply(source_match,function(x)x$formals_identical&&x$body_identical,logical(1))))

# B1 malformed layouts are constructed only as raw S4 slots and checked by the
# Matrix validity checker. NO arithmetic, conversion, cbind or native use follows.
raw_csc<-function(i,p,x,dims) {
 out<-Matrix::sparseMatrix(i=1,p=c(0,1),x=1,repr='C')
 out@i<-as.integer(i);out@p<-as.integer(p);out@x<-as.numeric(x);out@Dim<-as.integer(dims)
 out
}
raw_cases<-list(
 sorted_unique=raw_csc(c(0,1),c(0,2,2),c(2,3),c(2,2)),
 sorted_duplicate=raw_csc(c(0,0),c(0,2,2),c(2,3),c(2,2)),
 unsorted_unique=raw_csc(c(1,0),c(0,2,2),c(3,2),c(2,2)),
 ungrouped_cells_pointer=raw_csc(c(0,1,2),c(0,0,1),c(2,3,4),c(3,2)),
 interleaved_gene_offset=raw_csc(c(0,2,4),c(0,3),c(1,2,3),c(3,1)),
 interleaved_guide_offset=raw_csc(c(0,2),c(0,2),c(4,5),c(2,1)))
raw_records<-lapply(raw_cases,snapshot,dense=FALSE)
rm(raw_cases)
valid_dt<-data.table::data.table(cell_idx=c(0L,0L,1L,2L,2L),feature_idx=c(0L,2L,1L,0L,3L),x=c(0,2.5,4,1.2,7))
valid_before<-sha(valid_dt);csc<-original('create_csc_matrix')(clone(valid_dt),3L,4L)
stopifnot(isTRUE(methods::validObject(csc,test=TRUE)),identical(valid_before,sha(valid_dt)))
intended_interleaved<-list(global_feature_ids=1:5,modalities=c('Gene Expression','CRISPR Guide Capture','Gene Expression','CRISPR Guide Capture','Gene Expression'),gene_global_ids=c(1L,3L,5L),guide_global_ids=c(2L,4L),gene_match_rows=match(c(1L,3L,5L),c(1L,3L,5L))-1L,guide_match_rows=match(c(2L,4L),c(2L,4L))-1L)

# B6 fixtures: first-seen vector order differs from guide row order and cell order
# is intentionally nonalphabetical. Identical metadata duplicates count twice;
# a guide also contributes independently to two vectors; an unused row is ignored.
m<-matrix(c(0,.1,2,1.25,2,0,0,3,.125,4,0,2.5,9,8,7),5L,3L,byrow=TRUE,dimnames=list(c('g3','g1','g4','g2','unused'),c('cell-z','cell-a','cell-m')))
d<-data.frame(grna_id=c('g2','g1','g1','g3','g2','g4'),grna_target=c('geneX','geneX','geneX','geneX','geneX','non-targeting'),vector_id=c('vZ','vA','vA','vZ','vA','vB'))
d$annotation<-c('z','a','a','z','a','b')
positions<-which(m!=0,arr.ind=TRUE);ii<-c(positions[,1L]-1L,0L);jj<-c(positions[,2L]-1L,0L);xx<-c(m[positions],0)
perm<-rev(seq_along(xx))
t<-methods::new('dgTMatrix',i=as.integer(ii[perm]),j=as.integer(jj[perm]),x=as.numeric(xx[perm]),Dim=as.integer(dim(m)),Dimnames=dimnames(m))
lt<-methods::new('lgTMatrix',i=t@i,j=t@j,x=t@x>0,Dim=t@Dim,Dimnames=t@Dimnames)
fixtures<-list(numeric_dense=list(m=m,d=d),numeric_C=list(m=methods::as(t,'CsparseMatrix'),d=d),numeric_R=list(m=methods::as(t,'RsparseMatrix'),d=d),numeric_T_unsorted_stored_zero=list(m=t,d=d),logical_dense=list(m=m>0,d=d),logical_C_stored_false=list(m=methods::as(lt,'CsparseMatrix'),d=d),logical_R_stored_false=list(m=methods::as(lt,'RsparseMatrix'),d=d),logical_T_unsorted_stored_false=list(m=lt,d=d),integer_dense=list(m=matrix(as.integer(round(m)),nrow(m),ncol(m),dimnames=dimnames(m)),d=d),one_cell=list(m=t[,1L,drop=FALSE],d=d),zero_cells_helper_only=list(m=m[,integer(),drop=FALSE],d=d),zero_vectors_helper_only=list(m=m,d=d[integer(),,drop=FALSE]),one_vector=list(m=t,d=d[d$vector_id=='vA',,drop=FALSE]))
# Valid triplet duplicates are consumed by Matrix operations before vector sums.
# These are NOT malformed CSCs or repeated file-coordinate imports.
dup_numeric<-methods::new('dgTMatrix',i=c(0L,1L,0L,0L,1L),j=c(0L,0L,0L,1L,1L),x=c(1,2^53,1,0,.25),Dim=c(2L,2L),Dimnames=list(c('g1','g2'),c('cell-z','cell-a')))
dup_logical<-methods::new('lgTMatrix',i=c(0L,1L,0L,1L),j=c(0L,0L,0L,1L),x=c(TRUE,FALSE,TRUE,TRUE),Dim=c(2L,2L),Dimnames=dimnames(dup_numeric))
dup_metadata<-data.frame(grna_id=c('g1','g2','g1'),grna_target='geneX',vector_id='v')
simple_dup<-methods::new('dgTMatrix',i=c(0L,0L,1L),j=c(0L,0L,0L),x=c(2,3,7),Dim=c(2L,2L),Dimnames=list(c('g1','g2'),c('c2','c1')))
fixtures$simple_duplicate_numeric_and_metadata<-list(m=simple_dup,d=data.frame(grna_id=c('g2','g1','g1'),grna_target='target',vector_id=c('v2','v1','v1')))
fixtures$duplicate_numeric_T_rounding<-list(m=dup_numeric,d=dup_metadata)
fixtures$duplicate_logical_T<-list(m=dup_logical,d=dup_metadata)
for(label in c('large_plus_ones','fractional','cancellation_helper_only','overflow_helper_only')) {
 values<-switch(label,large_plus_ones=c(2^53,1,1),fractional=c(.1,.2,.3),cancellation_helper_only=c(1e16,1,-1e16),overflow_helper_only=c(1e308,1e308,0))
 dm<-matrix(values,3L,1L,dimnames=list(c('a','b','c'),'cell'))
 dd<-data.frame(grna_id=c('a','b','c'),grna_target='t',vector_id='v')
 fixtures[[paste0(label,'_dense')]]<-list(m=dm,d=dd)
 fixtures[[paste0(label,'_C')]]<-list(m=methods::as(Matrix::Matrix(dm,sparse=TRUE),'CsparseMatrix'),d=dd)
}
fixtures$factor_vector_ids_helper_only<-list(m=t,d=transform(d,vector_id=factor(vector_id,levels=c('vB','vA','vZ'))))
fixtures$factor_grna_ids_helper_only<-list(m=t,d=transform(d,grna_id=factor(grna_id)))
na_vectors<-d;na_vectors$vector_id[c(1L,4L)]<-NA_character_
fixtures$NA_vector_ids_helper_only<-list(m=t,d=na_vectors)
# Metadata deduplication is full-row distinct after dropping guide ID: different
# annotations for one vector remain separate metadata rows, while matrix sums
# still use every original metadata row.
metadata_extra<-d;metadata_extra$annotation[5L]<-'different-guide-annotation'
fixtures$metadata_distinct_extra_columns<-list(m=t,d=metadata_extra)

records<-lapply(names(fixtures),function(label) {
 f<-fixtures[[label]];stopifnot(!isS4(f$m)||isTRUE(methods::validObject(f$m,test=TRUE)))
 ref_input<-clone(f);alt_input<-clone(f);before_ref<-sha(ref_input);before_alt<-sha(alt_input)
 ref<-capture(original('collapse_grna_matrix')(ref_input$m,ref_input$d))
 alt<-capture(sparse_probe(alt_input$m,alt_input$d))
 metadata_out<-capture(original('collapse_grna_target_data_frame')(clone(f$d)))
 list(label=label,fixture_sha256=sha(f),input=snapshot(f$m),metadata=f$d,metadata_output=metadata_out$value,metadata_error=metadata_out$error,original=if(is.null(ref$error))snapshot(ref$value) else NULL,original_error=ref$error,original_warnings=ref$warnings,probe=if(is.null(alt$error))snapshot(alt$value) else NULL,probe_error=alt$error,probe_warnings=alt$warnings,exact_output_identical=identical(ref$value,alt$value),equal_at_1e_8=if(is.null(ref$error)&&is.null(alt$error))all.equal(ref$value,alt$value,tolerance=1e-8) else NULL,original_input_unchanged=identical(before_ref,sha(ref_input)),probe_input_unchanged=identical(before_alt,sha(alt_input)))
});names(records)<-names(fixtures)
saveRDS(fixtures,file.path(root,'fixtures.rds'),compress=FALSE)

# Small public imports verify the final sparse class and metadata/cell ordering.
response<-matrix(c(1,2,3,4,5,6,2,1,3),3L,3L,dimnames=list(c('rA','rB','rC'),colnames(m)))
public_aggregation<-lapply(c('numeric_C','numeric_R','numeric_T_unsorted_stored_zero','logical_T_unsorted_stored_false'),function(label) {
 f<-clone(fixtures[[label]]);before<-sha(f)
 z<-capture(sceptre::import_data(response_matrix=clone(response),grna_matrix=f$m,grna_target_data_frame=f$d,moi='low'))
 list(label=label,error=z$error,warnings=z$warnings,input_unchanged=identical(before,sha(f)),output=if(is.null(z$error))list(grna=snapshot(sceptre::get_grna_matrix(z$value)),response=snapshot(sceptre::get_response_matrix(z$value)),metadata=z$value@grna_target_data_frame,covariate_row_names=rownames(z$value@covariate_data_frame)) else NULL)
})
# Safe sorted ordinary Cell Ranger input: two directories in supplied batch order.
feature_lines<-c('rA\tMT-A\tGene Expression','rB\tB\tGene Expression','rC\tC\tGene Expression','gA\tA\tCRISPR Guide Capture','gB\tB\tCRISPR Guide Capture','gC\tC\tCRISPR Guide Capture')
dirs<-file.path(root,c('fixture-batch-second','fixture-batch-first'))
for(batch in seq_along(dirs)) {
 dir.create(dirs[batch],showWarnings=FALSE)
 writeLines(feature_lines,file.path(dirs[batch],'features.tsv'))
 writeLines(c('cell-z','cell-a'),file.path(dirs[batch],'barcodes.tsv'))
 vals<-as.numeric(seq_len(12L)+batch);vals[1L]<-0;vals[2L]<-2.5
 rows<-rep(seq_len(6L),2L);cols<-rep(seq_len(2L),each=6L)
 writeLines(c('%%MatrixMarket matrix coordinate real general','% bounded valid grouped fixture','6 2 12',paste(rows,cols,vals)),file.path(dirs[batch],'matrix.mtx'))
}
file_metadata<-data.frame(grna_id=c('gA','gB','gC'),grna_target=c('geneX','geneY','non-targeting'))
public_file<-capture(sceptre::import_data_from_cellranger(directories=dirs,moi='low',grna_target_data_frame=file_metadata,use_ondisc=FALSE))
file_record<-list(error=public_file$error,warnings=public_file$warnings,messages=public_file$messages,fixture_files=setNames(lapply(unlist(lapply(dirs,list.files,full.names=TRUE)),file_sha),unlist(lapply(dirs,list.files,full.names=TRUE))),output=if(is.null(public_file$error))list(response=snapshot(sceptre::get_response_matrix(public_file$value)),grna=snapshot(sceptre::get_grna_matrix(public_file$value)),metadata=public_file$value@grna_target_data_frame,covariate_row_names=rownames(public_file$value@covariate_data_frame),batch=as.character(public_file$value@covariate_data_frame$batch)) else NULL)

# Trace only the small control fixture to show where full-input as.matrix occurs.
trace_events<-list();colsum_events<-list()
record_as_matrix<-function(x) {
 trace_events[[length(trace_events)+1L]]<<-list(class=class(x),dim=dim(x))
}
record_colsums<-function(x) {
 colsum_events[[length(colsum_events)+1L]]<<-list(class=class(x),dim=dim(x))
}
lookup_env<-ns
while(!exists('as.matrix',envir=lookup_env,inherits=FALSE))lookup_env<-parent.env(lookup_env)
trace('as.matrix',where=lookup_env,tracer=quote(.GlobalEnv$record_as_matrix(x)),print=FALSE)
trace('colSums',where=asNamespace('Matrix'),tracer=quote(.GlobalEnv$record_colsums(x)),print=FALSE)
trace_input<-clone(fixtures$numeric_C)
trace_events<-list();colsum_events<-list();invisible(original('collapse_grna_matrix')(trace_input$m,trace_input$d));original_trace<-list(as_matrix=trace_events,colSums=colsum_events)
trace_events<-list();colsum_events<-list();invisible(sparse_probe(trace_input$m,trace_input$d));probe_trace<-list(as_matrix=trace_events,colSums=colsum_events)
untrace('as.matrix',where=lookup_env)
untrace('colSums',where=asNamespace('Matrix'))
sparse_colsum_methods<-setNames(lapply(c('dgCMatrix','dgRMatrix','dgTMatrix','lgCMatrix','lgRMatrix','lgTMatrix'),function(cls)deparse(body(methods::selectMethod('colSums',cls)))),c('dgCMatrix','dgRMatrix','dgTMatrix','lgCMatrix','lgRMatrix','lgTMatrix'))
# Stable intended scratch ordering, independent of any candidate import helper.
unsorted_rows<-data.frame(cell_idx=c(1L,0L,1L,0L),feature_idx=c(1L,1L,0L,0L),x=c(11,1,10,0),input_row=1:4)
ord<-order(unsorted_rows$cell_idx,unsorted_rows$feature_idx,unsorted_rows$input_row,method='radix')
ordered<-unsorted_rows[ord,,drop=FALSE]
intended_sorted<-raw_csc(ordered$feature_idx,c(0L,cumsum(tabulate(ordered$cell_idx+1L,nbins=2L))),ordered$x,c(2L,2L))
stopifnot(isTRUE(methods::validObject(intended_sorted,test=TRUE)))
rounding_sweep<-readRDS(file.path(root,'rounding-sweep.rds'))

report<-list(observed_at_utc=format(Sys.time(),tz='UTC',usetz=TRUE),reference_release_build=jsonlite::read_json('/private/tmp/sceptre-resource-reference-release-build.json',simplifyVector=FALSE),reference_release_build_sha256=file_sha('/private/tmp/sceptre-resource-reference-release-build.json'),status='READ-ONLY SMALL CHARACTERIZATION; NO PRODUCTION CHANGE OR RESOURCE QUALIFICATION',plan_revision='6a43fa5b2d32109a67e42138d3d845a6e5a2805a',original_source_revision='21f9ea098b69c444a884a49bc254de711968e3b5',reference_library=normalizePath(find.package('sceptre')),source_function_matches=source_match,source_manifest=jsonlite::read_json(file.path(root,'source-manifest.json'),simplifyVector=FALSE),runtime=list(R=R.version.string,platform=R.version$platform,session=capture.output(sessionInfo()),thread_env=as.list(Sys.getenv(c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'))),datatable_threads=data.table::getDTthreads(),dependencies=setNames(lapply(c('sceptre','Matrix','data.table','ondisc','dplyr','digest','jsonlite'),function(p)list(version=as.character(packageVersion(p)),path=find.package(p))),c('sceptre','Matrix','data.table','ondisc','dplyr','digest','jsonlite')),dlls=lapply(c('sceptre','Matrix','data_table','ondisc'),function(name){z<-getLoadedDLLs()[[name]];if(is.null(z))NULL else list(name=name,path=z[['path']],sha256=file_sha(z[['path']]))})),B1=list(raw_layout_validity_only=raw_records,grouped_unique_actual_create=snapshot(csc),grouped_input_unchanged=identical(valid_before,sha(valid_dt)),interleaved_intended_mapping=intended_interleaved,stable_order_intended=list(original_rows=unsorted_rows,stable_order=ord,ordered_rows=ordered,result=snapshot(intended_sorted)),public_sorted_cellranger=file_record,duplicate_policy_decision=list(status='EXPLICIT USER DECISION RELAYED BY PARENT DURING THIS CHARACTERIZATION',decision='Reject repeated file (cell,feature) coordinates with a clear import error; do not coalesce/sum. Existing valid dgT aggregation and repeated metadata multiplicity remain unchanged.',implementation='Not implemented here. No malformed CSC was sent to conversion, cbind, statistical or other downstream native kernels.')),B6=list(cases=records,public_aggregation=public_aggregation,full_input_as_matrix_trace=list(binding_environment=environmentName(lookup_env),original=original_trace,sparse_probe=probe_trace,limitation='Trace records this R-level full-input coercion and colSums argument representations; it is not allocation/RSS measurement.'),selected_sparse_colSums_methods=sparse_colsum_methods,rounding_sweep=list(seed=rounding_sweep$seed,cases=rounding_sweep$cases,exact_mismatches=length(rounding_sweep$mismatches),outside_1e_8=length(rounding_sweep$outside),artifact_sha256=file_sha(file.path(root,'rounding-sweep.rds'))),fixture_file=list(path=file.path(root,'fixtures.rds'),sha256=file_sha(file.path(root,'fixtures.rds'))),prototype=list(status='Characterization scratch only; not production/candidate qualification',body=deparse(sparse_probe),hash=sha(deparse(sparse_probe)),workspace='Sparse selected-group submatrix plus required dense vector-by-cell output; worst-case repeated metadata can expand sparse subgroup nnz. No whole-input as.matrix in prototype. No timing/RSS claims made.'),capability='No new backed aggregation. Original direct in-memory representation classes, output and metadata are characterized; full paired resource gates remain pending.'),summary=list(case_count=length(records),exact_equal=sum(vapply(records,function(x)isTRUE(x$exact_output_identical),logical(1))),original_inputs_unchanged=all(vapply(records,function(x)x$original_input_unchanged,logical(1))),probe_inputs_unchanged=all(vapply(records,function(x)x$probe_input_unchanged,logical(1)))))
jsonlite::write_json(report,file.path(root,'report.json'),pretty=TRUE,auto_unbox=TRUE,null='null',na='string',digits=17)
cat('B6 cases:',length(records),'exact matches:',report$summary$exact_equal,'\n')
for(r in records) if(!isTRUE(r$exact_output_identical))cat('MISMATCH:',r$label,'original error',r$original_error,'probe error',r$probe_error,'\n')
cat('Public file import error:',file_record$error,'\n')
cat('Public aggregation errors:',paste(vapply(public_aggregation,function(x)if(is.null(x$error))'none' else x$error,character(1)),collapse=' | '),'\n')
