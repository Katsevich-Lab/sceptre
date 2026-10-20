args<-commandArgs(TRUE)
stopifnot(length(args)==2L)
base<-normalizePath(args[1L])
root<-normalizePath(args[2L])
e<-pkgload::load_all(root,compile=FALSE,helpers=FALSE,attach_testthat=FALSE,quiet=TRUE)$env
plain<-function(x){
 if(is.character(x))return(gsub('[[:space:]]+',' ',trimws(x)))
 if(is.list(x))return(lapply(x,plain))
 x
}
tags<-function(block){
 result<-list()
 for(t in block$tags){
  if(t$tag %in% c('usage','.formals','backref'))next
  if(t$tag=='param'){
   names<-unlist(strsplit(t$val$name,',',fixed=TRUE))
   for(name in names)result[[paste0('param:',name)]]<-plain(t$val$description)
  }else{
   key<-t$tag;count<-if(length(result))sum(startsWith(names(result),paste0(key,':'))) else 0L
   result[[paste0(key,':',count+1L)]]<-plain(t$val)
  }
 }
 result[sort(names(result))]
}
files<-setdiff(list.files(file.path(base,'R'),pattern='[.]R$'),'s4_count_store.R')
rows<-list()
for(file in files){
 b<-roxygen2::parse_file(file.path(base,'R',file),env=e)
 a<-roxygen2::parse_file(file.path(root,'R',file),env=e)
 b<-lapply(b,tags);a<-lapply(a,tags)
 rows[[file]]<-data.frame(file=paste0('R/',file),blocks=length(a),roxygen_semantics_identical=identical(a,b))
}
result<-do.call(rbind,rows);rownames(result)<-NULL
write.csv(result,file.path(dirname(base),'roxygen-equivalence.csv'),row.names=FALSE)
print(result)
stopifnot(all(result$roxygen_semantics_identical))
cat('PASS parsed roxygen tag/parameter semantics in ',nrow(result),' files; parameter groups expanded before comparison\n',sep='')
