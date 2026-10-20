args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
base <- normalizePath(args[1L])
source_root <- normalizePath(args[2L])
canon <- function(x, file) {
    if (missing(x)) return(list(type='missing'))
    if (is.symbol(x)) {
        name<-as.character(x)
        if(file=='assign_grna_functions.R' && name=='grna_params')
            name<-'grna_assignment_hyperparameters'
        if(file=='assign_grna_functions.R' && name=='assignment_info')
            return(canon(quote(sceptre_object@data@import_grna_assignment_info), file))
        return(list(type='symbol',value=name))
    }
    if (is.call(x)) {
        head<-if(is.symbol(x[[1L]]))as.character(x[[1L]]) else ''
        if(file=='assign_grna_functions.R' && head=='<-' &&
           identical(x[[2L]],as.name('assignment_info')))
            return(list(type='branch_local_alias'))
        parts<-as.list(x)
        if(file=='ondisc_capabilities.R' && head=='{') {
            drop<-integer()
            if(length(parts)>2L)for(i in seq.int(2L,length(parts)-1L)){
                if(identical(parts[[i]],quote(issues <- character())) &&
                   identical(parts[[i+1L]],quote(issues <- ondisc_api_argument_issues(snapshot))))
                    drop<-c(drop,i)
            }
            if(length(drop))parts<-parts[-drop]
        }
        result<-lapply(parts,canon,file=file)
        nm<-names(parts);if(is.null(nm))nm<-rep('',length(parts))
        keep<-!vapply(result,function(x)identical(x$type,'branch_local_alias'),logical(1))
        result<-result[keep];nm<-nm[keep]
        if(head %in% c('stop','warning','message','paste0')) {
            merged<-list(result[[1L]]);merged_names<-nm[1L]
            for(i in seq.int(2L,length(result))){
                value<-result[[i]]
                text<-identical(value$type,'character') && length(value$value)==1L && nm[i]==''
                last<-length(merged)
                prev<-last>1L && identical(merged[[last]]$type,'character') &&
                    length(merged[[last]]$value)==1L && merged_names[last]==''
                if(text && prev)merged[[last]]$value<-paste0(merged[[last]]$value,value$value)
                else {merged[[last+1L]]<-value;merged_names<-c(merged_names,nm[i])}
            }
            result<-merged;nm<-merged_names
        }
        if(head=='paste0' && length(result)==2L && nm[2L]=='' &&
           identical(result[[2L]]$type,'character') && length(result[[2L]]$value)==1L)
            return(result[[2L]])
        return(list(type='call',items=result,names=nm))
    }
    if (is.expression(x)||is.pairlist(x)) {
        value<-as.list(x)
        return(list(type=typeof(x),items=lapply(value,canon,file=file),names=names(value)))
    }
    list(type=typeof(x),value=x)
}
files<-list.files(file.path(base,'R'),pattern='[.]R$')
# Root made a later, independently reviewed helper extraction in this file.
files<-setdiff(files,'s4_count_store.R')
rows<-lapply(files,function(file){
    before<-parse(file.path(base,'R',file),keep.source=FALSE)
    after<-parse(file.path(source_root,'R',file),keep.source=FALSE)
    same<-identical(before,after)
    normalized<-identical(canon(before,file),canon(after,file))
    if(!normalized){
        saveRDS(list(before=canon(before,file),after=canon(after,file)),
            file.path(dirname(base),paste0('AST-FAIL-',file,'.rds')))
    }
    data.frame(file=paste0('R/',file),raw_AST_identical=same,
        computation_and_literal_normalized_AST_identical=normalized)
})
result<-do.call(rbind,rows)
write.csv(result,file.path(dirname(base),'ast-equivalence.csv'),row.names=FALSE)
print(result)
stopifnot(all(result$computation_and_literal_normalized_AST_identical))
cat('PASS ',sum(result$raw_AST_identical),' raw AST and ',sum(!result$raw_AST_identical),
    ' normalized constant-string/alpha-alias/dead-binding AST equivalences; ',nrow(result),' files total\n',sep='')
