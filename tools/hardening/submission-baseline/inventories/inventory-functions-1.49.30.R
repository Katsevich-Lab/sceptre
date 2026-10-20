ns <- asNamespace('BiocCheck')
base <- '/private/tmp/sceptre-submission-matched-check'
for (role in c('submitted','qualified')) {
 root <- file.path(base,role,'sceptre')
 inventory <- list()
 for (path in sort(list.files(file.path(root,'R'),pattern='\\.[Rr]$',full.names=TRUE))) {
  data <- getParseData(parse(path,keep.source=TRUE),includeText=TRUE)
  res <- get('getFunctionLengths',ns)(data)
  if(!length(res) || !length(names(res))) next
  for(i in seq_along(res)) inventory[[length(inventory)+1L]] <- data.frame(file=paste0('R/',basename(path)),name=names(res)[i],as.list(res[[i]]),stringsAsFactors=FALSE)
 }
 out<-do.call(rbind,inventory)
 write.csv(out,file.path(base,paste0(role,'-functions-1.49.30.csv')),row.names=FALSE)
 cat(role,'long functions',sum(out$length>50),'\n')
 print(out[out$length>50,c('file','name','length','startLine','endLine')],row.names=FALSE)
}
