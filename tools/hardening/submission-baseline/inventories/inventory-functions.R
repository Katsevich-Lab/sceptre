ns <- asNamespace('BiocCheck')
base <- '/private/tmp/sceptre-submission-matched-check'
for (role in c('submitted','qualified')) {
 root <- file.path(base,role,'sceptre')
 inventory <- list()
 for (path in sort(list.files(file.path(root,'R'),pattern='\\.[Rr]$',full.names=TRUE))) {
  data <- getParseData(parse(path,keep.source=TRUE),includeText=TRUE)
  res <- get('getFunctionLengths',ns)(data)
  if(!length(res)) next
  for(n in names(res)) inventory[[length(inventory)+1L]] <- data.frame(file=paste0('R/',basename(path)),name=n,as.list(res[[n]]),stringsAsFactors=FALSE)
 }
 out<-do.call(rbind,inventory)
 write.csv(out,file.path(base,paste0(role,'-functions.csv')),row.names=FALSE)
 cat(role,'long functions',sum(out$length>50),'\n')
 print(out[out$length>50,c('file','name','length','startLine','endLine')],row.names=FALSE)
}
