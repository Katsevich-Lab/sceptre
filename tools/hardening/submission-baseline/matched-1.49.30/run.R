options(width=100)
stopifnot(as.character(packageVersion('BiocCheck'))=='1.49.30')
p<-'/private/tmp/sceptre-submission-checker-1.49.30/baseline'
writeLines(c(R.version.string,capture.output(sessionInfo())),file.path(p,'session.txt'))
result<-BiocCheck::BiocCheck(file.path(p,'sceptre_0.99.0.tar.gz'),debug=TRUE,callr=FALSE,`no-check-R-ver`=TRUE,`no-check-bioc-help`=TRUE,`quit-with-status`=FALSE)
saveRDS(result,file.path(p,'result.rds')); result$toJSON(file.path(p,'result.json'))
str(result,max.level=2)
