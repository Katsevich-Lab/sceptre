status <- 'FAILURE'
sourcepkg <- commandArgs(TRUE)[1]

write_output <- function(key, value){
  cat(paste0(key, '=', value, '\n'), file = Sys.getenv("GITHUB_OUTPUT"), append = TRUE)
}

# BiocCheck warns about non-standard repos
options(BiocManager.check_repositories = FALSE)
if(Sys.getenv("UNIVERSE_NAME") == 'bioc-release'){
  Sys.setenv(R_BIOC_VERSION = as.character(BiocManager:::.version_bioc("release")))
} else {
  Sys.setenv(R_BIOC_VERSION = as.character(BiocManager:::.version_bioc("devel")))
}

try({
	library(BiocManager)
	library(BiocCheck)
	# check-bioc-help seems to require some sort of admin access? Not really clear to me.
	results <- BiocCheck(sourcepkg, 'no-check-R-ver' = TRUE, 'no-check-bioc-help' = TRUE)
	if(length(results$error)){
		status <- 'ERROR'
	} else if(length(results$warning)){
		status <- 'WARNING'
	} else if(length(results$note)){
		status <- 'NOTE'
	} else {
		status <- 'OK'
	}
	restxt <- jsonlite::as_gzjson_b64(as.list(results$getNum()), auto_unbox = TRUE)
	write_output('CHECKRESULTS', restxt)
})

write_output('CHECKSTATUS', status)
