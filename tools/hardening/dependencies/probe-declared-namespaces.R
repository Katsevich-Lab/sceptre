# Run with R_DEFAULT_PACKAGES=NULL: no attached methods package may rescue a
# missing namespace import. This qualifies the final methods::slotNames fix.
stopifnot(!"package:methods" %in% search())
namespace <- loadNamespace("sceptre")
empty <- methods::new("AnalysisSpec")
guard <- get("analysis_spec_is_empty", envir = namespace)
stopifnot(guard(empty))
nonempty <- empty
nonempty@multiple_testing_method <- "BH"
stopifnot(!guard(nonempty), !"package:methods" %in% search())
cat(R.version.string, "\n")
cat("Empty/nonempty AnalysisSpec guard passes without attached methods.\n")
