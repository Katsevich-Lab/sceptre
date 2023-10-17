.onLoad <- function(libname, pkgname) {
  cat(paste0(crayon::blue("Welcome to sceptre.\nSubmit issues on the sceptre website: "),
             crayon::red("github.com/Katsevich-Lab/sceptre"),
             crayon::blue("\nRead the sceptre manual: "),
             crayon::red("timothy-barry.github.io/sceptre-book/")))
}
