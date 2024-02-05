# .onLoad <- function(libname, pkgname)
# {
#   library.dynam("batchCurve", pkgname, libname)
# }

.onAttach <- function(libname, pkgname) {

  t <- paste0(as.vector(utils::lsf.str(paste0("package:",pkgname))), collapse = "\n")

  msg <- paste0("_______________________________________\n  This is version ", utils::packageVersion(pkgname),
         " of ", pkgname,
         "\n____________________________________________\n",
         "Available functions:\n", t)
  packageStartupMessage(msg)
}
