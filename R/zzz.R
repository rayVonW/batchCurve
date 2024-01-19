# .onLoad <- function(libname, pkgname)
# {
#   library.dynam("gamma", pkgname, libname)
# }

.onAttach <- function(libname, pkgname) {

  t <- paste0(as.vector(lsf.str(paste0("package:",pkgname))), collapse = "\n")

  msg <- paste0("_______________________________\n  This is version ", packageVersion(pkgname),
         " of ", pkgname,
         "\n_________________________________\n",
         "Available functions:\n", t)
  packageStartupMessage(msg)
}
