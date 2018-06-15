.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("Package 'DatabionicSwarm' version ", packageVersion("DatabionicSwarm"),".\n","Type 'citation('DatabionicSwarm')' for citing this R package in publications.\n"))
}
.onLoad <- function(libname, pkgname) {
  library.dynam("DatabionicSwarm", pkgname, libname)
}