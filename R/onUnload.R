.onUnload <- function(libPath, pkgname) {

  library.dynam.unload('HivEstInfTime', libPath)

  invisible(NULL)
}
