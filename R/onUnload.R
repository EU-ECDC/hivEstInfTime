.onUnload <- function(libPath, pkgname) {

  library.dynam.unload('hivEstInfTime', libPath)

  invisible(NULL)
}
