GetSystemFile <- function(
  ...,
  package = 'hivEstInfTime'
) {
  return(system.file(..., package = package))
}
