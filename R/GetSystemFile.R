GetSystemFile <- function(
  ...,
  package = 'HivEstInfTime'
) {
  return(system.file(..., package = package))
}
