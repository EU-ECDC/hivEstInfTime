IsError <- function(
  x
) {
  return(inherits(x, 'try-error'))
}
