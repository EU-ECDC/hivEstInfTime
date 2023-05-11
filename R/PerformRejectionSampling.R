PerformRejectionSampling <- function(
  n,
  density,
  mode,
  lower,
  upper,
  verbose = FALSE,
  ...
) {
  out <- rep(0, n)
  k <- density(mode, ...)

  i <- 0
  acceptedCount <- 0
  while (acceptedCount < n) {
    i <- i + 1
    u <- runif(n - acceptedCount, lower, upper)
    v <- runif(n - acceptedCount)
    s <- density(u, ...)
    a <- v < s / k
    g <- sum(a)
    out[acceptedCount + seq_len(g)] <- u[a]
    acceptedCount <- acceptedCount + g
  }

  if (verbose) {
    PrintAlert('Iterations required: {.val {i}}')
    PrintAlert('Acceptance rate: {.val {n / i}}')
  }

  return(out)
}
