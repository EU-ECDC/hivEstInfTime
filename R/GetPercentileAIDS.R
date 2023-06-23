GetPercentileAIDS <- function(
  p = 0.5,
  v = NULL,
  upTime = NULL,
  lower = 0,
  x,
  dTime,
  betaAIDS,
  kappa
) {
  OptimFunc <- function(
    u,
    p,
    v,
    l,
    x,
    dTime,
    betaAIDS,
    kappa
    ) {
    res <- stats::integrate(
      VPostWAIDS,
      lower = l,
      upper = u,
      x = x,
      dTime = dTime,
      betaAIDS = betaAIDS,
      kappa = kappa
    )$value

    return(res / v - p)
  }

  fit <- try(stats::uniroot(
    OptimFunc,
    interval = c(0, upTime),
    p = p,
    v = v,
    l = lower,
    x = x,
    dTime = dTime,
    betaAIDS = betaAIDS,
    kappa = kappa
  ), silent = TRUE)

  percentile <- ifelse(
    IsError(fit) || is.na(fit$root),
    NA_real_,
    fit$root
  )

  return(percentile)
}
