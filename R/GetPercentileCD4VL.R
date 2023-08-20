GetPercentileCD4VL <- function(
  p = 0.5,
  v = NULL,
  upTime = NULL,
  lower = 0,
  y = y,
  xAIDS = xAIDS,
  maxDTime = maxDTime,
  betaAIDS = betaAIDS,
  kappa = kappa,
  bFE = bFE,
  varCovRE = varCovRE,
  baseCD4DM = baseCD4DM,
  fxCD4Data = fxCD4Data,
  baseVLDM = baseVLDM,
  fxVLData = fxVLData,
  baseRandEffDM = baseRandEffDM,
  fzData = fzData,
  err = err
) {
  OptimFunc <- function(
    p,
    v,
    upTime,
    lower = 0,
    y = y,
    xAIDS = xAIDS,
    maxDTime = maxDTime,
    betaAIDS = betaAIDS,
    kappa = kappa,
    bFE = bFE,
    varCovRE = varCovRE,
    baseCD4DM = baseCD4DM,
    fxCD4Data = fxCD4Data,
    baseVLDM = baseVLDM,
    fxVLData = fxVLData,
    baseRandEffDM = baseRandEffDM,
    fzData = fzData,
    err = err
    ) {
    res <- IntegratePostW(
      upper = upTime,
      lower = 0,
      y = y,
      xAIDS = xAIDS,
      maxDTime = maxDTime,
      betaAIDS = betaAIDS,
      kappa = kappa,
      bFE = bFE,
      varCovRE = varCovRE,
      baseCD4DM = baseCD4DM,
      fxCD4Data = fxCD4Data,
      baseVLDM = baseVLDM,
      fxVLData = fxVLData,
      baseRandEffDM = baseRandEffDM,
      fzData = fzData,
      err = err
    )$value

    return(res / v - p)
  }

  fit <- try(stats::uniroot(
    OptimFunc,
    interval = c(0, upTime),
    p = p,
    v = v,
    lower = lower,
    y = y,
    xAIDS = xAIDS,
    maxDTime = maxDTime,
    betaAIDS = betaAIDS,
    kappa = kappa,
    bFE = bFE,
    varCovRE = varCovRE,
    baseCD4DM = baseCD4DM,
    fxCD4Data = fxCD4Data,
    baseVLDM = baseVLDM,
    fxVLData = fxVLData,
    baseRandEffDM = baseRandEffDM,
    fzData = fzData,
    err = err
  ), silent = TRUE)

  percentile <- ifelse(
    IsError(fit) || is.na(fit$root),
    NA_real_,
    fit$root
  )

  return(percentile)
}
