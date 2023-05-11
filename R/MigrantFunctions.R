GetBaseCD4DesignMatrix <- function(
  data
) {
  dm <- model.matrix(
    formula(
      ~
      DTime * Gender +
        DTime * MigrantRegionOfOrigin +
        DTime * Transmission +
        DTime * Lspline(Age, knots = c(25, 35, 45)) +
        Lspline(Calendar, knots = c(16, 22))
    ),
    data = data
  )
  colMapping <- attr(dm, 'assign')

  return(list(
    dm = dm,
    colsDTime = which(colMapping == 1),
    colsGender = which(colMapping == 2),
    colsRegion = which(colMapping == 3),
    colsTrans = which(colMapping == 4),
    colsAge = which(colMapping == 5),
    colsCalendar = which(colMapping == 6),
    colsDTimeGender = which(colMapping == 7),
    colsDTimeRegion = which(colMapping == 8),
    colsDTimeTrans = which(colMapping == 9),
    colsDTimeAge = which(colMapping == 10)
  ))
}

GetBaseVLDesignMatrix <- function(
  data
) {
  dm <- model.matrix(
    formula(
      ~
        DTime * Gender +
          DTime * MigrantRegionOfOrigin +
          DTime * Transmission +
          DTime * Lspline(Age, knots = c(25, 35, 45)) +
          log(DTime) * Gender +
          log(DTime) * MigrantRegionOfOrigin +
          log(DTime) * Transmission +
          log(DTime) * Lspline(Age, knots = c(25, 35, 45)) +
          Lspline(Calendar, knots = c(16, 22))
    ),
    data = data
  )
  colMapping <- attr(dm, 'assign')

  return(list(
    dm = dm,
    colsDTime = which(colMapping == 1),
    colsGender = which(colMapping == 2),
    colsRegion = which(colMapping == 3),
    colsTrans = which(colMapping == 4),
    colsAge = which(colMapping == 5),
    colsLogDTime = which(colMapping == 6),
    colsCalendar = which(colMapping == 7),
    colsDTimeGender = which(colMapping == 8),
    colsDTimeRegion = which(colMapping == 9),
    colsDTimeTrans = which(colMapping == 10),
    colsDTimeAge = which(colMapping == 11),
    colsLogDTimeGender = which(colMapping == 12),
    colsLogDTimeRegion = which(colMapping == 13),
    colsLogDTimeTrans = which(colMapping == 14),
    colsLogDTimeAge = which(colMapping == 15)
  ))
}

GetBaseRandEffDesignMatrix <- function(
  data
) {
  dm <- model.matrix(
    formula(
      ~
        -1 +
          Consc + I(DTime * Consc) +
          Consr + I(DTime * Consr) + I(log(DTime + 1) * Consr)
    ),
    data = data
  )
  colMapping <- attr(dm, 'assign')

  return(list(
    dm = dm,
    colsConsc = which(colMapping == 1),
    colsDTimeConsc = which(colMapping == 2),
    colsConsr = which(colMapping == 3),
    colsDTimeConsr = which(colMapping == 4),
    colsLogDTimeConsr = which(colMapping == 5)
  ))
}

PostWAIDS <- function(
  w,
  x,
  dTime,
  betaAIDS,
  kappa
) {
  x[3] <- x[3] - w
  lambda <- exp(x %*% betaAIDS)[1, 1]

  p <- exp(log(kappa) + log(lambda) + (kappa - 1) * log(w + dTime) - lambda * (w + dTime)^kappa)
  return(p)
}

# Vectorize the functions as the "integrate" function works with vectorized functions
VPostWAIDS <- Vectorize(PostWAIDS, vectorize.args = c('w'))
