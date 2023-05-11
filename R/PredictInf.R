#' PredictInf
#'
#' Compute probability of infection prior to migration.
#'
#' @param input input
#' @param params params
#' @param sampleSize sampleSize
#'
#' @return data.table
#'
#' @examples
#' \dontrun{
#' PredictInf(input, params = GetMigrantParams(), sampleSize = 50)
#' }
#'
#' @export
PredictInf <- function( # nolint
  input,
  params = GetMigrantParams(),
  sampleSize = 50
) {
  outputAIDS <- data.table::copy(input$Data$AIDS)
  outputCD4VL <- data.table::copy(input$Data$CD4VL)
  sampleColNames <- paste0('SCtoDiag', seq_len(sampleSize))

  if (is.null(outputAIDS)) {
    outputAIDS <- data.table::data.table(
      Imputation = integer(),
      RecordId = character(),
      UniqueId = integer(),
      Ord = integer()
    )
  }

  if (is.null(outputCD4VL)) {
    outputCD4VL <- data.table::data.table(
      Imputation = integer(),
      RecordId = character(),
      UniqueId = integer(),
      Ord = integer()
    )
  }

  # AIDS -------------------------------------------------------------------------------------------
  countAIDS <- outputAIDS[, .N]
  countAIDSNChar <- nchar(as.character(countAIDS))
  outputAIDS[, ':='(
    ProbPre = NA_real_,
    EstSCtoDiag = NA_real_
  )]
  outputAIDS[, (sampleColNames) := NA_real_]

  xAIDS <- cbind(
    1,
    as.integer(outputAIDS$Gender == 'M'),
    outputAIDS$Age
  )
  startTime <- Sys.time()
  lastTime <- startTime
  PrintH1('Processing AIDS data')
  PrintAlert('Start time: {format(startTime)}')
  for (i in seq_len(countAIDS)) {
    if (i %% 1000 == 0) {
      currentTime <- Sys.time()
      percComplete <- stringi::stri_pad_left(sprintf('%0.2f%%', i / countAIDS * 100), width = 8)
      iterComplete <- stringi::stri_pad_left(sprintf('%d/%d', i, countAIDS), width = countAIDSNChar * 2 + 1) # nolint
      elapsedTime <- prettyunits::pretty_dt(currentTime - startTime)
      PrintAlert('{percComplete} ({iterComplete}) - {.timestamp {elapsedTime}}', type = 'success')
      lastTime <- currentTime
    }

    knownPrePost <- outputAIDS$KnownPrePost[i]
    if (knownPrePost != 'Unknown') {
      outputAIDS[i, ':='(
        ProbPre = data.table::fifelse(knownPrePost == 'Pre', 1.0, 0.0),
        EstSCtoDiag = NA_real_
      )]
      outputAIDS[i, (sampleColNames) := 0]
      next
    }

    intFit1 <- try(stats::integrate(
      VPostWAIDS,
      lower = 0,
      upper = outputAIDS$Mig[i],
      x = xAIDS[i, ],
      dTime = outputAIDS$DTime[i],
      betaAIDS = params$betaAIDS,
      kappa = params$kappa
    ), silent = TRUE)

    intFit2 <- try(stats::integrate(
      VPostWAIDS,
      lower = outputAIDS$Mig[i],
      upper = outputAIDS$U[i],
      x = xAIDS[i, ],
      dTime = outputAIDS$DTime[i],
      betaAIDS = params$betaAIDS,
      kappa = params$kappa
    ), silent = TRUE)

    if (IsError(intFit1) || IsError(intFit2) || intFit1$message != 'OK' || intFit2$message != 'OK') { # nolint
      next
    } else {
      outputAIDS[i, ProbPre := intFit2$value / (intFit1$value + intFit2$value)]

      modeFit <- try(stats::optim(
        outputAIDS$Mig[i],
        PostWAIDS,
        method = 'Brent',
        control = list(fnscale = -1),
        lower = 0,
        upper = outputAIDS$U[i],
        x = xAIDS[i, ],
        dTime = outputAIDS$DTime[i],
        betaAIDS = params$betaAIDS,
        kappa = params$kappa
      ), silent = TRUE)

      if (IsError(modeFit) || modeFit$convergence != 0) {
        next
      } else {
        outputAIDS[i, ModeSCtoDiag := modeFit$par]
        outputAIDS[i, (sampleColNames) := as.list(PerformRejectionSampling(
          n = sampleSize,
          density = VPostWAIDS,
          mode = outputAIDS$ModeSCtoDiag[i],
          lower = 0,
          upper = outputAIDS$U[i],
          x = xAIDS[i, ],
          dTime = outputAIDS$DTime[i],
          betaAIDS = params$betaAIDS,
          kappa = params$kappa
        ))]
      }
    }
  }
  endTime <- Sys.time()
  if (countAIDS > 0) {
    percComplete <- stringi::stri_pad_left(sprintf('%0.2f%%', i / countAIDS * 100), width = 8)
    iterComplete <- stringi::stri_pad_left(sprintf('%d/%d', i, countAIDS), width = countAIDSNChar * 2 + 1) # nolint
    elapsedTime <- prettyunits::pretty_dt(endTime - startTime)
    PrintAlert('{percComplete} ({iterComplete}) - {.timestamp {elapsedTime}}', type = 'success')
  } else {
    PrintAlert('No AIDS data to be processed')
  }
  PrintAlert('End time: {format(endTime)}')

  # CD4VL ------------------------------------------------------------------------------------------
  countCD4VL <- outputCD4VL[, data.table::uniqueN(UniqueId)]
  countCD4VLNChar <- nchar(as.character(countCD4VL))
  outputCD4VL[, ':='(
    ProbPre = NA_real_,
    EstSCtoDiag = NA_real_
  )]
  outputCD4VL[, (sampleColNames) := NA_real_]

  i <- 0
  startTime <- Sys.time()
  lastTime <- startTime
  PrintH1('Processing CD4VL data')
  PrintAlert('Start time: {format(startTime)}')
  for (uniqueId in outputCD4VL[, unique(UniqueId)]) {
    i <- i + 1
    if (i %% 1000 == 0) {
      currentTime <- Sys.time()
      percComplete <- stringi::stri_pad_left(sprintf('%0.2f%%', i / countCD4VL * 100), width = 8)
      iterComplete <- stringi::stri_pad_left(sprintf('%d/%d', i, countCD4VL), width = countCD4VLNChar * 2 + 1) # nolint
      elapsedTime <- prettyunits::pretty_dt(currentTime - startTime)
      PrintAlert('{percComplete} ({iterComplete}) - {.timestamp {elapsedTime}}', type = 'success')
      lastTime <- currentTime
    }

    dt <- outputCD4VL[UniqueId == uniqueId]

    # Predetermined status results in a fixed probability
    knownPrePost <- dt[Ord == 1, KnownPrePost]
    if (knownPrePost %chin% c('Pre', 'Post')) {
      outputCD4VL[UniqueId == uniqueId, ':='(
        ProbPre = data.table::fifelse(knownPrePost == 'Pre', 1.0, 0.0),
        EstSCtoDiag = 0
      )]
      outputCD4VL[UniqueId == uniqueId, (sampleColNames) := 0]
      next
    }

    migTime <- dt[Ord == 1, Mig]
    upTime <- dt[Ord == 1, U]
    y <- dt$YVar
    xAIDS <- as.matrix(dt[Ord == 1, .(1, as.integer(Gender == 'M'), Age)])
    maxDTime <- dt[, max(DTime)]
    betaAIDS <- matrix(params$betaAIDS, ncol = 1)
    kappa <- params$kappa
    bFE <- matrix(params$bFE, ncol = 1)
    varCovRE <- params$varCovRE
    formulaeData <- dt[, .(Gender, MigrantRegionOfOrigin, Transmission, Age, DTime, Calendar, Consc, Consr)] # nolint
    fxCD4Data <- formulaeData[Consc == 1]
    baseCD4DM <- GetBaseCD4DesignMatrix(fxCD4Data)
    fxVLData <- formulaeData[Consr == 1]
    baseVLDM <- GetBaseVLDesignMatrix(fxVLData)
    fzData <- dt[, .(Consc, CobsTime, Consr, RobsTime, RLogObsTime2, DTime)]
    baseRandEffDM <- GetBaseRandEffDesignMatrix(fzData)

    sigma2 <- params$sigma2
    errM <- dt$Consc * sigma2[1] + dt$Consr * sigma2[2]
    err <- diag(errM, nrow = length(errM))

    intFit1 <- try(IntegratePostW(
      lower = 0,
      upper = migTime,
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

    intFit2 <- try(IntegratePostW(
      lower = migTime,
      upper = upTime,
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

    if (
      IsError(intFit1) || IsError(intFit2) ||
        intFit1$errorCode != 0 || intFit2$errorCode != 0 ||
        is.na(intFit1$value) || is.na(intFit2$value)
    ) {
      next
    } else {
      outputCD4VL[UniqueId == uniqueId, ProbPre := intFit2$value / (intFit1$value + intFit2$value)]

      modeFit <- try(stats::optim(
        migTime,
        PostW,
        method = 'Brent',
        control = list(fnscale = -1),
        lower = 0,
        upper = upTime,
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

      if (IsError(modeFit) || modeFit$convergence != 0) {
        next
      } else {
        outputCD4VL[UniqueId == uniqueId, ModeSCtoDiag := modeFit$par]
        outputCD4VL[UniqueId == uniqueId, (sampleColNames) := as.list(PerformRejectionSampling(
          n = sampleSize,
          density = VPostW,
          mode = modeFit$par,
          lower = 0,
          upper = upTime,
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
        ))]
      }
    }
  }
  endTime <- Sys.time()
  if (countCD4VL > 0) {
      percComplete <- stringi::stri_pad_left(sprintf('%0.2f%%', i / countCD4VL * 100), width = 8)
      iterComplete <- stringi::stri_pad_left(sprintf('%d/%d', i, countCD4VL), width = countCD4VLNChar * 2 + 1) # nolint
      elapsedTime <- prettyunits::pretty_dt(endTime - startTime)
      PrintAlert('{percComplete} ({iterComplete}) - {.timestamp {elapsedTime}}', type = 'success')
  } else {
    PrintAlert('No CD4VL data to be processed')
  }
  PrintAlert('End time: {format(endTime)}')

  outputColNames <- union(c('UniqueId', 'ProbPre', 'Mig', 'ModeSCtoDiag'), sampleColNames)
  table <- list()
  if (nrow(outputAIDS) > 0) {
    table[['AIDS']] <- outputAIDS[, ..outputColNames]
  }
  if (nrow(outputCD4VL) > 0) {
    table[['CD4VL']] <- outputCD4VL[Ord == 1, ..outputColNames]
  }
  table <- rbindlist(table)

  if (nrow(table) == 0) {
    table <- data.table(
      UniqueId = integer(),
      ProbPre = numeric(),
      EstSCtoDiag = numeric()
    )
    table[, (sampleColNames) := NA_real_]
  }

  return(table)
}
