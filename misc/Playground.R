params <- GetMigrantParams()
names(params)
bFE <- params$bFE
bFE.CD4 <- params$bFECD4
bFEVL <- params$bFEVL
varCovRE <- params$varCovRE
VarCovRE.CD4 <- params$varCovRECD4
varCovREVL <- params$varCovREVL
sigma2 <- params$sigma2
sigma2.CD4 <- params$sigma2CD4
sigma2VL <- params$sigma2VL
beta.aids <- params$betaAIDS
kappa <- params$kappa

set.seed(10)
baseCD4VL$ord <- unlist(tapply(baseCD4VL$id, baseCD4VL$id, function(x) 1:length(x)))
# Get the design matrices
x <- split(baseCD4VL[, c("gender", "region_group", "mode", "ageDiag", "dtime", "calendar", "consc", "consr")], baseCD4VL$id)
z <- lapply(split(
  baseCD4VL[, c("consc", "Cobstime", "consr", "Robstime", "RlogObstime2", "dtime")],
  baseCD4VL$id
), function(x) as.matrix(x))
y <- split(baseCD4VL$yvar, baseCD4VL$id)
u <- split(baseCD4VL$u, baseCD4VL$id)
ids <- unique(baseCD4VL$id)
only <- split(baseCD4VL$only, baseCD4VL$id)
mig <- split(baseCD4VL$mig, baseCD4VL$id)
known <- split(baseCD4VL$KnownPrePost, baseCD4VL$id)

# One row per subject
baseCD4VL.id <- baseCD4VL[baseCD4VL$ord == 1, ]

xaids <- cbind(1, baseCD4VL.id$gender == "Male", baseCD4VL.id$ageDiag)
maxdtime <- tapply(baseCD4VL$dtime, baseCD4VL$id, max)
G <- length(unique(ids))
ind <- list()
for (i in 1:G) {
  ind[[i]] <- which(baseCD4VL$id == ids[i])
}

########################
### Posterior median ###
########################
baseCD4VL$estsctodiagMedian <- NA
baseCD4VL$estsctodiag2.5 <- NA
baseCD4VL$estsctodiag97.5 <- NA
baseCD4VL$estsctodiag30 <- NA
baseCD4VL$estsctodiag40 <- NA
baseCD4VL$estsctodiag60 <- NA

med_postW <- function(M) {
  integrate(VpostW, lower = 0, upper = M)$value / res1$value - 0.5
}

med_postWcd4 <- function(M) {
  integrate(VpostWcd4, lower = 0, upper = M)$value / res1$value - 0.5
}
med_postWcd4_2 <- function(M) {
  HivEstInfTime:::IntegratePostW(
    lower = 0,
    upper = M,
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
  )$value / res1$value - 0.5
}

med_postWvl <- function(M) {
  integrate(VpostWvl, lower = 0, upper = M)$value / res1$value - 0.5
}

#
ff_postW <- function(M, a = 0.5) {
  integrate(VpostW, lower = 0, upper = M)$value / res1$value - a
}

ff_postWcd4 <- function(M, a = 0.5) {
  integrate(VpostWcd4, lower = 0, upper = M)$value / res1$value - a
}

ff_postWvl <- function(M, a = 0.5) {
  integrate(VpostWvl, lower = 0, upper = M)$value / res1$value - a
}

i <- 1L
for (i in 1:G) {
  uptime <- u[[i]][1]

  if (only[[i]][1] == "Both") {
    res1 <- try(integrate(VpostW, lower = 0, upper = uptime), silent = TRUE)
    res2 <- try(uniroot(med_postW, interval = c(0, uptime)), silent = TRUE)

    if ('try-error' %in% class(res1) | 'try-error' %in% class(res2)) {
      next
    } else {
      Res <- res2$root
    }

    res3 <- try(uniroot(ff_postW, interval = c(0, uptime), a = 0.025),
      silent =
        TRUE
    )
    res4 <- try(uniroot(ff_postW, interval = c(0, uptime), a = 0.975),
      silent =
        TRUE
    )

    res30 <- try(uniroot(ff_postW, interval = c(0, uptime), a = 0.30),
      silent =
        TRUE
    )
    res40 <- try(uniroot(ff_postW, interval = c(0, uptime), a = 0.40),
      silent =
        TRUE
    )
    res60 <- try(uniroot(ff_postW, interval = c(0, uptime), a = 0.60),
      silent =
        TRUE
    )

    baseCD4VL$estsctodiagMedian[ind[[i]]] <- Res
    baseCD4VL$estsctodiag2.5[ind[[i]]] <- res3$root
    baseCD4VL$estsctodiag97.5[ind[[i]]] <- res4$root

    baseCD4VL$estsctodiag30[ind[[i]]] <- res30$root
    baseCD4VL$estsctodiag40[ind[[i]]] <- res40$root
    baseCD4VL$estsctodiag60[ind[[i]]] <- res60$root
  }

  if (only[[i]][1] == "CD4 only") {
    res1 <- try(integrate(VpostWcd4, lower = 0, upper = uptime), silent = TRUE)
    res1Mine <- try(HivEstInfTime:::IntegratePostW(
      lower = 0,
      upper = uptime,
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
    res1Mine <- intFit1$value + intFit2$value

    res2 <- try(stats::uniroot(med_postWcd4, interval = c(0, uptime)), silent = TRUE)
    res2Mine <- try(stats::uniroot(med_postWcd4_2, interval = c(0, upTime)), silent = TRUE)

    if ('try-error' %in% class(res1) | 'try-error' %in% class(res2)) {
      next
    } else {
      Res <- res2$root
    }

    res3 <- try(uniroot(ff_postWcd4, interval = c(0, uptime), a = 0.025),
      silent = TRUE
    )
    res4 <- try(uniroot(ff_postWcd4, interval = c(0, uptime), a = 0.975),
      silent = TRUE
    )
    res30 <- try(uniroot(ff_postWcd4, interval = c(0, uptime), a = 0.30),
      silent =
        TRUE
    )
    res40 <- try(uniroot(ff_postWcd4, interval = c(0, uptime), a = 0.40),
      silent =
        TRUE
    )
    res60 <- try(uniroot(ff_postWcd4, interval = c(0, uptime), a = 0.60),
      silent =
        TRUE
    )


    baseCD4VL$estsctodiagMedian[ind[[i]]] <- Res
    baseCD4VL$estsctodiag2.5[ind[[i]]] <- res3$root
    baseCD4VL$estsctodiag97.5[ind[[i]]] <- res4$root

    baseCD4VL$estsctodiag30[ind[[i]]] <- res30$root
    baseCD4VL$estsctodiag40[ind[[i]]] <- res40$root
    baseCD4VL$estsctodiag60[ind[[i]]] <- res60$root
  }

  if (only[[i]][1] == "VL only") {
    res1 <- try(integrate(VpostWvl, lower = 0, upper = uptime), silent = TRUE)
    res2 <- try(uniroot(med_postWvl, interval = c(0, uptime)), silent = TRUE)

    if ('try-error' %in% class(res1) | 'try-error' %in% class(res2)) {
      next
    } else {
      Res <- res2$root
    }

    res3 <- try(uniroot(ff_postWvl, interval = c(0, uptime), a = 0.025),
      silent =
        TRUE
    )
    res4 <- try(uniroot(ff_postWvl, interval = c(0, uptime), a = 0.975),
      silent =
        TRUE
    )

    res30 <- try(uniroot(ff_postWvl, interval = c(0, uptime), a = 0.30),
      silent =
        TRUE
    )
    res40 <- try(uniroot(ff_postWvl, interval = c(0, uptime), a = 0.40),
      silent =
        TRUE
    )
    res60 <- try(uniroot(ff_postWvl, interval = c(0, uptime), a = 0.60),
      silent =
        TRUE
    )


    baseCD4VL$estsctodiagMedian[ind[[i]]] <- Res
    baseCD4VL$estsctodiag2.5[ind[[i]]] <- res3$root
    baseCD4VL$estsctodiag97.5[ind[[i]]] <- res4$root

    baseCD4VL$estsctodiag30[ind[[i]]] <- res30$root
    baseCD4VL$estsctodiag40[ind[[i]]] <- res40$root
    baseCD4VL$estsctodiag60[ind[[i]]] <- res60$root
  }
  print(i)
}
