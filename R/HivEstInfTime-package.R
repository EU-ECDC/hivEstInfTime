#' @name hivEstInfTime
#'
#' @title
#' HIV Infection Time Estimation
#'
#' @description
#' Estimate HIV infection time
#'
#' @author
#' Author: Nikos Pantazis \email{npantaz@@med.uoa.gr}\cr
#' Author: Magdalena Rosinska \email{mrosinska@@pzh.gov.pl}\cr
#' Author: Ard van Sighem \email{a.i.vansighem@amsterdamumc.nl}\cr
#' Creator: Daniel Lewandowski \email{daniel@@nextpagesoft.net}
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom stats coef formula model.matrix runif
#' @importFrom utils read.csv
#' @import data.table
#'
#' @useDynLib hivEstInfTime, .registration = TRUE
"_PACKAGE"
