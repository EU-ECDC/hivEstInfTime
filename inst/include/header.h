#ifndef _hivPlatform_header_
#define _hivPlatform_header_

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#define NDEBUG
#define ARMA_NO_DEBUG
#define BOOST_DISABLE_ASSERTS
#define EIGEN_NO_DEBUG
#define EIGEN_NO_STATIC_ASSERT

#include <RcppArmadillo.h>
// RcppArmadillo undefines NDEBUG. Define it again.
#define NDEBUG
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppNumerical.h>

#include "Lspline.hpp"
#include "GetLogMVNPdf.hpp"
#include "PostW.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#endif // _hivPlatform_header_
