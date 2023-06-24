#ifndef _HivEstInfTime_UpdateCD4DesignMatrix_
#define _HivEstInfTime_UpdateCD4DesignMatrix_

#include "const.h"

namespace HivEstInfTime {

arma::dmat UpdateCD4DesignMatrix(
  const Rcpp::List& baseDM,
  const Rcpp::DataFrame& data,
  const double& w
) {

  arma::dmat dm = Rcpp::as<arma::dmat>(baseDM["dm"]);

  if (dm.n_rows == 0) {
    return dm;
  }

  const arma::uvec& colsDTime = Rcpp::as<arma::uvec>(baseDM["colsDTime"]);
  const arma::uvec& colsAge = Rcpp::as<arma::uvec>(baseDM["colsAge"]);
  const arma::uvec& colsCalendar = Rcpp::as<arma::uvec>(baseDM["colsCalendar"]);
  const arma::uvec& colsDTimeGender = Rcpp::as<arma::uvec>(baseDM["colsDTimeGender"]);
  const arma::uvec& colsDTimeRegion = Rcpp::as<arma::uvec>(baseDM["colsDTimeRegion"]);
  const arma::uvec& colsDTimeTrans = Rcpp::as<arma::uvec>(baseDM["colsDTimeTrans"]);
  const arma::uvec& colsDTimeAge = Rcpp::as<arma::uvec>(baseDM["colsDTimeAge"]);
  const arma::uvec& colsGender = Rcpp::as<arma::uvec>(baseDM["colsGender"]);
  const arma::uvec& colsRegion = Rcpp::as<arma::uvec>(baseDM["colsRegion"]);
  const arma::uvec& colsTrans = Rcpp::as<arma::uvec>(baseDM["colsTrans"]);

  const arma::dvec& DTime = Rcpp::as<arma::dvec>(data["DTime"]);
  const Rcpp::NumericVector& Age = Rcpp::as<Rcpp::NumericVector>(data["Age"]);
  const Rcpp::NumericVector& Calendar = Rcpp::as<Rcpp::NumericVector>(data["Calendar"]);

  dm.cols(colsDTime - 1) = DTime + w;
  dm.cols(colsAge - 1) = Rcpp::as<arma::dmat>(Lspline(Age - w, knotsAge));
  dm.cols(colsCalendar - 1) = Rcpp::as<arma::dmat>(Lspline(Calendar - w, knotsCalendar));
  dm.cols(colsDTimeGender - 1) = dm.cols(colsGender - 1).eval().each_col() % dm.cols(colsDTime - 1);
  dm.cols(colsDTimeRegion - 1) = dm.cols(colsRegion - 1).eval().each_col() % dm.cols(colsDTime - 1);
  dm.cols(colsDTimeTrans - 1) = dm.cols(colsTrans - 1).eval().each_col() % dm.cols(colsDTime - 1);
  dm.cols(colsDTimeAge - 1) = dm.cols(colsAge - 1).eval().each_col() % dm.cols(colsDTime - 1);

  return dm;
}

} // namespace

#endif // _HivEstInfTime_UpdateCD4DesignMatrix_
