#ifndef _hivEstInfTime_UpdateVLDesignMatrix_
#define _hivEstInfTime_UpdateVLDesignMatrix_

#include "const.h"

namespace hivEstInfTime {

arma::dmat UpdateVLDesignMatrix(
  const Rcpp::List& baseDM,
  const Rcpp::DataFrame& data,
  const double& w
) {
  arma::dmat dm = Rcpp::as<arma::dmat>(baseDM["dm"]);

  if (dm.n_rows == 0) {
    return dm;
  }

  const arma::uvec& colsDTime = Rcpp::as<arma::uvec>(baseDM["colsDTime"]);
  const arma::uvec& colsGender = Rcpp::as<arma::uvec>(baseDM["colsGender"]);
  const arma::uvec& colsRegion = Rcpp::as<arma::uvec>(baseDM["colsRegion"]);
  const arma::uvec& colsTrans = Rcpp::as<arma::uvec>(baseDM["colsTrans"]);
  const arma::uvec& colsAge = Rcpp::as<arma::uvec>(baseDM["colsAge"]);
  const arma::uvec& colsLogDTime = Rcpp::as<arma::uvec>(baseDM["colsLogDTime"]);
  const arma::uvec& colsCalendar = Rcpp::as<arma::uvec>(baseDM["colsCalendar"]);
  const arma::uvec& colsDTimeGender = Rcpp::as<arma::uvec>(baseDM["colsDTimeGender"]);
  const arma::uvec& colsDTimeRegion = Rcpp::as<arma::uvec>(baseDM["colsDTimeRegion"]);
  const arma::uvec& colsDTimeTrans = Rcpp::as<arma::uvec>(baseDM["colsDTimeTrans"]);
  const arma::uvec& colsDTimeAge = Rcpp::as<arma::uvec>(baseDM["colsDTimeAge"]);
  const arma::uvec& colsLogDTimeGender = Rcpp::as<arma::uvec>(baseDM["colsLogDTimeGender"]);
  const arma::uvec& colsLogDTimeRegion = Rcpp::as<arma::uvec>(baseDM["colsLogDTimeRegion"]);
  const arma::uvec& colsLogDTimeTrans = Rcpp::as<arma::uvec>(baseDM["colsLogDTimeTrans"]);
  const arma::uvec& colsLogDTimeAge = Rcpp::as<arma::uvec>(baseDM["colsLogDTimeAge"]);

  const arma::dvec& DTime = Rcpp::as<arma::dvec>(data["DTime"]);
  const Rcpp::NumericVector& Age = Rcpp::as<Rcpp::NumericVector>(data["Age"]);
  const Rcpp::NumericVector& Calendar = Rcpp::as<Rcpp::NumericVector>(data["Calendar"]);

  dm.cols(colsDTime - 1) = DTime + w;
  dm.cols(colsAge - 1) = Rcpp::as<arma::dmat>(Lspline(Age - w, knotsAge));
  dm.cols(colsLogDTime - 1) = arma::log(DTime + w + 0.013);
  dm.cols(colsCalendar - 1) = Rcpp::as<arma::dmat>(Lspline(Calendar - w, knotsCalendar));
  dm.cols(colsDTimeGender - 1) = dm.cols(colsGender - 1).eval().each_col() % dm.cols(colsDTime - 1);
  dm.cols(colsDTimeRegion - 1) = dm.cols(colsRegion - 1).eval().each_col() % dm.cols(colsDTime - 1);
  dm.cols(colsDTimeTrans - 1) = dm.cols(colsTrans - 1).eval().each_col() % dm.cols(colsDTime - 1);
  dm.cols(colsDTimeAge - 1) = dm.cols(colsAge - 1).eval().each_col() % dm.cols(colsDTime - 1);
  dm.cols(colsLogDTimeGender - 1) = dm.cols(colsGender - 1).eval().each_col() % dm.cols(colsLogDTime - 1);
  dm.cols(colsLogDTimeRegion - 1) = dm.cols(colsRegion - 1).eval().each_col() % dm.cols(colsLogDTime - 1);
  dm.cols(colsLogDTimeTrans - 1) = dm.cols(colsTrans - 1).eval().each_col() % dm.cols(colsLogDTime - 1);
  dm.cols(colsLogDTimeAge - 1) = dm.cols(colsAge - 1).eval().each_col() % dm.cols(colsLogDTime - 1);

  return dm;
}

} // namespace

#endif // _hivEstInfTime_UpdateVLDesignMatrix_
