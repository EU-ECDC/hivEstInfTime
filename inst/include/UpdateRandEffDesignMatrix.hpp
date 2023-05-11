#ifndef _hivPlatform_UpdateRandEffDesignMatrix_
#define _hivPlatform_UpdateRandEffDesignMatrix_

namespace hivPlatform {

arma::dmat UpdateRandEffDesignMatrix(
  const Rcpp::List& baseDM,
  const Rcpp::DataFrame& data,
  const double& w
) {
  arma::dmat dm = Rcpp::as<arma::dmat>(baseDM["dm"]);

  if (dm.n_rows == 0) {
    return dm;
  }

  const arma::uvec& colsDTimeConsc = Rcpp::as<arma::uvec>(baseDM["colsDTimeConsc"]);
  const arma::uvec& colsDTimeConsr = Rcpp::as<arma::uvec>(baseDM["colsDTimeConsr"]);
  const arma::uvec& colsLogDTimeConsr = Rcpp::as<arma::uvec>(baseDM["colsLogDTimeConsr"]);

  const arma::dvec& Consc = Rcpp::as<arma::dvec>(data["Consc"]);
  const arma::dvec& Consr = Rcpp::as<arma::dvec>(data["Consr"]);
  const arma::dvec& DTime = Rcpp::as<arma::dvec>(data["DTime"]);

  dm.cols(colsDTimeConsc - 1) = (DTime + w) % Consc;
  dm.cols(colsDTimeConsr - 1) = (DTime + w) % Consr;
  dm.cols(colsLogDTimeConsr - 1) = arma::log(DTime + w + 0.013) % Consr;

  return dm;
}

} // namespace

#endif // _hivPlatform_UpdateRandEffDesignMatrix_
