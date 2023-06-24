#ifndef _HivEstInfTime_GetLogMVNPdf_
#define _HivEstInfTime_GetLogMVNPdf_

#include "const.h"

namespace HivEstInfTime {

static double GetLogMVNPdf(
  const arma::dvec& x,
  const arma::dvec& mu,
  const arma::dmat& sigma
) {
  const arma::dmat rooti = arma::inv(arma::trimatu(arma::chol(sigma)));
  const double otherTerms = arma::sum(arma::log(rooti.diag())) + x.n_elem * halfLog2Pi;

  arma::dvec z = x - mu;
  for (size_t i = rooti.n_cols; i != -1; --i) {
    double tmp(0.0);
    for (size_t j = 0; j <= i; ++j) {
      tmp += rooti.at(j, i) * z[j];
    }
    z[i] = tmp;
  }

  return otherTerms - 0.5 * arma::dot(z, z);
}

} // namespace

#endif // _HivEstInfTime_GetLogMVNPdf_
