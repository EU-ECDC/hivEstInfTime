#ifndef _HivEstInfTime_Lspline_
#define _HivEstInfTime_Lspline_

namespace hivEstInfTime {

Rcpp::NumericVector Lspline(
  const Rcpp::NumericVector& x,
  const Rcpp::NumericVector& knots
) {

  const size_t n = x.length();
  const size_t nvars = knots.length() + 1;
  Rcpp::NumericMatrix rval(n, nvars);
  rval(Rcpp::_, 0) = Rcpp::pmin(x, knots[0]);
  rval(Rcpp::_, nvars - 1) = Rcpp::pmax(x, knots[nvars - 2]) - knots[nvars - 2];
  if (nvars > 2) {
    for (size_t i = 1; i != nvars - 1; ++i) {
      rval(Rcpp::_, i) = Rcpp::pmax(Rcpp::pmin(x, knots[i]), knots[i - 1]) - knots[i - 1];
    }
  }

  return rval;
}

} // namespace

#endif // _hivEstInfTime_Lspline_
