#ifndef _hivPlatform_consts_
#define _hivPlatform_consts_

static const double halfLog2Pi = -0.5 * std::log(2.0 * M_PI);
static const Rcpp::NumericVector knotsAge = Rcpp::NumericVector::create(25, 35, 45);
static const Rcpp::NumericVector knotsCalendar = Rcpp::NumericVector::create(16, 22);

#endif // _hivPlatform_consts_
