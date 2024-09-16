#include "header.h"

// [[Rcpp::export]]
Rcpp::NumericVector Lspline(
  const Rcpp::NumericVector& x,
  const Rcpp::NumericVector& knots
) {
  return hivEstInfTime::Lspline(x, knots);
};


// [[Rcpp::export]]
double GetLogMVNPdf(
  const arma::dvec& x,
  const arma::dvec& mu,
  const arma::dmat& sigma
) {
  return hivEstInfTime::GetLogMVNPdf(x, mu, sigma);
};

// [[Rcpp::export]]
double PostW(
  const double& w,
  const arma::dvec& y,
  const arma::dmat& xAIDS,
  const double& maxDTime,
  const arma::dmat& betaAIDS,
  const double& kappa,
  const arma::dmat& bFE,
  const arma::dmat& varCovRE,
  const Rcpp::List& baseCD4DM,
  const Rcpp::DataFrame& fxCD4Data,
  const Rcpp::List& baseVLDM,
  const Rcpp::DataFrame& fxVLData,
  const Rcpp::List& baseRandEffDM,
  const Rcpp::DataFrame& fzData,
  const arma::dmat& err
) {
  hivEstInfTime::PostW f(
    y, xAIDS, maxDTime, betaAIDS, kappa, bFE, varCovRE, baseCD4DM, fxCD4Data, baseVLDM, fxVLData,
    baseRandEffDM, fzData, err
  );

  return f(w);
};

// [[Rcpp::export]]
Rcpp::NumericVector VPostW(
  const arma::dvec& w,
  const arma::dvec& y,
  const arma::dmat& xAIDS,
  const double& maxDTime,
  const arma::dmat& betaAIDS,
  const double& kappa,
  const arma::dmat& bFE,
  const arma::dmat& varCovRE,
  const Rcpp::List& baseCD4DM,
  const Rcpp::DataFrame& fxCD4Data,
  const Rcpp::List& baseVLDM,
  const Rcpp::DataFrame& fxVLData,
  const Rcpp::List& baseRandEffDM,
  const Rcpp::DataFrame& fzData,
  const arma::dmat& err
) {
  hivEstInfTime::PostW f(
    y, xAIDS, maxDTime, betaAIDS, kappa, bFE, varCovRE, baseCD4DM, fxCD4Data, baseVLDM, fxVLData,
    baseRandEffDM, fzData, err
  );

  Rcpp::NumericVector out(w.n_elem);
  for (size_t i = 0; i != w.n_elem; ++i) {
    out[i] = f(w[i]);
  }
  out.attr("dim") = R_NilValue;

  return out;
};

// [[Rcpp::export]]
Rcpp::List IntegratePostW(
  const double& lower,
  const double& upper,
  const arma::dvec& y,
  const arma::dmat& xAIDS,
  const double& maxDTime,
  const arma::dmat& betaAIDS,
  const double& kappa,
  const arma::dmat& bFE,
  const arma::dmat& varCovRE,
  const Rcpp::List& baseCD4DM,
  const Rcpp::DataFrame& fxCD4Data,
  const Rcpp::List& baseVLDM,
  const Rcpp::DataFrame& fxVLData,
  const Rcpp::List& baseRandEffDM,
  const Rcpp::DataFrame& fzData,
  const arma::dmat& err
) {
  hivEstInfTime::PostW f(
    y, xAIDS, maxDTime, betaAIDS, kappa, bFE, varCovRE, baseCD4DM, fxCD4Data, baseVLDM, fxVLData,
    baseRandEffDM, fzData, err
  );

  double errEst;
  int errCode;
  const double res = Numer::integrate(
    f,
    lower,
    upper,
    errEst,
    errCode,
    100,
    0.0001220703,
    0.0001220703,
    Numer::Integrator<double>::GaussKronrod15
  );
  return Rcpp::List::create(
    Rcpp::Named("value") = res,
    Rcpp::Named("errorEstimate") = errEst,
    Rcpp::Named("errorCode") = errCode
  );
};


// [[Rcpp::export]]
Rcpp::List IntegrateMeanPostW(
  const double& lower,
  const double& upper,
  const arma::dvec& y,
  const arma::dmat& xAIDS,
  const double& maxDTime,
  const arma::dmat& betaAIDS,
  const double& kappa,
  const arma::dmat& bFE,
  const arma::dmat& varCovRE,
  const Rcpp::List& baseCD4DM,
  const Rcpp::DataFrame& fxCD4Data,
  const Rcpp::List& baseVLDM,
  const Rcpp::DataFrame& fxVLData,
  const Rcpp::List& baseRandEffDM,
  const Rcpp::DataFrame& fzData,
  const arma::dmat& err
) {
  hivEstInfTime::MeanPostW f(
    y, xAIDS, maxDTime, betaAIDS, kappa, bFE, varCovRE, baseCD4DM, fxCD4Data, baseVLDM, fxVLData,
    baseRandEffDM, fzData, err
  );

  double errEst;
  int errCode;
  const double res = Numer::integrate(
    f,
    lower,
    upper,
    errEst,
    errCode,
    100,
    0.0001220703,
    0.0001220703,
    Numer::Integrator<double>::GaussKronrod15
  );
  return Rcpp::List::create(
    Rcpp::Named("value") = res,
    Rcpp::Named("errorEstimate") = errEst,
    Rcpp::Named("errorCode") = errCode
  );
};
