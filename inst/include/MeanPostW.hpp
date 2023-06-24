#ifndef _HivEstInfTime_MeanPostW_
#define _HivEstInfTime_MeanPostW_

#include "Lspline.hpp"
#include "UpdateCD4DesignMatrix.hpp"
#include "UpdateVLDesignMatrix.hpp"
#include "UpdateRandEffDesignMatrix.hpp"
#include "GetLogMVNPdf.hpp"

namespace HivEstInfTime {

class MeanPostW: public Numer::Func
{
private:
  const arma::dvec& y;
  const arma::dmat& xAIDS;
  const double& maxDTime;
  const arma::dmat& betaAIDS;
  const double& kappa;
  const arma::dmat& bFE;
  const arma::dmat& varCovRE;
  const Rcpp::List& baseCD4DM;
  const Rcpp::DataFrame& fxCD4Data;
  const Rcpp::List& baseVLDM;
  const Rcpp::DataFrame& fxVLData;
  const Rcpp::List& baseRandEffDM;
  const Rcpp::DataFrame& fzData;
  const arma::dmat& err;
public:
  MeanPostW(
    const arma::dvec& y_,
    const arma::dmat& xAIDS_,
    const double& maxDTime_,
    const arma::dmat& betaAIDS_,
    const double& kappa_,
    const arma::dmat& bFE_,
    const arma::dmat& varCovRE_,
    const Rcpp::List& baseCD4DM_,
    const Rcpp::DataFrame& fxCD4Data_,
    const Rcpp::List& baseVLDM_,
    const Rcpp::DataFrame& fxVLData_,
    const Rcpp::List& baseRandEffDM_,
    const Rcpp::DataFrame& fzData_,
    const arma::dmat& err_
  ) :
    y(y_),
    xAIDS(xAIDS_),
    maxDTime(maxDTime_),
    betaAIDS(betaAIDS_),
    kappa(kappa_),
    bFE(bFE_),
    varCovRE(varCovRE_),
    baseCD4DM(baseCD4DM_),
    fxCD4Data(fxCD4Data_),
    baseVLDM(baseVLDM_),
    fxVLData(fxVLData_),
    baseRandEffDM(baseRandEffDM_),
    fzData(fzData_),
    err(err_)
    {}

  double operator()(const double& w) const
  {
    arma::dmat xAIDSnew(xAIDS);
    xAIDSnew(0, 2) -= w;

    const double lambda = arma::as_scalar(arma::exp(xAIDSnew * betaAIDS));

    const arma::dmat xCD4 = UpdateCD4DesignMatrix(baseCD4DM, fxCD4Data, w);
    const arma::dmat xVL = UpdateVLDesignMatrix(baseVLDM, fxVLData, w);
    const arma::dmat z = UpdateRandEffDesignMatrix(baseRandEffDM, fzData, w);

    const arma::dmat x = arma::join_cols(
      arma::join_rows(arma::zeros<arma::dmat>(xVL.n_rows, xCD4.n_cols), xVL),
      arma::join_rows(xCD4, arma::zeros<arma::dmat>(xCD4.n_rows, xVL.n_cols))
    );

    const arma::dvec mu = x * bFE;
    const arma::dmat var = z * (varCovRE * z.t()) + err;

    return w * std::exp(GetLogMVNPdf(y, mu, var) - lambda * std::pow(w + maxDTime, kappa));
  }
};

} // namespace

#endif // _HivEstInfTime_MeanPostW_
