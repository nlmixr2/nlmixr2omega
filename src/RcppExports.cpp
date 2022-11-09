// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "../inst/include/nlmixr2omegaArma.h"

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// nlmixr2omegaNew
arma::vec nlmixr2omegaNew(List omeList, int diagXform);
RcppExport SEXP _nlmixr2omega_nlmixr2omegaNew(SEXP omeListSEXP, SEXP diagXformSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type omeList(omeListSEXP);
    Rcpp::traits::input_parameter< int >::type diagXform(diagXformSEXP);
    rcpp_result_gen = Rcpp::wrap(nlmixr2omegaNew(omeList, diagXform));
    return rcpp_result_gen;
END_RCPP
}
// getTheta
NumericVector getTheta(SEXP inSEXP);
RcppExport SEXP _nlmixr2omega_getTheta(SEXP inSEXPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type inSEXP(inSEXPSEXP);
    rcpp_result_gen = Rcpp::wrap(getTheta(inSEXP));
    return rcpp_result_gen;
END_RCPP
}
// setTheta
RObject setTheta(SEXP inSEXP, arma::vec theta);
RcppExport SEXP _nlmixr2omega_setTheta(SEXP inSEXPSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type inSEXP(inSEXPSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(setTheta(inSEXP, theta));
    return rcpp_result_gen;
END_RCPP
}
// getCholOmegaInv
arma::mat getCholOmegaInv(Rcpp::XPtr<_nlmixr2omega_full_omega> p);
RcppExport SEXP _nlmixr2omega_getCholOmegaInv(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<_nlmixr2omega_full_omega> >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(getCholOmegaInv(p));
    return rcpp_result_gen;
END_RCPP
}
// getOmegaInv
arma::mat getOmegaInv(Rcpp::XPtr<_nlmixr2omega_full_omega> p);
RcppExport SEXP _nlmixr2omega_getOmegaInv(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<_nlmixr2omega_full_omega> >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(getOmegaInv(p));
    return rcpp_result_gen;
END_RCPP
}
// getdDomegaInv
arma::mat getdDomegaInv(Rcpp::XPtr<_nlmixr2omega_full_omega> p);
RcppExport SEXP _nlmixr2omega_getdDomegaInv(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<_nlmixr2omega_full_omega> >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(getdDomegaInv(p));
    return rcpp_result_gen;
END_RCPP
}
// getCholOmega1
arma::mat getCholOmega1(Rcpp::XPtr<_nlmixr2omega_full_omega> p);
RcppExport SEXP _nlmixr2omega_getCholOmega1(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<_nlmixr2omega_full_omega> >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(getCholOmega1(p));
    return rcpp_result_gen;
END_RCPP
}
// getOmegaR
RObject getOmegaR(SEXP inSEXP);
RcppExport SEXP _nlmixr2omega_getOmegaR(SEXP inSEXPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type inSEXP(inSEXPSEXP);
    rcpp_result_gen = Rcpp::wrap(getOmegaR(inSEXP));
    return rcpp_result_gen;
END_RCPP
}
// getCholOmega
arma::mat getCholOmega(Rcpp::XPtr<_nlmixr2omega_full_omega> p);
RcppExport SEXP _nlmixr2omega_getCholOmega(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<_nlmixr2omega_full_omega> >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(getCholOmega(p));
    return rcpp_result_gen;
END_RCPP
}
// getLogDetOMGAinv5
double getLogDetOMGAinv5(Rcpp::XPtr<_nlmixr2omega_full_omega> p);
RcppExport SEXP _nlmixr2omega_getLogDetOMGAinv5(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<_nlmixr2omega_full_omega> >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(getLogDetOMGAinv5(p));
    return rcpp_result_gen;
END_RCPP
}
// nlmixr2omega_tr28
arma::vec nlmixr2omega_tr28(Rcpp::XPtr<_nlmixr2omega_full_omega> p);
RcppExport SEXP _nlmixr2omega_nlmixr2omega_tr28(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<_nlmixr2omega_full_omega> >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(nlmixr2omega_tr28(p));
    return rcpp_result_gen;
END_RCPP
}
