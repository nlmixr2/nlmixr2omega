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

// omegaFromR
Rcpp::XPtr<_nlmixr2omega_full_omega> omegaFromR(List omeList, int diagXform);
RcppExport SEXP _nlmixr2omega_omegaFromR(SEXP omeListSEXP, SEXP diagXformSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type omeList(omeListSEXP);
    Rcpp::traits::input_parameter< int >::type diagXform(diagXformSEXP);
    rcpp_result_gen = Rcpp::wrap(omegaFromR(omeList, diagXform));
    return rcpp_result_gen;
END_RCPP
}
