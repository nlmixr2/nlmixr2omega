#define USE_FC_LEN_T
#define STRICT_R_HEADERS
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "omega.h"
#include "omegaR.h"

SEXP _nlmixr2omega_nlmixr2omegaNew(SEXP omeListSEXP, SEXP diagXformSEXP);
SEXP _nlmixr2omega_getTheta(SEXP pSEXP);
SEXP _nlmixr2omega_setTheta(SEXP pSEXP, SEXP thetaSEXP);
SEXP _nlmixr2omega_getCholOmegaInv(SEXP pSEXP);
SEXP _nlmixr2omega_getOmegaInv(SEXP pSEXP);
SEXP _nlmixr2omega_getdDomegaInv(SEXP pSEXP);
SEXP _nlmixr2omega_getCholOmega1(SEXP);
SEXP _nlmixr2omega_getOmegaR(SEXP pSEXP);
SEXP _nlmixr2omega_getCholOmega(SEXP);
SEXP _nlmixr2omega_getLogDetOMGAinv5(SEXP pSEXP);
SEXP _nlmixr2omega_nlmixr2omega_tr28(SEXP pSEXP);

void R_init_nlmixr2omega(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"_nlmixr2omega_getBuiltinSize", (DL_FUNC) &_nlmixr2omega_getBuiltinSize, 0},
    {"_nlmixr2omega_nlmixr2omegaNew", (DL_FUNC) &_nlmixr2omega_nlmixr2omegaNew, 2},
    {"_nlmixr2omega_getTheta", (DL_FUNC) &_nlmixr2omega_getTheta, 1},
    {"_nlmixr2omega_setTheta", (DL_FUNC) &_nlmixr2omega_setTheta, 2},
    {"_nlmixr2omega_getCholOmegaInv",
     (DL_FUNC) &_nlmixr2omega_getCholOmegaInv, 1},
    {"_nlmixr2omega_getOmegaInv",
     (DL_FUNC) &_nlmixr2omega_getOmegaInv, 1},
    {"_nlmixr2omega_getdDomegaInv",
     (DL_FUNC) &_nlmixr2omega_getdDomegaInv, 1},
    {"_nlmixr2omega_getOmegaR",
     (DL_FUNC) &_nlmixr2omega_getOmegaR, 1},
    {"_nlmixr2omega_getCholOmega",
     (DL_FUNC) &_nlmixr2omega_getCholOmega, 1},
    {"_nlmixr2omega_getLogDetOMGAinv5",
     (DL_FUNC) &_nlmixr2omega_getLogDetOMGAinv5, 1},
    {"_nlmixr2omega_nlmixr2omega_tr28",
     (DL_FUNC) &_nlmixr2omega_nlmixr2omega_tr28, 1},
    {NULL, NULL, 0} 
  };
  //R_RegisterCCallable("rxode2", "_rxode2_rxModelVars_", (DL_FUNC) &_rxode2_rxModelVars_);
  
  static const R_CMethodDef cMethods[] = {
    {NULL, NULL, 0, NULL}
  };

  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  
}

void R_unload_nlmixr2omega(DllInfo *info){
}
