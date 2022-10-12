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

void R_init_nlmixr2omega(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"_nlmixr2omega_getBuiltinSize", (DL_FUNC) &_nlmixr2omega_getBuiltinSize, 0},
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
