//Generated from ::document() for 12 dimensions
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>
#include <Rmath.h>
#include "omega.h"

SEXP _nlmixr2omega_getBuiltinSize(void) {
  SEXP ret = PROTECT(Rf_allocVector(INTSXP, 1));
  INTEGER(ret)[0] = _nlmixr2omega_matSize();
  UNPROTECT(1);
  return ret;
}


