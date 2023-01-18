#define STRICT_R_HEADER
#define ARMA_WARN_LEVEL 1
#define ARMA_DONT_USE_OPENMP // Known to cause speed problems
// #ifdef _OPENMP
// #include <omp.h>
// #endif
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <stdlib.h>

using namespace arma;
using namespace Rcpp;

#include "../inst/include/nlmixr2omegaArma.h"
#include "omegaR.h"
#include "omega.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext("nlmixr2omega", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

void _nlmixr2omegaAssignTheta(_nlmixr2omega_ind_omega *ome, arma::vec theta) {
  arma::vec cur = theta;
  ome->theta = cur;
  // sum_{k=1}^{n} k = n*(n+1)/2
  int d0 = theta.size();
  ome->dim = 0.5*sqrt(1.0 + d0 * 8.0) - 0.5;
  // n^2 + n = d0*2
  ome->cholOmegaInvBool = false;
  ome->omegaInvBool = false;
  ome->dOmegaInvBool = false;
  ome->dDomegaInvBool = false;
  ome->cholOmega1Bool = false;
  ome->omegaBool = false;
  ome->cholOmegaBool = false;
  ome->logDetOMGAinv5Bool = false;
  ome->tr28Bool = false;
  ome->omega47Bool = false;
}

arma::mat nlmixr2omega_cholOmegaInv(_nlmixr2omega_ind_omega *ome) {
  if (ome->cholOmegaInvBool && ome->cholOmegaInvMat.n_elem != 0) return ome->cholOmegaInvMat;
  int ts = ome->theta.size();
  int tn  = 0;
  arma::mat ret(ome->dim, ome->dim, arma::fill::zeros);
  ome->cFun(&(ome->dim), ome->theta.memptr(), &ts, &tn, ret.memptr());
  ome->cholOmegaInvMat = ret;
  ome->cholOmegaInvBool = true;
  return ome->cholOmegaInvMat;
}

arma::mat nlmixr2omega_omegaInv(_nlmixr2omega_ind_omega *ome) {
  if (ome->omegaInvBool && ome->omegaInvMat.n_elem != 0) return ome->omegaInvMat;
  int ts = ome->theta.size();
  int tn  = -1;
  arma::mat ret(ome->dim, ome->dim, arma::fill::zeros);
  ome->cFun(&(ome->dim), ome->theta.memptr(), &ts, &tn, ret.memptr());
  ome->omegaInvMat = ret;
  ome->omegaInvBool = true;
  return ome->omegaInvMat;
}

arma::cube nlmixr2omega_dOmegaInv(_nlmixr2omega_ind_omega *ome) {
  if (ome->dDomegaInvBool && ome->dOmegaInvCube.n_elem != 0) return ome->dOmegaInvCube;
  int ts = ome->theta.size();
  arma::cube ret(ome->dim, ome->dim, ts, arma::fill::zeros);
  int tn;
  for (int i = 0; i < ts; ++i) {
    tn = i+1;
    ome->cFun(&(ome->dim), ome->theta.memptr(), &ts, &tn, ret.slice(i).memptr());
  }
  ome->dOmegaInvCube = ret;
  ome->dDomegaInvBool = true;
  return ome->dOmegaInvCube;
}

arma::mat nlmixr2omega_dDomegaInv(_nlmixr2omega_ind_omega *ome) {
  if (ome->dDomegaInvBool && ome->dDomegaInvMat.n_elem != 0) return (ome->dDomegaInvMat);
  int ts = ome->theta.size();
  // dim  x ntheta ; "d.D.omegaInv"
  arma::mat ret(ome->dim, ts, arma::fill::zeros);
  int tn;
  for (int i = 0; i < ts; ++i) {
    tn = -2 - (i+1);
    ome->cFun(&(ome->dim), ome->theta.memptr(), &ts, &tn, ret.memptr() + i*ome->dim);
  }
  ome->dDomegaInvMat = ret;
  ome->dDomegaInvBool = true;
  return ome->dDomegaInvMat;
}


arma::mat rxToCholOmega(arma::mat cholMat){
  arma::mat cholO;
  bool success;
  try {
    success = inv(cholO, trimatu(cholMat));
    if (success) return cholO;
    success = inv(cholO, cholMat);
    if (success) return cholO;
    stop(_("can not invert in 'rxToCholOmega'"));
  } catch (...) {
    success = inv(cholO, cholMat);
    if (success) return cholO;
    stop(_("can not invert in 'rxToCholOmega'"));
  }
  // should not get here.
  return cholO;
}

arma::mat nlmixr2omega_cholOmega1(_nlmixr2omega_ind_omega *ome) {
  if (ome->cholOmega1Bool && ome->cholOmega1Mat.n_elem != 0) return ome->cholOmega1Mat;
  arma::mat coi =  nlmixr2omega_cholOmegaInv(ome);
  coi = rxToCholOmega(coi);
  ome->cholOmega1Bool = true;
  ome->cholOmega1Mat = coi;
  return ome->cholOmega1Mat;
}

arma::mat nlmixr2omega_omega(_nlmixr2omega_ind_omega *ome) {
  if (ome->omegaBool && ome->omegaMat.n_elem != 0) return ome->omegaMat;
  arma::mat u1 = nlmixr2omega_cholOmega1(ome);
  arma::mat omega = u1*trans(u1);
  ome->omegaBool = true;
  ome->omegaMat = omega;
  return ome->omegaMat;
}

arma::mat nlmixr2omega_cholOmega(_nlmixr2omega_ind_omega *ome) {
  if (ome->cholOmegaBool) return ome->cholOmegaMat;
  ome->cholOmegaBool = true;
  ome->cholOmegaMat = chol(nlmixr2omega_omega(ome));
  return ome->cholOmegaMat;
}

double nlmixr2omega_logDetOMGAinv5(_nlmixr2omega_ind_omega *ome) {
  if (ome->logDetOMGAinv5Bool) return ome->logDetOMGAinv5double;
  arma::mat coi =  nlmixr2omega_cholOmegaInv(ome);
  arma::vec diag = coi.diag();
  arma::vec ldiag = log(diag);
  ome->logDetOMGAinv5double = sum(ldiag);
  ome->logDetOMGAinv5Bool = true;
  return ome->logDetOMGAinv5double;
}

arma::vec nlmixr2omega_tr28(_nlmixr2omega_ind_omega *ome) {
  // 1/2*tr(d(Omega^-1)*Omega);
  if (ome->tr28Bool) return ome->tr28vec;
  unsigned int ts = ome->theta.size();
  arma::mat omega = nlmixr2omega_omega(ome);
  arma::cube dOmegaInv = nlmixr2omega_dOmegaInv(ome);
  arma::vec tr28(ts);
  for (unsigned int i = ts; i--;){
    arma::mat cur = (dOmegaInv.slice(i) * omega);
    tr28[i] =0.5*sum(cur.diag());
  }
  ome->tr28vec = tr28;
  ome->tr28Bool = true;
  return ome->tr28vec;
}

arma::cube nlmixr2omega_omega47(_nlmixr2omega_ind_omega *ome) {
  if (ome->omega47Bool) return ome->omega47Cube;
  arma::mat cholO =  nlmixr2omega_cholOmegaInv(ome);
  arma::cube dOmegaInv = nlmixr2omega_dOmegaInv(ome);
  arma::mat cEta = zeros(ome->dim,1);
  arma::mat c;
  unsigned int ts = ome->theta.size();
  arma::cube ret(ome->dim, ome->dim, ts, arma::fill::none);
  for (unsigned int i = ts; i--;){
    c = dOmegaInv.slice(i);
    arma::mat prodI(ome->dim, ome->dim, arma::fill::none);
    for (unsigned int j = ome->dim; j--;){
      cEta(j,0) = 1;
      prodI.col(j) = c*cEta;
      cEta(j,0) = 0;
    }
    ret.slice(i) = prodI;
  }
  ome->omega47Cube = ret;
  ome->omega47Bool = true;
  return ome->omega47Cube;
}

// Inverts the matrix
arma::mat nlmixr2omega_inv(arma::mat &smatrix) {
  // Invert matrix using RcppArmadillo.
  arma::mat imat;
  bool success;
  success = arma::inv(imat, smatrix);
  if (!success){
    imat = arma::pinv(smatrix);
    //REprintf(_("matrix seems singular; Using pseudo-inverse\n"));
  }
  return imat;
}

// Get the cholesky decomposition of the inverse
arma::mat nlmixr2omega_cholInv(arma::mat &mat) {
  return arma::chol(nlmixr2omega_inv(mat));
}

void nlmixr2omega_iniOmeStruct(_nlmixr2omega_ind_omega *ome,
                               arma::mat &mat, int diagXform) {
  arma::mat in =  nlmixr2omega_cholInv(mat);
  switch (diagXform) {
  case nlmixr2omega_sqrt:
    in.diag() = sqrt(in.diag());
    ome->cFun = _nlmixr2omega_mat_sqrt;
    break;
  case nlmixr2omega_log:
    in.diag() = log(in.diag());
    ome->cFun = _nlmixr2omega_mat_log;
  default:
    ome->cFun = _nlmixr2omega_mat;
    break;
  }
  arma::vec theta = in(trimatu_ind(size(in)));
  _nlmixr2omegaAssignTheta(ome, theta);
}

arma::vec _nlmixr2omega_full_getTheta_(_nlmixr2omega_full_omega *fome) {
  if (fome->nomes == 0) {
    arma::vec ret;
    return ret;
  } else if (fome->nomes == 1) {
    return fome->omes[0].theta;
  } else {
    arma::vec theta = fome->omes[0].theta;
    for (int i = 1; i < fome->nomes; ++i) {
      theta = join_cols(theta, fome->omes[i].theta);
    }
    return theta;
  }
}


arma::vec nlmixr2omegaNewVec(_nlmixr2omega_full_omega *fullPtr, int diagXform) {
  arma::vec ret(4);
  ret[0] = fullPtr->nomes;
  ret[1] = fullPtr->nTotTheta;
  ret[2] = fullPtr->nTotDim;
  ret[3] = diagXform;
  for (int i = 0; i < fullPtr->nomes; ++i) {
    _nlmixr2omega_ind_omega *ome = &(fullPtr->omes[i]);
    arma::vec n(1);
    n[0] = ome->dim;
    ret= join_cols(ret, n);
  }
  arma::vec theta = _nlmixr2omega_full_getTheta_(fullPtr);
  free(fullPtr->omes);
  fullPtr->omes = NULL;
  //arma::vec theta = _nlmixr2omega_full_getTheta_(fullPtr);
  return join_cols(ret, theta);
}

//[[Rcpp::export]]
arma::vec nlmixr2omegaNew(List omeList, int diagXform) {
  _nlmixr2omega_full_omega full;
  _nlmixr2omega_full_omega *fullPtr = &full;
  fullPtr->nomes = omeList.size();
  if (fullPtr->omes != NULL) free(fullPtr->omes);
  fullPtr->omes = (_nlmixr2omega_ind_omega*)malloc(fullPtr->nomes*sizeof(_nlmixr2omega_ind_omega));
  fullPtr->nTotTheta = 0;
  fullPtr->nTotDim = 0;
  for (int i = 0; i < fullPtr->nomes; ++i) {
    _nlmixr2omega_ind_omega *ome = &(fullPtr->omes[i]);
    arma::mat cur = as<arma::mat>(omeList[i]);
    nlmixr2omega_iniOmeStruct(ome, cur, diagXform);
    fullPtr->nTotTheta += ome->theta.size();
    fullPtr->nTotDim += ome->dim;
  }
  return nlmixr2omegaNewVec(fullPtr, diagXform);
}

_nlmixr2omega_full_omega nlmixr2omega_full_Create(arma::vec in) {
  _nlmixr2omega_full_omega full;
  _nlmixr2omega_full_omega *fullPtr = &full;
  double *ptr = in.memptr();
  fullPtr->nomes = (int)(in[0]);
  fullPtr->nTotTheta = (int)(in[1]);
  fullPtr->nTotDim = (int)(in[2]);
  int diagXform = (int)(in[3]);
  ptr += 4;
  fullPtr->omes = (_nlmixr2omega_ind_omega*)malloc(fullPtr->nomes*sizeof(_nlmixr2omega_ind_omega));
  for (int i = 0; i < fullPtr->nomes; ++i) {
    _nlmixr2omega_ind_omega *ome = &(fullPtr->omes[i]);
    ome->cholOmegaInvBool = false;
    ome->omegaInvBool = false;
    ome->dOmegaInvBool = false;
    ome->dDomegaInvBool = false;
    ome->cholOmega1Bool = false;
    ome->omegaBool = false;
    ome->cholOmegaBool = false;
    ome->logDetOMGAinv5Bool=false;
    ome->tr28Bool = false;
    ome->omega47Bool = false;
    switch (diagXform) {
    case nlmixr2omega_sqrt:
      ome->cFun = _nlmixr2omega_mat_sqrt;
      break;
    case nlmixr2omega_log:
      ome->cFun = _nlmixr2omega_mat_log;
    default:
      ome->cFun = _nlmixr2omega_mat;
      break;
    }
    ome->dim = ptr[0];
    ptr++;
  }
  for (int i = 0; i < fullPtr->nomes; ++i) {
    _nlmixr2omega_ind_omega *ome = &(fullPtr->omes[i]);
    int ntheta = 0.5*(ome->dim)*(ome->dim+1.0);
    arma::vec theta(ntheta);
    std::copy(ptr, ptr+ntheta, theta.begin());
    ptr += ntheta;
    _nlmixr2omegaAssignTheta(ome, theta);
  }
  return full;
}

_nlmixr2omega_full_omega omegaFromRgetFullOmegaFromSexp(RObject inSEXP) {
  if (!Rf_inherits(inSEXP, "nlmixr2omega")) {
    stop("needs to be class 'nlmixr2omega'");
  }
  List cur = as<List>(inSEXP);
  arma::vec in = as<arma::vec>(cur[0]);
  return nlmixr2omega_full_Create(in);
}

int omegaFromRgetDiagXfrom(RObject inSEXP) {
  if (!Rf_inherits(inSEXP, "nlmixr2omega")) {
    stop("needs to be class 'nlmixr2omega'");
  }
  List cur = as<List>(inSEXP);
  arma::vec in = as<arma::vec>(cur[0]);
  return (int)(in[3]);
}


//[[Rcpp::export]]
NumericVector getTheta(RObject inSEXP) {
  _nlmixr2omega_full_omega p = omegaFromRgetFullOmegaFromSexp(inSEXP);
  return wrap(_nlmixr2omega_full_getTheta_(&p));
}

void _nlmixr2omega_full_setTheta_(_nlmixr2omega_full_omega *fome,
                                  arma::vec theta) {
  if (fome->nTotTheta != theta.size()) {
    stop("incompatible size with this omega structure");
  }
  double *ptr = theta.memptr();
  for (int i = 0; i < fome->nomes; ++i) {
    _nlmixr2omega_ind_omega *ome = &(fome->omes[i]);
    arma::vec curTheta = arma::vec(ptr, ome->theta.size());
      _nlmixr2omegaAssignTheta(ome, curTheta);
      ptr += curTheta.size();
  }
}

//[[Rcpp::export]]
RObject setTheta(RObject inSEXP, arma::vec theta) {
  _nlmixr2omega_full_omega p = omegaFromRgetFullOmegaFromSexp(inSEXP);
  _nlmixr2omega_full_setTheta_(&p, theta);
  List ret = List::create(nlmixr2omegaNewVec(&p, omegaFromRgetDiagXfrom(inSEXP)),
                          VECTOR_ELT(inSEXP, 1));
  ret.attr("class") = "nlmixr2omega";
  return ret;
}

arma::mat _nlmixr2omega_full_cholOmegaInv(_nlmixr2omega_full_omega *fome) {
  arma::mat ret(fome->nTotDim, fome->nTotDim, arma::fill::zeros);
  int curBlock = 0;
  for (int i = 0; i < fome->nomes; ++i) {
    _nlmixr2omega_ind_omega *ome = &(fome->omes[i]);
    int curDim = ome->dim;
    ret.submat(curBlock, curBlock,
               curBlock+curDim-1, curBlock+curDim-1) =
      nlmixr2omega_cholOmegaInv(ome);
    curBlock += curDim;
  }
  return ret;
}

//[[Rcpp::export]]
RObject getCholOmegaInv(RObject inSEXP) {
  _nlmixr2omega_full_omega p = omegaFromRgetFullOmegaFromSexp(inSEXP);
  NumericMatrix ret = wrap(_nlmixr2omega_full_cholOmegaInv(&p));
  List v = as<List>(inSEXP);
  ret.attr("dimnames") = v[1];
  return ret;
}

arma::mat _nlmixr2omega_full_omegaInv(_nlmixr2omega_full_omega *fome) {
  arma::mat ret(fome->nTotDim, fome->nTotDim, arma::fill::zeros);
  int curBlock = 0;
  for (int i = 0; i < fome->nomes; ++i) {
    _nlmixr2omega_ind_omega *ome = &(fome->omes[i]);
    int curDim = ome->dim;
    ret.submat(curBlock, curBlock,
               curBlock+curDim-1, curBlock+curDim-1) =
      nlmixr2omega_omegaInv(ome);
    curBlock += curDim;
  }
  return ret;
}

//[[Rcpp::export]]
RObject getOmegaInv(RObject inSEXP) {
  _nlmixr2omega_full_omega p = omegaFromRgetFullOmegaFromSexp(inSEXP);
  NumericMatrix ret = wrap(_nlmixr2omega_full_omegaInv(&p));
  List v = as<List>(inSEXP);
  ret.attr("dimnames") = v[1];
  return ret;
}

arma::mat _nlmixr2omega_full_dDomegaInv(_nlmixr2omega_full_omega *fome) {
  arma::mat ret(fome->nTotDim, fome->nTotDim, arma::fill::zeros);
  int curBlock = 0;
  for (int i = 0; i < fome->nomes; ++i) {
    _nlmixr2omega_ind_omega *ome = &(fome->omes[i]);
    int curDim = ome->dim;
    ret.submat(curBlock, curBlock,
               curBlock+curDim-1, curBlock+curDim-1) =
      nlmixr2omega_dDomegaInv(ome);
    curBlock += curDim;
  }
  print(wrap(ret));
  return ret;
}

//[[Rcpp::export]]
RObject getdDomegaInv(RObject inSEXP) {
  _nlmixr2omega_full_omega p = omegaFromRgetFullOmegaFromSexp(inSEXP);
  RObject ret = wrap(_nlmixr2omega_full_dDomegaInv(&p));
  //ret.attr("dimnames") = PROTECT(VECTOR_ELT(inSEXP, 1));
  //UNPROTECT(1);
  return ret;
}

arma::mat _nlmixr2omega_full_cholOmega1(_nlmixr2omega_full_omega *fome) {
  arma::mat ret(fome->nTotDim, fome->nTotDim, arma::fill::zeros);
  int curBlock = 0;
  for (int i = 0; i < fome->nomes; ++i) {
    _nlmixr2omega_ind_omega *ome = &(fome->omes[i]);
    int curDim = ome->dim;
    ret.submat(curBlock, curBlock,
               curBlock+curDim-1, curBlock+curDim-1) =
      nlmixr2omega_cholOmega1(ome);
    curBlock += curDim;
  }
  return ret;
}

//[[Rcpp::export]]
RObject getCholOmega1(RObject inSEXP) {
  _nlmixr2omega_full_omega p = omegaFromRgetFullOmegaFromSexp(inSEXP);
  RObject ret = wrap(_nlmixr2omega_full_cholOmega1(&p));
  List v = as<List>(inSEXP);
  ret.attr("dimnames") = v[1];
  return ret;
}

arma::mat _nlmixr2omega_full_omegaR(_nlmixr2omega_full_omega *fome) {
  arma::mat ret(fome->nTotDim, fome->nTotDim, arma::fill::zeros);
  int curBlock = 0;
  for (int i = 0; i < fome->nomes; ++i) {
    _nlmixr2omega_ind_omega *ome = &(fome->omes[i]);
    int curDim = ome->dim;
    arma::mat cur = nlmixr2omega_omega(ome);
    ret.submat(curBlock, curBlock,
               curBlock+curDim-1, curBlock+curDim-1) = cur;
    curBlock += curDim;
  }
  return ret;
}

//[[Rcpp::export]]
RObject getOmegaR(RObject inSEXP) {
  _nlmixr2omega_full_omega p = omegaFromRgetFullOmegaFromSexp(inSEXP);
  NumericMatrix ret = wrap(_nlmixr2omega_full_omegaR(&p));
  List v = as<List>(inSEXP);
  ret.attr("dimnames") = v[1];
  return ret;
}

arma::mat _nlmixr2omega_full_cholOmega(_nlmixr2omega_full_omega *fome) {
  arma::mat ret(fome->nTotDim, fome->nTotDim, arma::fill::zeros);
  int curBlock = 0;
  for (int i = 0; i < fome->nomes; ++i) {
    _nlmixr2omega_ind_omega *ome = &(fome->omes[i]);
    int curDim = ome->dim;
    ret.submat(curBlock, curBlock,
               curBlock+curDim, curBlock+curDim) =
      nlmixr2omega_cholOmega(ome);
    curBlock += curDim;
  }
  return ret;
}

//[[Rcpp::export]]
arma::mat getCholOmega(Rcpp::XPtr<_nlmixr2omega_full_omega> p) {
  _nlmixr2omega_full_omega* v = p.get();
  return _nlmixr2omega_full_omegaR(v);
}

double _nlmixr2omega_full_logDetOMGAinv5(_nlmixr2omega_full_omega *fome) {
  double ret = 0.0;
  for (int i = 0; i < fome->nomes; ++i) {
    _nlmixr2omega_ind_omega *ome = &(fome->omes[i]);
    ret += nlmixr2omega_logDetOMGAinv5(ome);
  }
  return ret;
}

//[[Rcpp::export]]
double getLogDetOMGAinv5(Rcpp::XPtr<_nlmixr2omega_full_omega> p) {
  _nlmixr2omega_full_omega* v = p.get();
  return _nlmixr2omega_full_logDetOMGAinv5(v);
}

arma::vec _nlmixr2omega_full_tr28(_nlmixr2omega_full_omega *fome) {
  arma::vec ret(fome->nTotTheta, arma::fill::none);
  double *ptr = ret.memptr();
  for (int i = 0; i < fome->nomes; ++i) {
    _nlmixr2omega_ind_omega *ome = &(fome->omes[i]);
    arma::vec tr28 = nlmixr2omega_tr28(ome);
    std::copy(tr28.begin(), tr28.end(), ptr);
    ptr += tr28.size();
  }
  return ret;
}

//[[Rcpp::export]]
arma::vec nlmixr2omega_tr28(Rcpp::XPtr<_nlmixr2omega_full_omega> p) {
  _nlmixr2omega_full_omega* v = p.get();
  return _nlmixr2omega_full_tr28(v);
}
