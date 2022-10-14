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
  ome->theta = theta;
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
  if (ome->cholOmegaInvBool) return ome->cholOmegaInvMat;
  int ts = ome->theta.size();
  int tn  = 0;
  arma::mat ret(ome->dim, ome->dim, arma::fill::zeros);
  ome->cFun(&(ome->dim), ome->theta.memptr(), &ts, &tn, ret.memptr());
  ome->cholOmegaInvMat = ret;
  ome->cholOmegaInvBool = true;
  return ome->cholOmegaInvMat;
}

arma::mat nlmixr2omega_omegaInv(_nlmixr2omega_ind_omega *ome) {
  if (ome->omegaInvBool) return ome->omegaInvMat;
  int ts = ome->theta.size();
  int tn  = -1;
  arma::mat ret(ome->dim, ome->dim, arma::fill::zeros);
  ome->cFun(&(ome->dim), ome->theta.memptr(), &ts, &tn, ret.memptr());
  ome->omegaInvMat = ret;
  ome->omegaInvBool = true;
  return ome->omegaInvMat;
}

arma::cube nlmixr2omega_dOmegaInv(_nlmixr2omega_ind_omega *ome) {
  if (ome->dDomegaInvBool) return ome->dOmegaInvCube;
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

arma::mat nlmixr2omegt_dDomegaInv(_nlmixr2omega_ind_omega *ome) {
  if (ome->dDomegaInvBool) return (ome->dDomegaInvMat);
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
  return (ome->dDomegaInvMat);
}


arma::mat rxToCholOmega(arma::mat cholMat){
  // Only the cholesky is needed for the liklihood calculation
  // trimatu is faster, but it seems to have problems sometimes with certain BLAS combinations:
  // See https://github.com/nlmixrdevelopment/rxode2/issues/84
  // Only the cholesky is needed for the liklihood calculation
  // trimatu is faster, but it seems to have problems sometimes with certain BLAS combinations:
  // See https://github.com/nlmixrdevelopment/rxode2/issues/84
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
  if (ome->cholOmega1Bool) return ome->cholOmega1Mat;
  arma::mat coi =  nlmixr2omega_cholOmegaInv(ome);
  coi = rxToCholOmega(coi);
  ome->cholOmega1Bool = true;
  ome->cholOmega1Mat = coi;
  return ome->cholOmega1Mat;
}

arma::mat nlmixr2omega_omega(_nlmixr2omega_ind_omega *ome) {
  if (ome->omegaBool) return ome->omegaMat;
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
  for (unsigned int i = ts;i--;){
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
  // sum_{k=1}^{n} k = n*(n+1)/2
  arma::vec theta(0.5 * in.n_rows * (in.n_rows + 1.0));
  int j = 0;
  for (unsigned int j = 0; j < in.n_rows; ++j) {
    for (unsigned int i = j; i < in.n_rows; ++i) {
      theta(j)  = in(i, j);
    }
  }
  _nlmixr2omegaAssignTheta(ome, theta);
}

//[[Rcpp::export]]
Rcpp::XPtr<_nlmixr2omega_full_omega> nlmixr2omegaNew(List omeList, int diagXform) {
  _nlmixr2omega_full_omega full;
  _nlmixr2omega_full_omega *fullPtr = &full;
  fullPtr->nomes = omeList.size();
  if (fullPtr->omes != NULL) R_Free(fullPtr->omes);
  fullPtr->omes = R_Calloc(fullPtr->nomes,_nlmixr2omega_ind_omega);
  fullPtr->nTotTheta = 0;
  for (int i = 0; i < fullPtr->nomes; ++i) {
    _nlmixr2omega_ind_omega *ome = &(fullPtr->omes[i]);
    arma::mat cur = as<arma::mat>(omeList[i]);
    nlmixr2omega_iniOmeStruct(ome, cur, diagXform);
    fullPtr->nTotTheta += ome->theta.size();
  }
  Rcpp::XPtr<_nlmixr2omega_full_omega> ptr(fullPtr);
  return ptr;
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

//[[Rcpp::export]]
arma::vec getTheta(Rcpp::XPtr<_nlmixr2omega_full_omega> p) {
  _nlmixr2omega_full_omega* v = p.get();
  return  _nlmixr2omega_full_getTheta_(v);
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
RObject setTheta(Rcpp::XPtr<_nlmixr2omega_full_omega> p, arma::vec theta) {
  _nlmixr2omega_full_omega* v = p.get();
  _nlmixr2omega_full_setTheta_(v, theta);  
  return R_NilValue;
}
