#ifndef __nlmixr2omegaArma_h__
#define __nlmixr2omegaArma_h__

extern "C" {
  typedef void (*_nlmixr2omega_mat_t)(int *dm, double *_t, int *length_theta, int  *_tn, double *ret);
}

struct _nlmixr2omega_ind_omega {
  _nlmixr2omega_mat_t cFun = NULL;
  arma::vec theta;
  int dim;
  
  bool cholOmegaInvBool = false;
  arma::mat cholOmegaInvMat; // "chol.omegaInv"
  
  bool omegaInvBool = false;
  arma::mat omegaInvMat; // "omegaInv"
  
  bool dOmegaInvBool = false;
  arma::cube dOmegaInvCube; // dim x dim x ntheta "d.omegaInv"
  
  bool dDomegaInvBool = false;
  arma::mat dDomegaInvMat; // ntheta x dim; "d.D.omegaInv"
  
  bool cholOmega1Bool = false;
  arma::mat cholOmega1Mat; // "chol.omega1"
  
  bool omegaBool = false;
  arma::mat omegaMat; // "omega"
  
  bool cholOmegaBool = false;
  arma::mat cholOmegaMat; //"chol.omega"
  
  bool logDetOMGAinv5Bool = false;
  double logDetOMGAinv5double; //"log.det.OMGAinv.5"
  
  bool tr28Bool = false;
  arma::vec tr28vec; // "tr.28"
  
  // could be neta x neta x ntheta or dim x dim x ntheta
  bool omega47Bool = false;
  arma::cube omega47Cube; 
};

class _nlmixr2omega_full_omega {
public:
  _nlmixr2omega_full_omega() {
    
  }
  ~_nlmixr2omega_full_omega() {
    if (omes != NULL) R_Free(omes);
  }
  _nlmixr2omega_ind_omega *omes = NULL;
  int nomes = 0;
  int nTotTheta = 0;
  int nTotDim = 0;
};


#define nlmixr2omega_sqrt 1
#define nlmixr2omega_log 2

#endif
