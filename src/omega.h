#ifndef __omega_h__
#define __omega_h__
#if defined(__cplusplus)
extern "C" {
#endif

  int _nlmixr2omega_matSize(void);
  void _nlmixr2omega_mat(int *dm, double *_t, int *length_theta, int  *_tn, double *ret);
  
#if defined(__cplusplus)
}
#endif

#endif
