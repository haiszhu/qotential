// gcc -mavx512f -mfma -c -o ckernels_expts.o ckernels_expts.c
#ifdef __cplusplus
extern "C" {
    void clap3ddlpmat_c_( int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *A);
    void csimd128lap3ddlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *A);
    void csimd256lap3ddlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *A);
    void csimd512lap3ddlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *A);
    void clap3dslpmat_c_( int64_t *M, double *r0, int64_t *N, double *r, double *w, double *A);
    void csimd256lap3dslpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *w, double *A);
    void csimd512lap3dslpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *w, double *A);
    void csimd256lap3dsdlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *As, double *Ad);
    void csimd512lap3dsdlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *As, double *Ad);
    void cline3omega0123vr_c_(double *r0, double *m_all, int64_t *nbd, int64_t *order, double *onm0, double *onm1, double *onm2, double *onm3, int64_t *h_dim, double *omega0, double *omega1, double *omega2, double *omega3);
    void csimd128line3omega0123vr_c_(double *r0, double *m_all, int64_t *nbd, int64_t *order, double *onm0, double *onm1, double *onm2, double *onm3, int64_t *h_dim, double *omega0, double *omega1, double *omega2, double *omega3);
    void csimd256line3omega0123vr_c_(double *r0, double *m_all, int64_t *nbd, int64_t *order, double *onm0, double *onm1, double *onm2, double *onm3, int64_t *h_dim, double *omega0, double *omega1, double *omega2, double *omega3);
    void csimd512line3omega0123vr_c_(double *r0, double *m_all, int64_t *nbd, int64_t *order, double *onm0, double *onm1, double *onm2, double *onm3, int64_t *h_dim, double *omega0, double *omega1, double *omega2, double *omega3);
    void cmomentsad_vr_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk);
    void csimd128momentsad_vr_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk);
    void csimd256momentsad_vr_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk);
    void csimd512momentsad_vr_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk);
    // version2 of c simd momentsad_vr...
    // we do not make use of the fast log kernel evaluation since that depends on intel 
    // this one is ok to use, avx512 only gives 50% speedup when n = n4 or n8 in standalone test but not so much (maybe 5-10%) in the line3adpmoments subroutine... not sure why?
    // look at assist_simd_momentsad_vr, returns Nk, Mk
    // maybe should do another attemp when understand how to vectorize efficiently better... 
    void csimd256momentsad_vrv2_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk);
    void csimd512momentsad_vrv2_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk);
    }
#endif

// detect platform
#if defined(__arm__) || defined(__aarch64__)
  #define ARCH_ARM
#elif defined(__x86_64__) || defined(_M_X64)
  #if defined(__AVX2__)
    #define ARCH_AVX2
  #endif
  #define ARCH_X86
  #if defined(__AVX512F__)
    #define ARCH_AVX512
  #endif
#endif

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>


/*=====================================================================================================*/
// general C
void clap3ddlpmat_c_(  int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *A )
{
  int64_t m,n,i,j;
  double dx, dy, dz, rr;
  double pi4inv = 1.0/(16.0*atan(1.0));
  double threshsq = 1.0e-28;
  m = *M;
  n = *N;
  // 
  double *r0x = (double*)aligned_alloc(8, m * sizeof(double));
  double *r0y = (double*)aligned_alloc(8, m * sizeof(double));
  double *r0z = (double*)aligned_alloc(8, m * sizeof(double));
  // 
  for (i = 0; i < m; ++i) {
    r0x[i] = r0[3*i];
    r0y[i] = r0[3*i+1];
    r0z[i] = r0[3*i+2];
  }
  // 
  for (j=0; j<n; j++){
    for (i=0; i<m; i++){
      dx = r0x[i] - r[3*j];
      dy = r0y[i] - r[3*j+1];
      dz = r0z[i] - r[3*j+2];
      rr = dx*dx + dy*dy + dz*dz;
      if (rr > threshsq) {
        A[m*j+i] = pi4inv*w[j]*(dx*rn[3*j]+dy*rn[3*j+1]+dz*rn[3*j+2])/(sqrt(rr)*rr);
      } else {
        A[m*j+i] = 0.0;
      }
    }
  }
  // Free allocated memory
  free(r0x);
  free(r0y);
  free(r0z);
}
void clap3dslpmat_c_(  int64_t *M, double *r0, int64_t *N, double *r, double *w, double *A )
{
  int64_t m,n,i,j;
  double dx, dy, dz, rr;
  double pi4inv = 1.0/(16.0*atan(1.0));
  double threshsq = 1.0e-28;
  m = *M;
  n = *N;
  // 
  double *r0x = (double*)aligned_alloc(8, m * sizeof(double));
  double *r0y = (double*)aligned_alloc(8, m * sizeof(double));
  double *r0z = (double*)aligned_alloc(8, m * sizeof(double));
  // 
  for (i = 0; i < m; ++i) {
    r0x[i] = r0[3*i];
    r0y[i] = r0[3*i+1];
    r0z[i] = r0[3*i+2];
  }
  // 
  for (j=0; j<n; j++){
    for (i=0; i<m; i++){
      dx = r0x[i] - r[3*j];
      dy = r0y[i] - r[3*j+1];
      dz = r0z[i] - r[3*j+2];
      rr = dx*dx + dy*dy + dz*dz;
      if (rr > threshsq) {
        A[m*j+i] = pi4inv*w[j]/sqrt(rr);
      } else {
        A[m*j+i] = 0.0;
      }
    }
  }
  // Free allocated memory
  free(r0x);
  free(r0y);
  free(r0z);
}
void cline3omega0123vr_c_(double *r0, double *m_all, int64_t *Nbd, int64_t *Order, double *onm0, double *onm1, double *onm2, double *onm3, int64_t *H_dim, double *omega0, double *omega1, double *omega2, double *omega3)
{
  int64_t nbd, order, h_dim, i, ii, j, l, nbd4;
  double pi, tmpval, r00, r01, r02;
  nbd = *Nbd;
  order = *Order;
  h_dim = *H_dim;
  nbd4 = 4*nbd;
  // 
  double *mk_j_ii = (double*)aligned_alloc(8, 4 * nbd * sizeof(double));
  pi = 4.0*atan(1.0);
  tmpval = -1.0/(4.0*pi);
  r00 = tmpval*r0[0];
  r01 = tmpval*r0[1];
  r02 = tmpval*r0[2];
  // 
  for (i = 0; i < h_dim; ++i){
    omega0[i] = 0.0;
    omega1[i] = 0.0;
    omega2[i] = 0.0;
    omega3[i] = 0.0;
  }
  // loop over moment order
  i = 0;
  for (ii = 0; ii < order+1; ++ii){
    // 
    for (j = 0; j < nbd; ++j){
      mk_j_ii[4*j]   = tmpval*m_all[(ii+1)*nbd+j];
      mk_j_ii[4*j+1] = r00*m_all[ii*nbd+j];
      mk_j_ii[4*j+2] = r01*m_all[ii*nbd+j];
      mk_j_ii[4*j+3] = r02*m_all[ii*nbd+j];
    }
    // 
    for (l = 0; l < ii; ++l){
      for (j = 0; j < nbd4; ++j){
        omega0[i] +=  mk_j_ii[j] * onm0[nbd4*i+j]; 
        omega1[i] +=  mk_j_ii[j] * onm1[nbd4*i+j]; 
        omega2[i] +=  mk_j_ii[j] * onm2[nbd4*i+j]; 
        omega3[i] +=  mk_j_ii[j] * onm3[nbd4*i+j]; 
      }
      ++i;
    }
  }
}
void cmomentsad_vr_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk)
{
  int64_t n, order, k, j;
  n = *N;
  order = *Order;
  double r0norm, r0norm_inv, tmp, rx, ry, rz, rmr0x, rmr0y, rmr0z;
  double LMNcommon, rnorm2, rnorm, rnorm_inv, r0normplusrnorm, r0mr;
  double *r0dotr = (double*)aligned_alloc(8, n * sizeof(double));
  double *rnorm2_inv = (double*)aligned_alloc(8, n * sizeof(double));
  double *r0dotr_over_rnorm2 = (double*)aligned_alloc(8, n * sizeof(double));
  double *r0norm2_over_rnorm2 = (double*)aligned_alloc(8, n * sizeof(double));
  double *r0mr_inv = (double*)aligned_alloc(8, n * sizeof(double));
  double *r0mr_over_rnorm2 = (double*)aligned_alloc(8, n * sizeof(double));
  // 
  r0norm = sqrt(r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2]);
  r0norm_inv = 1.0/r0norm;
  // 
  for (k = 0; k < n; ++k){
    rx = r[3*k];
    ry = r[3*k+1];
    rz = r[3*k+2];
    rmr0x = r0[0] - rx; 
    rmr0y = r0[1] - ry; 
    rmr0z = r0[2] - rz;
    rnorm2 = rx*rx + ry*ry + rz*rz;
    rnorm = sqrt(rnorm2);
    rnorm_inv = 1.0/rnorm;
    r0normplusrnorm = r0norm + rnorm;
    r0mr = sqrt(rmr0x*rmr0x + rmr0y*rmr0y + rmr0z*rmr0z);
    // to be used later
    r0dotr[k] = r0[0]*rx + r0[1]*ry + r0[2]*rz;
    rnorm2_inv[k] = rnorm_inv*rnorm_inv;
    r0dotr_over_rnorm2[k] = r0dotr[k]*rnorm2_inv[k];
    r0norm2_over_rnorm2[k] = (r0norm*r0norm)*rnorm2_inv[k];
    r0mr_inv[k] = 1.0/r0mr;
    r0mr_over_rnorm2[k] = r0mr*rnorm2_inv[k];
    // 
    LMNcommon = r0normplusrnorm/(r0normplusrnorm
                  + r0mr_inv[k]*( r0norm*r0norm - rnorm2));
    tmp = (r0normplusrnorm+r0mr)*LMNcommon*r0mr_inv[k];
    Nk[k] = log(tmp)*rnorm_inv;
    Nk[n+k] = Nk[k]*r0dotr_over_rnorm2[k] 
              + (r0mr-r0norm)*rnorm2_inv[k];
    // 
    LMNcommon = LMNcommon/(r0norm*rnorm+r0dotr[k]);
    Mk[k] = LMNcommon*((r0normplusrnorm* 
        r0mr_inv[k]+rnorm*r0norm_inv)*r0mr_inv[k]-r0norm_inv);
    Mk[n+k] = r0norm*Mk[k]/(r0norm+r0mr);
  }
  for (k = 1; k < order; ++k){
    for (j = 0; j < n; ++j){
      Nk[(k+1)*n+j] = ((2*k+1)*r0dotr_over_rnorm2[j]*Nk[k*n+j] -
                        k*r0norm2_over_rnorm2[j]*Nk[(k-1)*n+j] +
                      r0mr_over_rnorm2[j])/(k+1);
      Mk[(k+1)*n+j] = (r0dotr[j]*Mk[k*n+j]+k*Nk[(k-1)*n+j]-r0mr_inv[j])*rnorm2_inv[j];
    }
  }
}



#ifdef ARCH_X86
#include <immintrin.h>


/*=====================================================================================================*/
// is there such a thing for X86?
void csimd128lap3ddlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *A) {}
void csimd128line3omega0123vr_c_(double *r0, double *m_all, int64_t *Nbd, int64_t *Order, double *onm0, double *onm1, double *onm2, double *onm3, int64_t *H_dim, double *omega0, double *omega1, double *omega2, double *omega3) {}
void csimd128momentsad_vr_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk) {}

/*=====================================================================================================*/
#if defined(ARCH_AVX2)
void csimd256lap3ddlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *A) 
{
  int64_t m, n, i, j, m4;
  double pi4inv = 1.0/(16.0*atan(1.0));
  m = *M;
  n = *N;
  m4 = (m/4)*4;
  // 
  double *r0x = (double*)aligned_alloc(32, m4 * sizeof(double));
  double *r0y = (double*)aligned_alloc(32, m4 * sizeof(double));
  double *r0z = (double*)aligned_alloc(32, m4 * sizeof(double));
  // 
  for (i = 0; i < m4; ++i) {
    r0x[i] = r0[3*i];
    r0y[i] = r0[3*i+1];
    r0z[i] = r0[3*i+2];
  }
  //   
  __m256d half = _mm256_set1_pd(0.5);
  __m256d three = _mm256_set1_pd(3.0); // for newton iteration
  for (j = 0; j < n; ++j) {
    __m256d rx_v  = _mm256_set1_pd(r[3*j]);
    __m256d ry_v  = _mm256_set1_pd(r[3*j+1]);
    __m256d rz_v  = _mm256_set1_pd(r[3*j+2]);
    __m256d rnx_v = _mm256_set1_pd(rn[3*j]);
    __m256d rny_v = _mm256_set1_pd(rn[3*j+1]);
    __m256d rnz_v = _mm256_set1_pd(rn[3*j+2]);
    __m256d w_v   = _mm256_set1_pd(pi4inv*w[j]);
    for (i = 0; i <= m4 - 4; i += 4) {
      __m256d r0x_v    = _mm256_load_pd(&r0x[i]);
      __m256d r0y_v    = _mm256_load_pd(&r0y[i]);
      __m256d r0z_v    = _mm256_load_pd(&r0z[i]);
      __m256d dx_v     = _mm256_sub_pd(r0x_v, rx_v);
      __m256d dy_v     = _mm256_sub_pd(r0y_v, ry_v);
      __m256d dz_v     = _mm256_sub_pd(r0z_v, rz_v);
      __m256d rr_v     = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(dx_v, dx_v), _mm256_mul_pd(dy_v, dy_v)), _mm256_mul_pd(dz_v, dz_v));
      __m256d rdotrn_v = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(dx_v, rnx_v), _mm256_mul_pd(dy_v, rny_v)), _mm256_mul_pd(dz_v, rnz_v));
      // copied from wen yan's DemoSIMD: https://github.com/wenyan4work/DemoSIMD/blob/master/src/rsqrt.cpp
      __m256d ri_v = _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(rr_v))); // inverse sqrt approximate
      __m256d muls = _mm256_mul_pd(_mm256_mul_pd(rr_v, ri_v), ri_v);
      ri_v = _mm256_mul_pd(_mm256_mul_pd(half, ri_v), _mm256_sub_pd(three, muls)); // newton iteration
      muls = _mm256_mul_pd(_mm256_mul_pd(rr_v, ri_v), ri_v);
      ri_v = _mm256_mul_pd(_mm256_mul_pd(half, ri_v), _mm256_sub_pd(three, muls)); // newton iteration
      rr_v = _mm256_mul_pd(_mm256_mul_pd(ri_v,ri_v),ri_v); // ri_v**3 misuse of variable name
      __m256d dlpij_v  = _mm256_mul_pd(_mm256_mul_pd(w_v,rdotrn_v), rr_v);
      // __m256d dlpij_v  = _mm256_div_pd(_mm256_mul_pd(w_v,rdotrn_v), _mm256_mul_pd(_mm256_sqrt_pd(rr_v),rr_v));
      // _mm256_store_pd(&A[i+j*m], inv_sqrt_rr_v); // this one crashes
      _mm256_storeu_pd(&A[i+j*m], dlpij_v); // Moves to unaligned memory location

    }
    // compute the rest, simd over source? or not worth it
    for (i=m4; i<m; i++){
      double dx = r0[3*i]   - r[3*j];
      double dy = r0[3*i+1] - r[3*j+1];
      double dz = r0[3*i+2] - r[3*j+2];
      double rr = dx*dx + dy*dy + dz*dz;
      A[m*j+i] = pi4inv*w[j]*(dx*rn[3*j]+dy*rn[3*j+1]+dz*rn[3*j+2])/(sqrt(rr)*rr);
    }
  }
  free(r0x);
  free(r0y);
  free(r0z);
}
void csimd256lap3dslpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *w, double *A) 
{
  int64_t m, n, i, j, m4;
  double pi4inv = 1.0/(16.0*atan(1.0));
  m = *M;
  n = *N;
  m4 = (m/4)*4;
  // 
  double *r0x = (double*)aligned_alloc(32, m4 * sizeof(double));
  double *r0y = (double*)aligned_alloc(32, m4 * sizeof(double));
  double *r0z = (double*)aligned_alloc(32, m4 * sizeof(double));
  // 
  for (i = 0; i < m4; ++i) {
    r0x[i] = r0[3*i];
    r0y[i] = r0[3*i+1];
    r0z[i] = r0[3*i+2];
  }
  //   
  __m256d half = _mm256_set1_pd(0.5);
  __m256d three = _mm256_set1_pd(3.0); // for newton iteration
  for (j = 0; j < n; ++j) {
    __m256d rx_v  = _mm256_set1_pd(r[3*j]);
    __m256d ry_v  = _mm256_set1_pd(r[3*j+1]);
    __m256d rz_v  = _mm256_set1_pd(r[3*j+2]);
    __m256d w_v   = _mm256_set1_pd(pi4inv*w[j]);
    for (i = 0; i <= m4 - 4; i += 4) {
      __m256d r0x_v    = _mm256_load_pd(&r0x[i]);
      __m256d r0y_v    = _mm256_load_pd(&r0y[i]);
      __m256d r0z_v    = _mm256_load_pd(&r0z[i]);
      __m256d dx_v     = _mm256_sub_pd(r0x_v, rx_v);
      __m256d dy_v     = _mm256_sub_pd(r0y_v, ry_v);
      __m256d dz_v     = _mm256_sub_pd(r0z_v, rz_v);
      __m256d rr_v     = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(dx_v, dx_v), _mm256_mul_pd(dy_v, dy_v)), _mm256_mul_pd(dz_v, dz_v));
      // copied from wen yan's DemoSIMD: https://github.com/wenyan4work/DemoSIMD/blob/master/src/rsqrt.cpp
      __m256d ri_v = _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(rr_v))); // inverse sqrt approximate
      __m256d muls = _mm256_mul_pd(_mm256_mul_pd(rr_v, ri_v), ri_v);
      ri_v = _mm256_mul_pd(_mm256_mul_pd(half, ri_v), _mm256_sub_pd(three, muls)); // newton iteration
      muls = _mm256_mul_pd(_mm256_mul_pd(rr_v, ri_v), ri_v);
      ri_v = _mm256_mul_pd(_mm256_mul_pd(half, ri_v), _mm256_sub_pd(three, muls)); // newton iteration
      __m256d slpij_v  = _mm256_mul_pd(w_v, ri_v);
      // __m256d dlpij_v  = _mm256_div_pd(_mm256_mul_pd(w_v,rdotrn_v), _mm256_mul_pd(_mm256_sqrt_pd(rr_v),rr_v));
      // _mm256_store_pd(&A[i+j*m], inv_sqrt_rr_v); // this one crashes
      _mm256_storeu_pd(&A[i+j*m], slpij_v); // Moves to unaligned memory location

    }
    // compute the rest, simd over source? or not worth it
    for (i=m4; i<m; i++){
      double dx = r0[3*i]   - r[3*j];
      double dy = r0[3*i+1] - r[3*j+1];
      double dz = r0[3*i+2] - r[3*j+2];
      double rr = dx*dx + dy*dy + dz*dz;
      A[m*j+i] = pi4inv*w[j]/sqrt(rr);
    }
  }
  free(r0x);
  free(r0y);
  free(r0z);
}
void csimd256lap3dsdlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *As, double *Ad) 
{
  int64_t m, n, i, j, m4;
  double pi4inv = 1.0/(16.0*atan(1.0));
  m = *M;
  n = *N;
  m4 = (m/4)*4;
  // 
  double *r0x = (double*)aligned_alloc(32, m4 * sizeof(double));
  double *r0y = (double*)aligned_alloc(32, m4 * sizeof(double));
  double *r0z = (double*)aligned_alloc(32, m4 * sizeof(double));
  // 
  for (i = 0; i < m4; ++i) {
    r0x[i] = r0[3*i];
    r0y[i] = r0[3*i+1];
    r0z[i] = r0[3*i+2];
  }
  //   
  __m256d half = _mm256_set1_pd(0.5);
  __m256d three = _mm256_set1_pd(3.0); // for newton iteration
  for (j = 0; j < n; ++j) {
    __m256d rx_v  = _mm256_set1_pd(r[3*j]);
    __m256d ry_v  = _mm256_set1_pd(r[3*j+1]);
    __m256d rz_v  = _mm256_set1_pd(r[3*j+2]);
    __m256d rnx_v = _mm256_set1_pd(rn[3*j]);
    __m256d rny_v = _mm256_set1_pd(rn[3*j+1]);
    __m256d rnz_v = _mm256_set1_pd(rn[3*j+2]);
    __m256d w_v   = _mm256_set1_pd(pi4inv*w[j]);
    for (i = 0; i <= m4 - 4; i += 4) {
      __m256d r0x_v    = _mm256_load_pd(&r0x[i]);
      __m256d r0y_v    = _mm256_load_pd(&r0y[i]);
      __m256d r0z_v    = _mm256_load_pd(&r0z[i]);
      __m256d dx_v     = _mm256_sub_pd(r0x_v, rx_v);
      __m256d dy_v     = _mm256_sub_pd(r0y_v, ry_v);
      __m256d dz_v     = _mm256_sub_pd(r0z_v, rz_v);
      __m256d rr_v     = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(dx_v, dx_v), _mm256_mul_pd(dy_v, dy_v)), _mm256_mul_pd(dz_v, dz_v));
      __m256d rdotrn_v = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(dx_v, rnx_v), _mm256_mul_pd(dy_v, rny_v)), _mm256_mul_pd(dz_v, rnz_v));
      // copied from wen yan's DemoSIMD: https://github.com/wenyan4work/DemoSIMD/blob/master/src/rsqrt.cpp
      __m256d ri_v = _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(rr_v))); // inverse sqrt approximate
      __m256d muls = _mm256_mul_pd(_mm256_mul_pd(rr_v, ri_v), ri_v);
      ri_v = _mm256_mul_pd(_mm256_mul_pd(half, ri_v), _mm256_sub_pd(three, muls)); // newton iteration
      muls = _mm256_mul_pd(_mm256_mul_pd(rr_v, ri_v), ri_v);
      ri_v = _mm256_mul_pd(_mm256_mul_pd(half, ri_v), _mm256_sub_pd(three, muls)); // newton iteration
      rr_v = _mm256_mul_pd(_mm256_mul_pd(ri_v,ri_v),ri_v); // ri_v**3 misuse of variable name
      __m256d dlpij_v  = _mm256_mul_pd(_mm256_mul_pd(w_v,rdotrn_v), rr_v);
      __m256d slpij_v  = _mm256_mul_pd(w_v, ri_v);
      // __m256d dlpij_v  = _mm256_div_pd(_mm256_mul_pd(w_v,rdotrn_v), _mm256_mul_pd(_mm256_sqrt_pd(rr_v),rr_v));
      // _mm256_store_pd(&A[i+j*m], inv_sqrt_rr_v); // this one crashes
      _mm256_storeu_pd(&Ad[i+j*m], dlpij_v); // Moves to unaligned memory location
      _mm256_storeu_pd(&As[i+j*m], slpij_v); 
    }
    // compute the rest, simd over source? or not worth it
    for (i=m4; i<m; i++){
      double dx = r0[3*i]   - r[3*j];
      double dy = r0[3*i+1] - r[3*j+1];
      double dz = r0[3*i+2] - r[3*j+2];
      double rr = dx*dx + dy*dy + dz*dz;
      Ad[m*j+i] = pi4inv*w[j]*(dx*rn[3*j]+dy*rn[3*j+1]+dz*rn[3*j+2])/(sqrt(rr)*rr);
      As[m*j+i] = pi4inv*w[j]/sqrt(rr);
    }
  }
  free(r0x);
  free(r0y);
  free(r0z);
}
void csimd256line3omega0123vr_c_(double *r0, double *m_all, int64_t *Nbd, int64_t *Order, double *onm0, double *onm1, double *onm2, double *onm3, int64_t *H_dim, double *omega0, double *omega1, double *omega2, double *omega3) 
{
  int64_t nbd, order, h_dim, i, ii, j, l, nbd4, nbd44;
  double pi, tmpval, r00, r01, r02;
  nbd = *Nbd;
  order = *Order;
  h_dim = *H_dim;
  nbd4 = 4*nbd;
  nbd44 = (nbd4/4)*4;
  // 
  double *mk_j_ii = (double*)aligned_alloc(32, 4 * nbd * sizeof(double));
  pi = 4.0*atan(1.0);
  tmpval = -1.0/(4.0*pi);
  r00 = tmpval*r0[0];
  r01 = tmpval*r0[1];
  r02 = tmpval*r0[2];
  // 
  for (i = 0; i < h_dim; ++i){
    omega0[i] = 0.0;
    omega1[i] = 0.0;
    omega2[i] = 0.0;
    omega3[i] = 0.0;
  }
  // loop over moment order
  i = 0;
  for (ii = 0; ii < order+1; ++ii){
    // 
    for (j = 0; j < nbd; ++j){
      mk_j_ii[4*j]   = tmpval*m_all[(ii+1)*nbd+j];
      mk_j_ii[4*j+1] = r00*m_all[ii*nbd+j];
      mk_j_ii[4*j+2] = r01*m_all[ii*nbd+j];
      mk_j_ii[4*j+3] = r02*m_all[ii*nbd+j];
    }
    // 
    for (l = 0; l < ii; ++l){
      __m256d omega0_v = _mm256_set1_pd(0.0);
      __m256d omega1_v = _mm256_set1_pd(0.0);
      __m256d omega2_v = _mm256_set1_pd(0.0);
      __m256d omega3_v = _mm256_set1_pd(0.0);
      for (j = 0; j <= nbd44 - 4; j += 4){
        __m256d mk = _mm256_load_pd(&mk_j_ii[j]);
        __m256d onm0_v = _mm256_load_pd(&onm0[nbd4*i+j]);
        __m256d onm1_v = _mm256_load_pd(&onm1[nbd4*i+j]);
        __m256d onm2_v = _mm256_load_pd(&onm2[nbd4*i+j]);
        __m256d onm3_v = _mm256_load_pd(&onm3[nbd4*i+j]);
        omega0_v = _mm256_fmadd_pd(mk, onm0_v, omega0_v);
        omega1_v = _mm256_fmadd_pd(mk, onm1_v, omega1_v);
        omega2_v = _mm256_fmadd_pd(mk, onm2_v, omega2_v);
        omega3_v = _mm256_fmadd_pd(mk, onm3_v, omega3_v);
      }
      // omega0[i] = _mm256_reduce_add_pd(omega0_v); // no such thing...
      __m256d temp = _mm256_hadd_pd(omega0_v, omega0_v);
      omega0[i] = ((double*)&temp)[0] + ((double*)&temp)[2];
      temp = _mm256_hadd_pd(omega1_v, omega1_v);
      omega1[i] = ((double*)&temp)[0] + ((double*)&temp)[2];
      temp = _mm256_hadd_pd(omega2_v, omega2_v);
      omega2[i] = ((double*)&temp)[0] + ((double*)&temp)[2];
      temp = _mm256_hadd_pd(omega3_v, omega3_v);
      omega3[i] = ((double*)&temp)[0] + ((double*)&temp)[2];
      // compute the rest if needed
      for (j=nbd44; j<nbd4; j++){
        omega0[i] +=  mk_j_ii[j] * onm0[nbd4*i+j]; 
        omega1[i] +=  mk_j_ii[j] * onm1[nbd4*i+j]; 
        omega2[i] +=  mk_j_ii[j] * onm2[nbd4*i+j]; 
        omega3[i] +=  mk_j_ii[j] * onm3[nbd4*i+j]; 
      }
      ++i;
    }
  }
  free(mk_j_ii);
}
void csimd256momentsad_vr_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk)
{
// to be implemented, look at csimd512momentsad_vr_c_
}
// version 2
void csimd256momentsad_vrv2_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk)
{
  int64_t n = *N, order = *Order, k, j, n4, nmn4, idx;
  double tmpk1, tmpk2;
  double r0norm2, r0norm, r0norm_inv;
  n4 = (n/4)*4;
  double *r4x = (double*)aligned_alloc(32, n4 * sizeof(double));
  double *r4y = (double*)aligned_alloc(32, n4 * sizeof(double));
  double *r4z = (double*)aligned_alloc(32, n4 * sizeof(double));
  double *tmpk3 = (double*)aligned_alloc(32, order * sizeof(double));
  r0norm2 = (r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2]);
  r0norm  = sqrt(r0norm2);
  r0norm_inv = 1.0/r0norm;
  __m256d r0x_v = _mm256_set1_pd(r0[0]);
  __m256d r0y_v = _mm256_set1_pd(r0[1]);
  __m256d r0z_v = _mm256_set1_pd(r0[2]);
  __m256d r0norm_v = _mm256_set1_pd(r0norm);
  __m256d r0norm2_v = _mm256_set1_pd(r0norm2);
  __m256d r0norm_inv_v = _mm256_set1_pd(r0norm_inv);
  __m256d one_v = _mm256_set1_pd(1.0);
  for (j = 0; j < n4; ++j) {
    r4x[j] = r[3*j];
    r4y[j] = r[3*j+1];
    r4z[j] = r[3*j+2];
  }
  for (k = 1; k < order; ++k) {
    tmpk3[k] = 1.0/(k+1);
  }
  idx = 0;
  __m256d Nnext_v = _mm256_set1_pd(0.0);
  __m256d Mnext_v = _mm256_set1_pd(0.0);
  for (j = 0; j < n4; j += 4){
    __m256d rx_v = _mm256_load_pd(&r4x[j]);
    __m256d ry_v = _mm256_load_pd(&r4y[j]);
    __m256d rz_v = _mm256_load_pd(&r4z[j]);
    __m256d rmr0x_v = _mm256_sub_pd(r0x_v, rx_v);
    __m256d rmr0y_v = _mm256_sub_pd(r0y_v, ry_v);
    __m256d rmr0z_v = _mm256_sub_pd(r0z_v, rz_v);
    // 
    __m256d rnorm2_v = _mm256_fmadd_pd(rx_v, rx_v, _mm256_fmadd_pd(ry_v, ry_v, _mm256_mul_pd(rz_v, rz_v)));
    __m256d rnorm_v = _mm256_sqrt_pd(rnorm2_v);
    __m256d rnorm_inv_v = _mm256_div_pd(one_v, rnorm_v);
    __m256d rnorm2_inv_v = _mm256_mul_pd(rnorm_inv_v, rnorm_inv_v);
    __m256d r0mr_v = _mm256_sqrt_pd(_mm256_fmadd_pd(rmr0x_v, rmr0x_v, _mm256_fmadd_pd(rmr0y_v, rmr0y_v, _mm256_mul_pd(rmr0z_v, rmr0z_v))));
    __m256d r0mr_inv_v = _mm256_div_pd(one_v, r0mr_v);
    __m256d r0dotr_v = _mm256_fmadd_pd(r0x_v, rx_v, _mm256_fmadd_pd(r0y_v, ry_v, _mm256_mul_pd(r0z_v, rz_v)));
    // 
    __m256d rnorm2mr0dotr_v = _mm256_sub_pd(rnorm2_v, r0dotr_v);
    __m256d r0rmr0dotr_inv_v = _mm256_div_pd(one_v, _mm256_fmsub_pd(r0norm_v, rnorm_v, r0dotr_v));
    __m256d tmp_v = _mm256_mul_pd(_mm256_fmadd_pd(rnorm_v, r0mr_v, rnorm2mr0dotr_v), r0rmr0dotr_inv_v);
    __m256d Nprev_v = {log(tmp_v[0]), log(tmp_v[1]), log(tmp_v[2]), log(tmp_v[3])};
    Nprev_v = _mm256_mul_pd(rnorm_inv_v, Nprev_v);
    __m256d Ncurr_v = _mm256_mul_pd(rnorm2_inv_v, _mm256_fmadd_pd(r0dotr_v, Nprev_v, _mm256_sub_pd(r0mr_v, r0norm_v)));
    // 
    _mm256_storeu_pd(&Nk[j], Nprev_v);
    _mm256_storeu_pd(&Nk[n + j], Ncurr_v);
    // 
    __m256d r0rpr0dotr_inv_v = _mm256_div_pd(one_v,_mm256_fmadd_pd(r0norm_v, rnorm_v, r0dotr_v));
    __m256d term1_v = _mm256_mul_pd(r0rpr0dotr_inv_v, r0rmr0dotr_inv_v);
    __m256d Mprev_v = _mm256_mul_pd(term1_v,_mm256_fmadd_pd(rnorm2mr0dotr_v, r0mr_inv_v, _mm256_mul_pd(r0dotr_v,r0norm_inv_v)));
    __m256d Mcurr_v = _mm256_mul_pd(term1_v, _mm256_fmadd_pd(_mm256_sub_pd(r0dotr_v, r0norm2_v), r0mr_inv_v, r0norm_v));
    // 
    _mm256_storeu_pd(&Mk[j], Mprev_v);
    _mm256_storeu_pd(&Mk[n + j], Mcurr_v);
    for (k = 1; k < order; ++k) {
      Nnext_v = _mm256_mul_pd(_mm256_fmadd_pd( _mm256_set1_pd(2*k+1), _mm256_mul_pd(r0dotr_v, Ncurr_v), 
                              _mm256_fmadd_pd( _mm256_set1_pd(-k), _mm256_mul_pd(r0norm2_v, Nprev_v), r0mr_v)),
                              _mm256_mul_pd(rnorm2_inv_v, _mm256_set1_pd(tmpk3[k])));
      Mnext_v = _mm256_mul_pd(_mm256_fmadd_pd(r0dotr_v, Mcurr_v, 
                              _mm256_fmsub_pd(_mm256_set1_pd(k), Nprev_v, r0mr_inv_v)), 
                                rnorm2_inv_v);
      _mm256_storeu_pd(&Nk[(k+1)*n + j], Nnext_v);
      _mm256_storeu_pd(&Mk[(k+1)*n + j], Mnext_v);
      Nprev_v = Ncurr_v;
      Ncurr_v = Nnext_v;
      Mcurr_v = Mnext_v;
    }
  }
  nmn4 = n - n4;
  if (nmn4 > 0) {
    double tmp, rx, ry, rz, rmr0x, rmr0y, rmr0z;
    double LMNcommon, rnorm2, rnorm, rnorm_inv, r0normplusrnorm, r0mr;
    double *r0dotr = (double*)aligned_alloc(32, nmn4 * sizeof(double));
    double *rnorm2_inv = (double*)aligned_alloc(32, nmn4 * sizeof(double));
    double *r0dotr_over_rnorm2 = (double*)aligned_alloc(32, nmn4 * sizeof(double));
    double *r0norm2_over_rnorm2 = (double*)aligned_alloc(32, nmn4 * sizeof(double));
    double *r0mr_inv = (double*)aligned_alloc(32, nmn4 * sizeof(double));
    double *r0mr_over_rnorm2 = (double*)aligned_alloc(32, nmn4 * sizeof(double));
    for (k=n4; k<n; k++){
      rx = r[3*k];
      ry = r[3*k+1];
      rz = r[3*k+2];
      rmr0x = r0[0] - rx; 
      rmr0y = r0[1] - ry; 
      rmr0z = r0[2] - rz;
      rnorm2 = rx*rx + ry*ry + rz*rz;
      rnorm = sqrt(rnorm2);
      rnorm_inv = 1.0/rnorm;
      r0normplusrnorm = r0norm + rnorm;
      r0mr = sqrt(rmr0x*rmr0x + rmr0y*rmr0y + rmr0z*rmr0z);
      // to be used later
      r0dotr[k-n4] = r0[0]*rx + r0[1]*ry + r0[2]*rz;
      rnorm2_inv[k-n4] = rnorm_inv*rnorm_inv;
      r0dotr_over_rnorm2[k-n4] = r0dotr[k-n4]*rnorm2_inv[k-n4];
      r0norm2_over_rnorm2[k-n4] = (r0norm*r0norm)*rnorm2_inv[k-n4];
      r0mr_inv[k-n4] = 1.0/r0mr;
      r0mr_over_rnorm2[k-n4] = r0mr*rnorm2_inv[k-n4];
      // 
      LMNcommon = r0normplusrnorm/(r0normplusrnorm
                    + r0mr_inv[k-n4]*( r0norm*r0norm - rnorm2));
      tmp = (r0normplusrnorm+r0mr)*LMNcommon*r0mr_inv[k-n4];
      Nk[k] = log(tmp)*rnorm_inv;
      Nk[n+k] = Nk[k]*r0dotr_over_rnorm2[k-n4] 
                + (r0mr-r0norm)*rnorm2_inv[k-n4];
      // 
      LMNcommon = LMNcommon/(r0norm*rnorm+r0dotr[k-n4]);
      Mk[k] = LMNcommon*((r0normplusrnorm* 
        r0mr_inv[k-n4]+rnorm*r0norm_inv)*r0mr_inv[k-n4]-r0norm_inv);
      Mk[n+k] = r0norm*Mk[k]/(r0norm+r0mr);
    }
    for (k = 1; k < order; ++k){
      tmpk1 = (2*k+1)*tmpk3[k];
      tmpk2 = k*tmpk3[k];
      for (j = n4; j < n; ++j){
        Nk[(k+1)*n+j] = tmpk1*r0dotr_over_rnorm2[j-n4]*Nk[k*n+j] -
                        tmpk2*r0norm2_over_rnorm2[j-n4]*Nk[(k-1)*n+j] +
                        tmpk3[k]*r0mr_over_rnorm2[j-n4];
        Mk[(k+1)*n+j] = (r0dotr[j-n4]*Mk[k*n+j]+k*Nk[(k-1)*n+j]-r0mr_inv[j-n4])*rnorm2_inv[j-n4];
      }
    }
  }  
}
#else
void csimd256lap3ddlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *A) {}
void csimd256line3omega0123vr_c_(double *r0, double *m_all, int64_t *Nbd, int64_t *Order, double *onm0, double *onm1, double *onm2, double *onm3, int64_t *H_dim, double *omega0, double *omega1, double *omega2, double *omega3) {}
void csimd256momentsad_vr_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk) {}
void csimd256momentsad_vrv2_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk) {}
#endif



/*=====================================================================================================*/
#ifdef ARCH_AVX512
void csimd512lap3ddlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *A) 
{
  int64_t m, n, i, j, m8;
  double pi4inv = 1.0/(16.0*atan(1.0));
  m = *M;
  n = *N;
  m8 = (m/8)*8;
  // 
  double *r0x = (double*)aligned_alloc(64, m8 * sizeof(double));
  double *r0y = (double*)aligned_alloc(64, m8 * sizeof(double));
  double *r0z = (double*)aligned_alloc(64, m8 * sizeof(double));
  // 
  for (i = 0; i < m8; ++i) {
    r0x[i] = r0[3*i];
    r0y[i] = r0[3*i+1];
    r0z[i] = r0[3*i+2];
  }
  //   
  __m512d half = _mm512_set1_pd(0.5);
  __m512d three = _mm512_set1_pd(3.0); // for newton iteration
  for (j = 0; j < n; ++j) {
    __m512d rx_v  = _mm512_set1_pd(r[3*j]);
    __m512d ry_v  = _mm512_set1_pd(r[3*j+1]);
    __m512d rz_v  = _mm512_set1_pd(r[3*j+2]);
    __m512d rnx_v = _mm512_set1_pd(rn[3*j]);
    __m512d rny_v = _mm512_set1_pd(rn[3*j+1]);
    __m512d rnz_v = _mm512_set1_pd(rn[3*j+2]);
    __m512d w_v   = _mm512_set1_pd(pi4inv*w[j]);
    for (i = 0; i <= m8 - 8; i += 8) {
      __m512d r0x_v = _mm512_load_pd(&r0x[i]);
      __m512d r0y_v = _mm512_load_pd(&r0y[i]);
      __m512d r0z_v = _mm512_load_pd(&r0z[i]);
      __m512d dx_v = _mm512_sub_pd(r0x_v, rx_v);
      __m512d dy_v = _mm512_sub_pd(r0y_v, ry_v);
      __m512d dz_v = _mm512_sub_pd(r0z_v, rz_v);
      __m512d rr_v = _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(dx_v, dx_v), _mm512_mul_pd(dy_v, dy_v)), _mm512_mul_pd(dz_v, dz_v));
      __m512d rdotrn_v = _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(dx_v, rnx_v), _mm512_mul_pd(dy_v, rny_v)), _mm512_mul_pd(dz_v, rnz_v));
      // modified based on wen yan's DemoSIMD: https://github.com/wenyan4work/DemoSIMD/blob/master/src/rsqrt.cpp
      __m512d ri_v = _mm512_cvtps_pd(_mm256_rsqrt_ps(_mm512_cvtpd_ps(rr_v))); // inverse sqrt approximate
      __m512d muls = _mm512_mul_pd(_mm512_mul_pd(rr_v, ri_v), ri_v);
      ri_v = _mm512_mul_pd(_mm512_mul_pd(half, ri_v), _mm512_sub_pd(three, muls)); // newton iteration
      muls = _mm512_mul_pd(_mm512_mul_pd(rr_v, ri_v), ri_v);
      ri_v = _mm512_mul_pd(_mm512_mul_pd(half, ri_v), _mm512_sub_pd(three, muls)); // newton iteration
      rr_v = _mm512_mul_pd(_mm512_mul_pd(ri_v,ri_v),ri_v); // ri_v**3 misuse of variable name
      __m512d dlpij_v  = _mm512_mul_pd(_mm512_mul_pd(w_v,rdotrn_v), rr_v);
      // __m512d dlpij_v  = _mm512_div_pd(_mm512_mul_pd(w_v,rdotrn_v), _mm512_mul_pd(_mm512_sqrt_pd(rr_v),rr_v));
      // _mm512_store_pd(&A[i+j*m], inv_sqrt_rr_v); // this one crashes
      _mm512_storeu_pd(&A[i+j*m], dlpij_v); // Moves to unaligned memory location
    }
    // compute the rest, simd over source? or not worth it
    for (i=m8; i<m; i++){
      double dx = r0[3*i]   - r[3*j];
      double dy = r0[3*i+1] - r[3*j+1];
      double dz = r0[3*i+2] - r[3*j+2];
      double rr = dx*dx + dy*dy + dz*dz;
      A[m*j+i] = pi4inv*w[j]*(dx*rn[3*j]+dy*rn[3*j+1]+dz*rn[3*j+2])/(sqrt(rr)*rr);
    }
  }
  free(r0x);
  free(r0y);
  free(r0z);
}
void csimd512lap3dslpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *w, double *A) 
{
  int64_t m, n, i, j, m8;
  double pi4inv = 1.0/(16.0*atan(1.0));
  m = *M;
  n = *N;
  m8 = (m/8)*8;
  // 
  double *r0x = (double*)aligned_alloc(64, m8 * sizeof(double));
  double *r0y = (double*)aligned_alloc(64, m8 * sizeof(double));
  double *r0z = (double*)aligned_alloc(64, m8 * sizeof(double));
  // 
  for (i = 0; i < m8; ++i) {
    r0x[i] = r0[3*i];
    r0y[i] = r0[3*i+1];
    r0z[i] = r0[3*i+2];
  }
  //   
  __m512d half = _mm512_set1_pd(0.5);
  __m512d three = _mm512_set1_pd(3.0); // for newton iteration
  for (j = 0; j < n; ++j) {
    __m512d rx_v  = _mm512_set1_pd(r[3*j]);
    __m512d ry_v  = _mm512_set1_pd(r[3*j+1]);
    __m512d rz_v  = _mm512_set1_pd(r[3*j+2]);
    __m512d w_v   = _mm512_set1_pd(pi4inv*w[j]);
    for (i = 0; i <= m8 - 8; i += 8) {
      __m512d r0x_v = _mm512_load_pd(&r0x[i]);
      __m512d r0y_v = _mm512_load_pd(&r0y[i]);
      __m512d r0z_v = _mm512_load_pd(&r0z[i]);
      __m512d dx_v = _mm512_sub_pd(r0x_v, rx_v);
      __m512d dy_v = _mm512_sub_pd(r0y_v, ry_v);
      __m512d dz_v = _mm512_sub_pd(r0z_v, rz_v);
      __m512d rr_v = _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(dx_v, dx_v), _mm512_mul_pd(dy_v, dy_v)), _mm512_mul_pd(dz_v, dz_v));
      // modified based on wen yan's DemoSIMD: https://github.com/wenyan4work/DemoSIMD/blob/master/src/rsqrt.cpp
      __m512d ri_v = _mm512_cvtps_pd(_mm256_rsqrt_ps(_mm512_cvtpd_ps(rr_v))); // inverse sqrt approximate
      __m512d muls = _mm512_mul_pd(_mm512_mul_pd(rr_v, ri_v), ri_v);
      ri_v = _mm512_mul_pd(_mm512_mul_pd(half, ri_v), _mm512_sub_pd(three, muls)); // newton iteration
      muls = _mm512_mul_pd(_mm512_mul_pd(rr_v, ri_v), ri_v);
      ri_v = _mm512_mul_pd(_mm512_mul_pd(half, ri_v), _mm512_sub_pd(three, muls)); // newton iteration
      __m512d slpij_v  = _mm512_mul_pd(w_v, ri_v);
      _mm512_storeu_pd(&A[i+j*m], slpij_v); // Moves to unaligned memory location
    }
    // compute the rest, simd over source? or not worth it
    for (i=m8; i<m; i++){
      double dx = r0[3*i]   - r[3*j];
      double dy = r0[3*i+1] - r[3*j+1];
      double dz = r0[3*i+2] - r[3*j+2];
      double rr = dx*dx + dy*dy + dz*dz;
      A[m*j+i] = pi4inv*w[j]/sqrt(rr);
    }
  }
  free(r0x);
  free(r0y);
  free(r0z);
}
void csimd512lap3dsdlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *As, double *Ad) 
{
  int64_t m, n, i, j, m8;
  double pi4inv = 1.0/(16.0*atan(1.0));
  m = *M;
  n = *N;
  m8 = (m/8)*8;
  // 
  double *r0x = (double*)aligned_alloc(64, m8 * sizeof(double));
  double *r0y = (double*)aligned_alloc(64, m8 * sizeof(double));
  double *r0z = (double*)aligned_alloc(64, m8 * sizeof(double));
  // 
  for (i = 0; i < m8; ++i) {
    r0x[i] = r0[3*i];
    r0y[i] = r0[3*i+1];
    r0z[i] = r0[3*i+2];
  }
  //   
  __m512d half = _mm512_set1_pd(0.5);
  __m512d three = _mm512_set1_pd(3.0); // for newton iteration
  for (j = 0; j < n; ++j) {
    __m512d rx_v  = _mm512_set1_pd(r[3*j]);
    __m512d ry_v  = _mm512_set1_pd(r[3*j+1]);
    __m512d rz_v  = _mm512_set1_pd(r[3*j+2]);
    __m512d rnx_v = _mm512_set1_pd(rn[3*j]);
    __m512d rny_v = _mm512_set1_pd(rn[3*j+1]);
    __m512d rnz_v = _mm512_set1_pd(rn[3*j+2]);
    __m512d w_v   = _mm512_set1_pd(pi4inv*w[j]);
    for (i = 0; i <= m8 - 8; i += 8) {
      __m512d r0x_v = _mm512_load_pd(&r0x[i]);
      __m512d r0y_v = _mm512_load_pd(&r0y[i]);
      __m512d r0z_v = _mm512_load_pd(&r0z[i]);
      __m512d dx_v = _mm512_sub_pd(r0x_v, rx_v);
      __m512d dy_v = _mm512_sub_pd(r0y_v, ry_v);
      __m512d dz_v = _mm512_sub_pd(r0z_v, rz_v);
      __m512d rr_v = _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(dx_v, dx_v), _mm512_mul_pd(dy_v, dy_v)), _mm512_mul_pd(dz_v, dz_v));
      __m512d rdotrn_v = _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(dx_v, rnx_v), _mm512_mul_pd(dy_v, rny_v)), _mm512_mul_pd(dz_v, rnz_v));
      // modified based on wen yan's DemoSIMD: https://github.com/wenyan4work/DemoSIMD/blob/master/src/rsqrt.cpp
      __m512d ri_v = _mm512_cvtps_pd(_mm256_rsqrt_ps(_mm512_cvtpd_ps(rr_v))); // inverse sqrt approximate
      __m512d muls = _mm512_mul_pd(_mm512_mul_pd(rr_v, ri_v), ri_v);
      ri_v = _mm512_mul_pd(_mm512_mul_pd(half, ri_v), _mm512_sub_pd(three, muls)); // newton iteration
      muls = _mm512_mul_pd(_mm512_mul_pd(rr_v, ri_v), ri_v);
      ri_v = _mm512_mul_pd(_mm512_mul_pd(half, ri_v), _mm512_sub_pd(three, muls)); // newton iteration
      rr_v = _mm512_mul_pd(_mm512_mul_pd(ri_v,ri_v),ri_v); // ri_v**3 misuse of variable name
      __m512d dlpij_v  = _mm512_mul_pd(_mm512_mul_pd(w_v,rdotrn_v), rr_v);
      __m512d slpij_v  = _mm512_mul_pd(w_v, ri_v);
      // __m512d dlpij_v  = _mm512_div_pd(_mm512_mul_pd(w_v,rdotrn_v), _mm512_mul_pd(_mm512_sqrt_pd(rr_v),rr_v));
      // _mm512_store_pd(&A[i+j*m], inv_sqrt_rr_v); // this one crashes
      _mm512_storeu_pd(&Ad[i+j*m], dlpij_v); // Moves to unaligned memory location
      _mm512_storeu_pd(&As[i+j*m], slpij_v);
    }
    // compute the rest, simd over source? or not worth it
    for (i=m8; i<m; i++){
      double dx = r0[3*i]   - r[3*j];
      double dy = r0[3*i+1] - r[3*j+1];
      double dz = r0[3*i+2] - r[3*j+2];
      double rr = dx*dx + dy*dy + dz*dz;
      Ad[m*j+i] = pi4inv*w[j]*(dx*rn[3*j]+dy*rn[3*j+1]+dz*rn[3*j+2])/(sqrt(rr)*rr);
      As[m*j+i] = pi4inv*w[j]/sqrt(rr);
    }
  }
  free(r0x);
  free(r0y);
  free(r0z);
}
void csimd512line3omega0123vr_c_(double *r0, double *m_all, int64_t *Nbd, int64_t *Order, double *onm0, double *onm1, double *onm2, double *onm3, int64_t *H_dim, double *omega0, double *omega1, double *omega2, double *omega3)
{
  int64_t nbd, order, h_dim, i, ii, j, l, nbd4, nbd48;
  double pi, tmpval, r00, r01, r02;
  nbd = *Nbd;
  order = *Order;
  h_dim = *H_dim;
  nbd4 = 4*nbd;
  nbd48 = (nbd4/8)*8;
  // 
  double *mk_j_ii = (double*)aligned_alloc(64, 4 * nbd * sizeof(double));
  pi = 4.0*atan(1.0);
  tmpval = -1.0/(4.0*pi);
  r00 = tmpval*r0[0];
  r01 = tmpval*r0[1];
  r02 = tmpval*r0[2];
  // 
  for (i = 0; i < h_dim; ++i){
    omega0[i] = 0.0;
    omega1[i] = 0.0;
    omega2[i] = 0.0;
    omega3[i] = 0.0;
  }
  // loop over moment order
  i = 0;
  for (ii = 0; ii < order+1; ++ii){
    // 
    for (j = 0; j < nbd; ++j){
      mk_j_ii[4*j]   = tmpval*m_all[(ii+1)*nbd+j];
      mk_j_ii[4*j+1] = r00*m_all[ii*nbd+j];
      mk_j_ii[4*j+2] = r01*m_all[ii*nbd+j];
      mk_j_ii[4*j+3] = r02*m_all[ii*nbd+j];
    }
    // 
    for (l = 0; l < ii; ++l){
      __m512d omega0_v = _mm512_set1_pd(0.0);
      __m512d omega1_v = _mm512_set1_pd(0.0);
      __m512d omega2_v = _mm512_set1_pd(0.0);
      __m512d omega3_v = _mm512_set1_pd(0.0);
      for (j = 0; j <= nbd48 - 8; j += 8){
        __m512d mk = _mm512_load_pd(&mk_j_ii[j]);
        __m512d onm0_v = _mm512_load_pd(&onm0[nbd4*i+j]);
        __m512d onm1_v = _mm512_load_pd(&onm1[nbd4*i+j]);
        __m512d onm2_v = _mm512_load_pd(&onm2[nbd4*i+j]);
        __m512d onm3_v = _mm512_load_pd(&onm3[nbd4*i+j]);
        omega0_v = _mm512_fmadd_pd(mk, onm0_v, omega0_v);
        omega1_v = _mm512_fmadd_pd(mk, onm1_v, omega1_v);
        omega2_v = _mm512_fmadd_pd(mk, onm2_v, omega2_v);
        omega3_v = _mm512_fmadd_pd(mk, onm3_v, omega3_v);
      }
      omega0[i] = _mm512_reduce_add_pd(omega0_v);
      omega1[i] = _mm512_reduce_add_pd(omega1_v);
      omega2[i] = _mm512_reduce_add_pd(omega2_v);
      omega3[i] = _mm512_reduce_add_pd(omega3_v);
      // compute the rest if needed
      for (j=nbd48; j<nbd4; j++){
        omega0[i] +=  mk_j_ii[j] * onm0[nbd4*i+j]; 
        omega1[i] +=  mk_j_ii[j] * onm1[nbd4*i+j]; 
        omega2[i] +=  mk_j_ii[j] * onm2[nbd4*i+j]; 
        omega3[i] +=  mk_j_ii[j] * onm3[nbd4*i+j]; 
      }
      ++i;
    }
  }
  free(mk_j_ii);
}
void csimd512momentsad_vr_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk)
{
  int64_t n = *N, order = *Order, k, j, n8;
  double r0norm, r0norm_inv, tmp, rkx, rky, rkz, rmr0x, rmr0y, rmr0z;
  double LMNcommon, rnorm2, rnorm, rnorm_inv, r0normplusrnorm, r0mr;
  double *r0dotr = (double*)aligned_alloc(64, n * sizeof(double));
  double *rnorm2_inv = (double*)aligned_alloc(64, n * sizeof(double));
  double *r0dotr_over_rnorm2 = (double*)aligned_alloc(64, n * sizeof(double));
  double *r0norm2_over_rnorm2 = (double*)aligned_alloc(64, n * sizeof(double));
  double *r0mr_inv = (double*)aligned_alloc(64, n * sizeof(double));
  double *r0mr_over_rnorm2 = (double*)aligned_alloc(64, n * sizeof(double));
  // 
  n8 = (n/8)*8;
  r0norm = sqrt(r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2]);
  r0norm_inv = 1.0/r0norm;
  // 
  double *rx = (double*)aligned_alloc(64, n8 * sizeof(double));
  double *ry = (double*)aligned_alloc(64, n8 * sizeof(double));
  double *rz = (double*)aligned_alloc(64, n8 * sizeof(double));
  // to call _mm512_load_pd, otherwise call _mm512_loadu_pd...
  for (j = 0; j < n8; ++j) {
    rx[j] = r[3*j];
    ry[j] = r[3*j+1];
    rz[j] = r[3*j+2];
  }
  // 
  __m512d r0x_v = _mm512_set1_pd(r0[0]);
  __m512d r0y_v = _mm512_set1_pd(r0[1]);
  __m512d r0z_v = _mm512_set1_pd(r0[2]);
  __m512d r0norm_v = _mm512_set1_pd(r0norm);
  __m512d r0norm_inv_v = _mm512_set1_pd(r0norm_inv);
  __m512d one_v = _mm512_set1_pd(1.0);
  // 
  for (k = 0; k < n8; k += 8) {
    // 
    __m512d rx_v = _mm512_load_pd(&rx[k]);
    __m512d ry_v = _mm512_load_pd(&ry[k]);
    __m512d rz_v = _mm512_load_pd(&rz[k]);
    __m512d rmr0x_v = _mm512_sub_pd(r0x_v, rx_v);
    __m512d rmr0y_v = _mm512_sub_pd(r0y_v, ry_v);
    __m512d rmr0z_v = _mm512_sub_pd(r0z_v, rz_v);
    __m512d rnorm2_v = _mm512_fmadd_pd(rx_v, rx_v, _mm512_fmadd_pd(ry_v, ry_v, _mm512_mul_pd(rz_v, rz_v)));
    __m512d rnorm_v = _mm512_sqrt_pd(rnorm2_v);
    __m512d rnorm_inv_v = _mm512_div_pd(one_v, rnorm_v);
    __m512d r0normplusrnorm_v = _mm512_add_pd(r0norm_v, rnorm_v);
    __m512d r0mr_v = _mm512_sqrt_pd(_mm512_fmadd_pd(rmr0x_v, rmr0x_v, _mm512_fmadd_pd(rmr0y_v, rmr0y_v, _mm512_mul_pd(rmr0z_v, rmr0z_v))));
    __m512d r0dotr_v = _mm512_fmadd_pd(r0x_v, rx_v, _mm512_fmadd_pd(r0y_v, ry_v, _mm512_mul_pd(r0z_v, rz_v)));
    __m512d rnorm2_inv_v = _mm512_mul_pd(rnorm_inv_v, rnorm_inv_v);
    __m512d r0dotr_over_rnorm2_v = _mm512_mul_pd(r0dotr_v, rnorm2_inv_v);
    __m512d r0norm2_over_rnorm2_v = _mm512_mul_pd(_mm512_mul_pd(r0norm_v, r0norm_v), rnorm2_inv_v);
    __m512d r0mr_inv_v = _mm512_div_pd(one_v, r0mr_v);
    __m512d r0mr_over_rnorm2_v = _mm512_mul_pd(r0mr_v, rnorm2_inv_v);
    // 
    _mm512_store_pd(&r0dotr[k], r0dotr_v);
    _mm512_store_pd(&rnorm2_inv[k], rnorm2_inv_v);
    _mm512_store_pd(&r0dotr_over_rnorm2[k], r0dotr_over_rnorm2_v);
    _mm512_store_pd(&r0norm2_over_rnorm2[k], r0norm2_over_rnorm2_v);
    _mm512_store_pd(&r0mr_inv[k], r0mr_inv_v);
    _mm512_store_pd(&r0mr_over_rnorm2[k], r0mr_over_rnorm2_v);
    // 
    __m512d LMNcommon_v = _mm512_div_pd(r0normplusrnorm_v, _mm512_add_pd(r0normplusrnorm_v, _mm512_mul_pd(r0mr_inv_v, _mm512_sub_pd(_mm512_mul_pd(r0norm_v, r0norm_v), rnorm2_v))));
    __m512d tmp_v = _mm512_mul_pd(_mm512_mul_pd(_mm512_add_pd(r0normplusrnorm_v, r0mr_v), LMNcommon_v), r0mr_inv_v);
    //// __m512d log_tmp_v = _mm512_log_pd(tmp_v); // no such thing in gcc ... look it up... 07/20/24 Hai
    __m512d log_tmp_v = {log(tmp_v[0]), log(tmp_v[1]), log(tmp_v[2]), log(tmp_v[3]),
                         log(tmp_v[4]), log(tmp_v[5]), log(tmp_v[6]), log(tmp_v[7])};
    __m512d Nk_v = _mm512_mul_pd(log_tmp_v, rnorm_inv_v);
    _mm512_storeu_pd(&Nk[k], Nk_v);
    __m512d Nk_next_v = _mm512_fmadd_pd(Nk_v, r0dotr_over_rnorm2_v, _mm512_mul_pd(_mm512_sub_pd(r0mr_v, r0norm_v), rnorm2_inv_v));
    _mm512_storeu_pd(&Nk[n + k], Nk_next_v);
    // 
    LMNcommon_v = _mm512_div_pd(LMNcommon_v, _mm512_fmadd_pd(r0norm_v, rnorm_v, r0dotr_v));
    __m512d Mk_v = _mm512_mul_pd(LMNcommon_v,_mm512_sub_pd(_mm512_mul_pd(_mm512_add_pd(_mm512_mul_pd(r0normplusrnorm_v, r0mr_inv_v), _mm512_mul_pd(rnorm_v, r0norm_inv_v)), r0mr_inv_v), r0norm_inv_v));
    _mm512_storeu_pd(&Mk[k], Mk_v);
    __m512d Mk_next_v = _mm512_div_pd(_mm512_mul_pd(r0norm_v, Mk_v), _mm512_add_pd(r0norm_v, r0mr_v));
    _mm512_storeu_pd(&Mk[n + k], Mk_next_v);
  }
  // fix the rest
  for (k=n8; k<n; k++){
    rkx = r[3*k];
    rky = r[3*k+1];
    rkz = r[3*k+2];
    rmr0x = r0[0] - rkx; 
    rmr0y = r0[1] - rky; 
    rmr0z = r0[2] - rkz;
    rnorm2 = rkx*rkx + rky*rky + rkz*rkz;
    rnorm = sqrt(rnorm2);
    rnorm_inv = 1.0/rnorm;
    r0normplusrnorm = r0norm + rnorm;
    r0mr = sqrt(rmr0x*rmr0x + rmr0y*rmr0y + rmr0z*rmr0z);
    // to be used later
    r0dotr[k] = r0[0]*rkx + r0[1]*rky + r0[2]*rkz;
    rnorm2_inv[k] = rnorm_inv*rnorm_inv;
    r0dotr_over_rnorm2[k] = r0dotr[k]*rnorm2_inv[k];
    r0norm2_over_rnorm2[k] = (r0norm*r0norm)*rnorm2_inv[k];
    r0mr_inv[k] = 1.0/r0mr;
    r0mr_over_rnorm2[k] = r0mr*rnorm2_inv[k];
    // 
    LMNcommon = r0normplusrnorm/(r0normplusrnorm
                  + r0mr_inv[k]*( r0norm*r0norm - rnorm2));
    tmp = (r0normplusrnorm+r0mr)*LMNcommon*r0mr_inv[k];
    Nk[k] = log(tmp)*rnorm_inv;
    Nk[n+k] = Nk[k]*r0dotr_over_rnorm2[k] 
              + (r0mr-r0norm)*rnorm2_inv[k];
    // 
    LMNcommon = LMNcommon/(r0norm*rnorm+r0dotr[k]);
    Mk[k] = LMNcommon*((r0normplusrnorm* 
        r0mr_inv[k]+rnorm*r0norm_inv)*r0mr_inv[k]-r0norm_inv);
    Mk[n+k] = r0norm*Mk[k]/(r0norm+r0mr);
  }
  for (k = 1; k < order; ++k) {
    for (j = 0; j < n8; j += 8) {
      __m512d Nk_v = _mm512_loadu_pd(&Nk[k*n + j]);
      __m512d Nk_prev_v = _mm512_loadu_pd(&Nk[(k-1)*n + j]);
      __m512d Mk_v = _mm512_loadu_pd(&Mk[k*n + j]);
      __m512d r0dotr_v = _mm512_load_pd(&r0dotr[j]);
      __m512d rnorm2_inv_v = _mm512_load_pd(&rnorm2_inv[j]);
      __m512d r0dotr_over_rnorm2_v = _mm512_load_pd(&r0dotr_over_rnorm2[j]);
      __m512d r0norm2_over_rnorm2_v = _mm512_load_pd(&r0norm2_over_rnorm2[j]);
      __m512d r0mr_inv_v = _mm512_load_pd(&r0mr_inv[j]);
      __m512d r0mr_over_rnorm2_v = _mm512_load_pd(&r0mr_over_rnorm2[j]);
      // 
      __m512d Nk_next_v = _mm512_div_pd(_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(2*k+1), _mm512_mul_pd(r0dotr_over_rnorm2_v, Nk_v)), _mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(k), _mm512_mul_pd(r0norm2_over_rnorm2_v, Nk_prev_v)), r0mr_over_rnorm2_v)), _mm512_set1_pd(k+1));
      _mm512_storeu_pd(&Nk[(k+1)*n + j], Nk_next_v);
      __m512d Mk_next_v = _mm512_mul_pd(_mm512_sub_pd(_mm512_fmadd_pd(_mm512_set1_pd(k), Nk_prev_v, _mm512_mul_pd(r0dotr_v, Mk_v)), r0mr_inv_v), rnorm2_inv_v);
      _mm512_storeu_pd(&Mk[(k+1)*n + j], Mk_next_v);
    }
    // fix the rest
    for (j = n8; j < n; ++j){
      Nk[(k+1)*n+j] = ((2*k+1)*r0dotr_over_rnorm2[j]*Nk[k*n+j] -
                        k*r0norm2_over_rnorm2[j]*Nk[(k-1)*n+j] +
                      r0mr_over_rnorm2[j])/(k+1);
      Mk[(k+1)*n+j] = (r0dotr[j]*Mk[k*n+j]+k*Nk[(k-1)*n+j]-r0mr_inv[j])*rnorm2_inv[j];
    }
  }

  free(r0dotr);
  free(rnorm2_inv);
  free(r0dotr_over_rnorm2);
  free(r0norm2_over_rnorm2);
  free(r0mr_inv);
  free(r0mr_over_rnorm2);
}
// this one is ok to use, avx512 only gives 50% speedup when n = n8, look at assist_simd_momentsad_vr, returns Nk8, Mk8
// otherwise, speed is about the same...
void csimd512momentsad_vrv2_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk)
{
  int64_t n = *N, order = *Order, k, j, n8, nmn8, idx;
  double tmpk1, tmpk2;
  double r0norm2, r0norm, r0norm_inv;
  n8 = (n/8)*8;
  double *r8x = (double*)aligned_alloc(64, n8 * sizeof(double));
  double *r8y = (double*)aligned_alloc(64, n8 * sizeof(double));
  double *r8z = (double*)aligned_alloc(64, n8 * sizeof(double));
  double *tmpk3 = (double*)aligned_alloc(64, order * sizeof(double));
  r0norm2 = (r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2]);
  r0norm  = sqrt(r0norm2);
  r0norm_inv = 1.0/r0norm;
  __m512d r0x_v = _mm512_set1_pd(r0[0]);
  __m512d r0y_v = _mm512_set1_pd(r0[1]);
  __m512d r0z_v = _mm512_set1_pd(r0[2]);
  __m512d r0norm_v = _mm512_set1_pd(r0norm);
  __m512d r0norm2_v = _mm512_set1_pd(r0norm2);
  __m512d r0norm_inv_v = _mm512_set1_pd(r0norm_inv);
  __m512d one_v = _mm512_set1_pd(1.0);
  for (j = 0; j < n8; ++j) {
    r8x[j] = r[3*j];
    r8y[j] = r[3*j+1];
    r8z[j] = r[3*j+2];
  }
  for (k = 1; k < order; ++k) {
    tmpk3[k] = 1.0/(k+1);
  }
  idx = 0;
  __m512d Nnext_v = _mm512_set1_pd(0.0);
  __m512d Mnext_v = _mm512_set1_pd(0.0);
  for (j = 0; j < n8; j += 8){
    __m512d rx_v = _mm512_load_pd(&r8x[j]);
    __m512d ry_v = _mm512_load_pd(&r8y[j]);
    __m512d rz_v = _mm512_load_pd(&r8z[j]);
    __m512d rmr0x_v = _mm512_sub_pd(r0x_v, rx_v);
    __m512d rmr0y_v = _mm512_sub_pd(r0y_v, ry_v);
    __m512d rmr0z_v = _mm512_sub_pd(r0z_v, rz_v);
    // 
    __m512d rnorm2_v = _mm512_fmadd_pd(rx_v, rx_v, _mm512_fmadd_pd(ry_v, ry_v, _mm512_mul_pd(rz_v, rz_v)));
    __m512d rnorm_v = _mm512_sqrt_pd(rnorm2_v);
    __m512d rnorm_inv_v = _mm512_div_pd(one_v, rnorm_v);
    __m512d rnorm2_inv_v = _mm512_mul_pd(rnorm_inv_v, rnorm_inv_v);
    __m512d r0mr_v = _mm512_sqrt_pd(_mm512_fmadd_pd(rmr0x_v, rmr0x_v, _mm512_fmadd_pd(rmr0y_v, rmr0y_v, _mm512_mul_pd(rmr0z_v, rmr0z_v))));
    __m512d r0mr_inv_v = _mm512_div_pd(one_v, r0mr_v);
    __m512d r0dotr_v = _mm512_fmadd_pd(r0x_v, rx_v, _mm512_fmadd_pd(r0y_v, ry_v, _mm512_mul_pd(r0z_v, rz_v)));
    // 
    __m512d rnorm2mr0dotr_v = _mm512_sub_pd(rnorm2_v, r0dotr_v);
    __m512d r0rmr0dotr_inv_v = _mm512_div_pd(one_v, _mm512_fmsub_pd(r0norm_v, rnorm_v, r0dotr_v));
    __m512d tmp_v = _mm512_mul_pd(_mm512_fmadd_pd(rnorm_v, r0mr_v, rnorm2mr0dotr_v), r0rmr0dotr_inv_v);
    __m512d Nprev_v = {log(tmp_v[0]), log(tmp_v[1]), log(tmp_v[2]), log(tmp_v[3]),
                         log(tmp_v[4]), log(tmp_v[5]), log(tmp_v[6]), log(tmp_v[7])};
    Nprev_v = _mm512_mul_pd(rnorm_inv_v, Nprev_v);
    __m512d Ncurr_v = _mm512_mul_pd(rnorm2_inv_v, _mm512_fmadd_pd(r0dotr_v, Nprev_v, _mm512_sub_pd(r0mr_v, r0norm_v)));
    // 
    _mm512_storeu_pd(&Nk[j], Nprev_v);
    _mm512_storeu_pd(&Nk[n + j], Ncurr_v);
    // 
    __m512d r0rpr0dotr_inv_v = _mm512_div_pd(one_v,_mm512_fmadd_pd(r0norm_v, rnorm_v, r0dotr_v));
    __m512d term1_v = _mm512_mul_pd(r0rpr0dotr_inv_v, r0rmr0dotr_inv_v);
    __m512d Mprev_v = _mm512_mul_pd(term1_v,_mm512_fmadd_pd(rnorm2mr0dotr_v, r0mr_inv_v, _mm512_mul_pd(r0dotr_v,r0norm_inv_v)));
    __m512d Mcurr_v = _mm512_mul_pd(term1_v, _mm512_fmadd_pd(_mm512_sub_pd(r0dotr_v, r0norm2_v), r0mr_inv_v, r0norm_v));
    // 
    _mm512_storeu_pd(&Mk[j], Mprev_v);
    _mm512_storeu_pd(&Mk[n + j], Mcurr_v);
    for (k = 1; k < order; ++k) {
      Nnext_v = _mm512_mul_pd(_mm512_fmadd_pd( _mm512_set1_pd(2*k+1), _mm512_mul_pd(r0dotr_v, Ncurr_v), 
                              _mm512_fmadd_pd( _mm512_set1_pd(-k), _mm512_mul_pd(r0norm2_v, Nprev_v), r0mr_v)),
                              _mm512_mul_pd(rnorm2_inv_v, _mm512_set1_pd(tmpk3[k])));
      Mnext_v = _mm512_mul_pd(_mm512_fmadd_pd(r0dotr_v, Mcurr_v, 
                              _mm512_fmsub_pd(_mm512_set1_pd(k), Nprev_v, r0mr_inv_v)), 
                                rnorm2_inv_v);
      _mm512_storeu_pd(&Nk[(k+1)*n + j], Nnext_v);
      _mm512_storeu_pd(&Mk[(k+1)*n + j], Mnext_v);
      Nprev_v = Ncurr_v;
      Ncurr_v = Nnext_v;
      Mcurr_v = Mnext_v;
    }
  }
  nmn8 = n - n8;
  if (nmn8 > 0) {
    double tmp, rx, ry, rz, rmr0x, rmr0y, rmr0z;
    double LMNcommon, rnorm2, rnorm, rnorm_inv, r0normplusrnorm, r0mr;
    double *r0dotr = (double*)aligned_alloc(64, nmn8 * sizeof(double));
    double *rnorm2_inv = (double*)aligned_alloc(64, nmn8 * sizeof(double));
    double *r0dotr_over_rnorm2 = (double*)aligned_alloc(64, nmn8 * sizeof(double));
    double *r0norm2_over_rnorm2 = (double*)aligned_alloc(64, nmn8 * sizeof(double));
    double *r0mr_inv = (double*)aligned_alloc(64, nmn8 * sizeof(double));
    double *r0mr_over_rnorm2 = (double*)aligned_alloc(64, nmn8 * sizeof(double));
    for (k=n8; k<n; k++){
      rx = r[3*k];
      ry = r[3*k+1];
      rz = r[3*k+2];
      rmr0x = r0[0] - rx; 
      rmr0y = r0[1] - ry; 
      rmr0z = r0[2] - rz;
      rnorm2 = rx*rx + ry*ry + rz*rz;
      rnorm = sqrt(rnorm2);
      rnorm_inv = 1.0/rnorm;
      r0normplusrnorm = r0norm + rnorm;
      r0mr = sqrt(rmr0x*rmr0x + rmr0y*rmr0y + rmr0z*rmr0z);
      // to be used later
      r0dotr[k-n8] = r0[0]*rx + r0[1]*ry + r0[2]*rz;
      rnorm2_inv[k-n8] = rnorm_inv*rnorm_inv;
      r0dotr_over_rnorm2[k-n8] = r0dotr[k-n8]*rnorm2_inv[k-n8];
      r0norm2_over_rnorm2[k-n8] = (r0norm*r0norm)*rnorm2_inv[k-n8];
      r0mr_inv[k-n8] = 1.0/r0mr;
      r0mr_over_rnorm2[k-n8] = r0mr*rnorm2_inv[k-n8];
      // 
      LMNcommon = r0normplusrnorm/(r0normplusrnorm
                    + r0mr_inv[k-n8]*( r0norm*r0norm - rnorm2));
      tmp = (r0normplusrnorm+r0mr)*LMNcommon*r0mr_inv[k-n8];
      Nk[k] = log(tmp)*rnorm_inv;
      Nk[n+k] = Nk[k]*r0dotr_over_rnorm2[k-n8] 
                + (r0mr-r0norm)*rnorm2_inv[k-n8];
      // 
      LMNcommon = LMNcommon/(r0norm*rnorm+r0dotr[k-n8]);
      Mk[k] = LMNcommon*((r0normplusrnorm* 
        r0mr_inv[k-n8]+rnorm*r0norm_inv)*r0mr_inv[k-n8]-r0norm_inv);
      Mk[n+k] = r0norm*Mk[k]/(r0norm+r0mr);
    }
    for (k = 1; k < order; ++k){
      tmpk1 = (2*k+1)*tmpk3[k];
      tmpk2 = k*tmpk3[k];
      for (j = n8; j < n; ++j){
        Nk[(k+1)*n+j] = tmpk1*r0dotr_over_rnorm2[j-n8]*Nk[k*n+j] -
                        tmpk2*r0norm2_over_rnorm2[j-n8]*Nk[(k-1)*n+j] +
                        tmpk3[k]*r0mr_over_rnorm2[j-n8];
        Mk[(k+1)*n+j] = (r0dotr[j-n8]*Mk[k*n+j]+k*Nk[(k-1)*n+j]-r0mr_inv[j-n8])*rnorm2_inv[j-n8];
      }
    }
  }
}
#else
void csimd512lap3ddlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *A) {}
void csimd512lap3dslpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *w, double *A) {}
void csimd512lap3dsdlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *As, double *Ad) {}
void csimd512line3omega0123vr_c_(double *r0, double *m_all, int64_t *Nbd, int64_t *Order, double *onm0, double *onm1, double *onm2, double *onm3, int64_t *H_dim, double *omega0, double *omega1, double *omega2, double *omega3) {}
void csimd512momentsad_vr_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk) {}
void csimd512momentsad_vrv2_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk8, double *Mk8) {}
#endif

#endif

/*=====================================================================================================*/
#ifdef ARCH_ARM
#include <arm_neon.h>
void csimd128lap3ddlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *A) 
{
  int64_t m, n, i, j, m2;
  double pi4inv = 1.0/(16.0*atan(1.0));
  m = *M;
  n = *N;
  m2 = (m/2)*2;

  double *r0x = (double*)aligned_alloc(16, m2 * sizeof(double));
  double *r0y = (double*)aligned_alloc(16, m2 * sizeof(double));
  double *r0z = (double*)aligned_alloc(16, m2 * sizeof(double));

  for (i = 0; i < m2; ++i) {
    r0x[i] = r0[3*i];
    r0y[i] = r0[3*i+1];
    r0z[i] = r0[3*i+2];
  }
  // 
  float64x2_t half = vdupq_n_f64(0.5);
  float64x2_t three = vdupq_n_f64(3.0); // for newton iteration
  for (j = 0; j < n; ++j) {
    float64x2_t rx_v = vdupq_n_f64(r[3*j]);
    float64x2_t ry_v = vdupq_n_f64(r[3*j+1]);
    float64x2_t rz_v = vdupq_n_f64(r[3*j+2]);
    float64x2_t rnx_v = vdupq_n_f64(rn[3*j]);
    float64x2_t rny_v = vdupq_n_f64(rn[3*j+1]);
    float64x2_t rnz_v = vdupq_n_f64(rn[3*j+2]);
    float64x2_t w_v   = vdupq_n_f64(pi4inv*w[j]);
    for (i = 0; i <= m2 - 2; i += 2) {
      float64x2_t r0x_v = vld1q_f64(&r0x[i]);
      float64x2_t r0y_v = vld1q_f64(&r0y[i]);
      float64x2_t r0z_v = vld1q_f64(&r0z[i]);
      float64x2_t dx_v = vsubq_f64(r0x_v, rx_v);
      float64x2_t dy_v = vsubq_f64(r0y_v, ry_v);
      float64x2_t dz_v = vsubq_f64(r0z_v, rz_v);
      float64x2_t rr_v = vaddq_f64(vaddq_f64(vmulq_f64(dx_v, dx_v), vmulq_f64(dy_v, dy_v)), vmulq_f64(dz_v, dz_v));
      float64x2_t rdotrn_v = vaddq_f64(vaddq_f64(vmulq_f64(dx_v, rnx_v), vmulq_f64(dy_v, rny_v)), vmulq_f64(dz_v, rnz_v));
      // 1
      // float64x2_t ri_v = vrsqrteq_f64(rr_v); // this estimate is not very accurate
      // 2
      // // modified based on wen yan's DemoSIMD (seems slower on arm float64x2_t): https://github.com/wenyan4work/DemoSIMD/blob/master/src/rsqrt.cpp
      // float64x2_t ri_v = vcvt_f64_f32(vrsqrte_f32(vcvt_f32_f64(rr_v))); // inverse sqrt approximate
      // float64x2_t muls = vmulq_f64(vmulq_f64(rr_v, ri_v), ri_v);
      // ri_v = vmulq_f64(vmulq_f64(half, ri_v), vsubq_f64(three, muls)); // newton iteration
      // muls = vmulq_f64(vmulq_f64(rr_v, ri_v), ri_v);
      // ri_v = vmulq_f64(vmulq_f64(half, ri_v), vsubq_f64(three, muls)); // newton iteration
      // rr_v = vmulq_f64(vmulq_f64(ri_v,ri_v),ri_v); // ri_v**3 misuse of variable name
      // float64x2_t dlpij_v  = vmulq_f64(vmulq_f64(w_v,rdotrn_v), rr_v);
      // 3
      float64x2_t dlpij_v  = vdivq_f64(vmulq_f64(w_v,rdotrn_v), vmulq_f64(vsqrtq_f64(rr_v),rr_v));
      vst1q_f64(&A[i+j*m], dlpij_v);
    }
    for (i = m2; i < m; i++) {
      double dx = r0[3*i]   - r[3*j];
      double dy = r0[3*i+1] - r[3*j+1];
      double dz = r0[3*i+2] - r[3*j+2];
      double rr = dx*dx + dy*dy + dz*dz;
      A[m*j+i] = pi4inv*w[j]*(dx*rn[3*j]+dy*rn[3*j+1]+dz*rn[3*j+2])/(sqrt(rr)*rr);
    }
  }
  free(r0x);
  free(r0y);
  free(r0z);
}
// is there such a thing for ARM?
void csimd256lap3ddlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *A) {}
void csimd512lap3ddlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *A) {}
void csimd256lap3dslpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *w, double *A) {}
void csimd512lap3dslpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *w, double *A) {}
void csimd256lap3dsdlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *As, double *Ad) {}
void csimd512lap3dsdlpmat_c_(int64_t *M, double *r0, int64_t *N, double *r, double *rn, double *w, double *As, double *Ad) {}
void csimd128line3omega0123vr_c_(double *r0, double *m_all, int64_t *nbd, int64_t *order, double *onm0, double *onm1, double *onm2, double *onm3, int64_t *h_dim, double *omega0, double *omega1, double *omega2, double *omega3) 
{
// to be implemented
}
void csimd256line3omega0123vr_c_(double *r0, double *m_all, int64_t *nbd, int64_t *order, double *onm0, double *onm1, double *onm2, double *onm3, int64_t *h_dim, double *omega0, double *omega1, double *omega2, double *omega3) {}
void csimd512line3omega0123vr_c_(double *r0, double *m_all, int64_t *nbd, int64_t *order, double *onm0, double *onm1, double *onm2, double *onm3, int64_t *h_dim, double *omega0, double *omega1, double *omega2, double *omega3) {}

void csimd128momentsad_vr_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk) 
{
// to be implemented
}
void csimd256momentsad_vr_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk) {}
void csimd256momentsad_vrv2_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk) {}
void csimd512momentsad_vr_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk, double *Mk) {}
void csimd512momentsad_vrv2_c_(double *r0, int64_t *N, double *r, int64_t *Order, int64_t *Flag, double *Nk8, double *Mk8) {}

#endif
