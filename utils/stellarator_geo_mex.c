#define PSPLIT_MATCH_MATLAB 0

#include <math.h>
#include <stdlib.h>

#define ASPECT 1.65

typedef long long i8;

static void stellaratorparam(double t, double p, double *x)
{
  double ct = cos(t), st = sin(t);
  double c2pt = cos(2.0*p - t), s2pt = sin(2.0*p - t);
  double c2p  = cos(2.0*p),     s2p  = sin(2.0*p);
  double cp   = cos(p),         sp   = sin(p);
  double cpt  = cos(p + t),     spt  = sin(p + t);
  double cmpt = cos(-p + t),    smpt = sin(-p + t);

  x[0] =  0.17*ct*c2pt + 0.11*ct*c2p + 1.0*ct*cp + 4.5*ct
        - 0.25*ct*cp   + 0.07*ct*cpt - 0.45*ct*cmpt;
  x[1] =  0.17*st*c2pt + 0.11*st*c2p + 1.0*st*cp + 4.5*st
        - 0.25*st*cp   + 0.07*st*cpt - 0.45*st*cmpt;
  x[2] =  0.17*s2pt    + 0.11*s2p    + 1.0*sp    + 0.0
        + 0.25*sp      + 0.07*spt    - 0.45*smpt;
}

typedef struct {
  double a1, b1, a2, b2;
  double jj, ii;
} chart_t;

static void chart_eval(const chart_t *c, double ts1, double ts2,
                       double t1, double t2, double *x)
{
  double u = ((c->a1*t1 + c->b1) + 2.0*c->jj) * ts1;
  double v = ((c->a2*t2 + c->b2) + 2.0*c->ii) * ts2;
  stellaratorparam(u, v, x);
}

static i8 build_charts(i8 mp, i8 np, i8 p, const double *x0,
                       double ts1, double ts2, chart_t *cl)
{
  chart_t c0;
  double xll[3], xul[3], xlr[3], dt, dp;
  i8 i, j, nq = 0;

  c0.a1 = 1.0; c0.b1 = 0.0; c0.a2 = 1.0; c0.b2 = 0.0;
  for (j = 1; j <= np; j++) {
    for (i = 1; i <= mp; i++) {
      c0.jj = (double)j;  c0.ii = (double)i;
      chart_eval(&c0, ts1, ts2, x0[0],   x0[0],   xll);
      chart_eval(&c0, ts1, ts2, x0[0],   x0[p-1], xul);
      chart_eval(&c0, ts1, ts2, x0[p-1], x0[0],   xlr);
      dt = sqrt((xll[0]-xlr[0])*(xll[0]-xlr[0]) + (xll[1]-xlr[1])*(xll[1]-xlr[1])
              + (xll[2]-xlr[2])*(xll[2]-xlr[2]));
      dp = sqrt((xll[0]-xul[0])*(xll[0]-xul[0]) + (xll[1]-xul[1])*(xll[1]-xul[1])
              + (xll[2]-xul[2])*(xll[2]-xul[2]));
      if (dt/dp > ASPECT) {
        if (cl) { cl[nq] = c0; cl[nq].a1 = 0.5; cl[nq].b1 = -0.5; } nq++;
        if (cl) { cl[nq] = c0; cl[nq].a1 = 0.5; cl[nq].b1 =  0.5; } nq++;
      } else if (dp/dt > ASPECT) {
        if (cl) { cl[nq] = c0; cl[nq].a2 = 0.5; cl[nq].b2 = -0.5; } nq++;
#if PSPLIT_MATCH_MATLAB
        if (cl) { cl[nq] = c0; cl[nq].a2 = 0.5; cl[nq].b2 = -0.5; } nq++;
#else
        if (cl) { cl[nq] = c0; cl[nq].a2 = 0.5; cl[nq].b2 =  0.5; } nq++;
#endif
      } else {
        if (cl) cl[nq] = c0; nq++;
      }
    }
  }
  return nq;
}

#define BARYCENTRIC 0

#if BARYCENTRIC

static void interpmat_1d(const double *t, i8 q, const double *s, i8 p, double *L)
{
  i8 i, j, k, hit;
  double *bw = (double*) malloc(sizeof(double)*(size_t)p);
  double den, c;

  for (j = 0; j < p; j++) {
    double prod = 1.0;
    for (k = 0; k < p; k++) if (k != j) prod *= (s[j] - s[k]);
    bw[j] = 1.0/prod;
  }

  for (i = 0; i < q; i++) {
    hit = -1;
    for (j = 0; j < p; j++) if (t[i] == s[j]) { hit = j; break; }
    if (hit >= 0) {
      for (j = 0; j < p; j++) L[i + j*q] = (j == hit) ? 1.0 : 0.0;
      continue;
    }
    den = 0.0;
    for (j = 0; j < p; j++) {
      c = bw[j]/(t[i] - s[j]);
      L[i + j*q] = c;
      den += c;
    }
    for (j = 0; j < p; j++) L[i + j*q] /= den;
  }

  free(bw);
}

#else

static void lu_solve(double *A, i8 p, double *B, i8 q)
{
  i8 k, i, j, ipv;
  double pmax, t;

  for (k = 0; k < p; k++) {
    ipv = k;  pmax = fabs(A[k + k*p]);
    for (i = k+1; i < p; i++)
      if (fabs(A[i + k*p]) > pmax) { pmax = fabs(A[i + k*p]);  ipv = i; }
    if (ipv != k) {
      for (j = 0; j < p; j++) { t = A[k + j*p]; A[k + j*p] = A[ipv + j*p]; A[ipv + j*p] = t; }
      for (j = 0; j < q; j++) { t = B[k + j*p]; B[k + j*p] = B[ipv + j*p]; B[ipv + j*p] = t; }
    }
    for (i = k+1; i < p; i++) {
      double mlt = A[i + k*p]/A[k + k*p];
      A[i + k*p] = mlt;
      for (j = k+1; j < p; j++) A[i + j*p] -= mlt*A[k + j*p];
      for (j = 0;   j < q; j++) B[i + j*p] -= mlt*B[k + j*p];
    }
  }
  for (j = 0; j < q; j++)
    for (i = p-1; i >= 0; i--) {
      double sm = B[i + j*p];
      for (k = i+1; k < p; k++) sm -= A[i + k*p]*B[k + j*p];
      B[i + j*p] = sm/A[i + i*p];
    }
}

static void interpmat_1d(const double *t, i8 q, const double *s, i8 p, double *L)
{
  i8 i, j;
  double *A = (double*) malloc(sizeof(double)*(size_t)(p*p));
  double *B = (double*) malloc(sizeof(double)*(size_t)(p*q));

  for (i = 0; i < p; i++) {
    A[0 + i*p] = 1.0;
    for (j = 1; j < p; j++) A[j + i*p] = A[j-1 + i*p]*s[i];
  }
  for (i = 0; i < q; i++) {
    B[0 + i*p] = 1.0;
    for (j = 1; j < p; j++) B[j + i*p] = B[j-1 + i*p]*t[i];
  }

  lu_solve(A, p, B, q);

  for (i = 0; i < q; i++)
    for (j = 0; j < p; j++)
      L[i + j*q] = B[j + i*p];

  free(B);  free(A);
}

#endif

static void build_Mt(const double *uv, i8 hdim, const double *x0, i8 p, double *Mt)
{
  i8 i, a, b;
  double *t1 = (double*) malloc(sizeof(double)*(size_t)hdim);
  double *t2 = (double*) malloc(sizeof(double)*(size_t)hdim);
  double *L1 = (double*) malloc(sizeof(double)*(size_t)(hdim*p));
  double *L2 = (double*) malloc(sizeof(double)*(size_t)(hdim*p));

  for (i = 0; i < hdim; i++) { t1[i] = uv[2*i]; t2[i] = uv[2*i + 1]; }
  interpmat_1d(t1, hdim, x0, p, L1);
  interpmat_1d(t2, hdim, x0, p, L2);

  for (i = 0; i < hdim; i++)
    for (a = 0; a < p; a++)
      for (b = 0; b < p; b++)
        Mt[(b + a*p) + i*(p*p)] = L1[i + a*hdim] * L2[i + b*hdim];

  free(L2); free(L1); free(t2); free(t1);
}

i8 stellarator_geo_ntri(const i8 *mp, const i8 *np, const i8 *p, const double *x0)
{
  return 2*build_charts(*mp, *np, *p, x0,
                        M_PI/(double)(*np), M_PI/(double)(*mp), NULL);
}

void stellarator_geo_uv2x(const i8 *mp_, const i8 *np_, const i8 *p_,
                          const double *x0, const i8 *itri_, const i8 *nuv_,
                          const double *uv, double *x)
{
  i8 mp = *mp_, np = *np_, p = *p_, itri = *itri_, nuv = *nuv_;
  double ts1 = M_PI/(double)np, ts2 = M_PI/(double)mp;
  i8 ipan = (itri - 1)/2, side = (itri - 1) % 2, i, d;
  chart_t *cl, c;
  double xll[3], xul[3], xlr[3], xur[3], d1 = 0.0, d2 = 0.0;
  double V[3][2];

  cl = (chart_t*) malloc(sizeof(chart_t)*(size_t)(2*mp*np));
  build_charts(mp, np, p, x0, ts1, ts2, cl);
  c = cl[ipan];
  free(cl);

  chart_eval(&c, ts1, ts2, x0[0],   x0[0],   xll);
  chart_eval(&c, ts1, ts2, x0[p-1], x0[p-1], xur);
  chart_eval(&c, ts1, ts2, x0[0],   x0[p-1], xul);
  chart_eval(&c, ts1, ts2, x0[p-1], x0[0],   xlr);
  for (d = 0; d < 3; d++) {
    d1 += (xll[d]-xur[d])*(xll[d]-xur[d]);
    d2 += (xul[d]-xlr[d])*(xul[d]-xlr[d]);
  }

  if (d1 + 1e-13 > d2) {
    if (side == 0) { V[0][0]=-1; V[0][1]=-1;  V[1][0]= 1; V[1][1]=-1;  V[2][0]=-1; V[2][1]= 1; }
    else           { V[0][0]= 1; V[0][1]= 1;  V[1][0]=-1; V[1][1]= 1;  V[2][0]= 1; V[2][1]=-1; }
  } else {
    if (side == 0) { V[0][0]=-1; V[0][1]= 1;  V[1][0]=-1; V[1][1]=-1;  V[2][0]= 1; V[2][1]= 1; }
    else           { V[0][0]= 1; V[0][1]=-1;  V[1][0]= 1; V[1][1]= 1;  V[2][0]=-1; V[2][1]=-1; }
  }

  for (i = 0; i < nuv; i++) {
    double u = uv[2*i], v = uv[2*i + 1];
    double t1 = V[0][0] + (V[1][0]-V[0][0])*u + (V[2][0]-V[0][0])*v;
    double t2 = V[0][1] + (V[1][1]-V[0][1])*u + (V[2][1]-V[0][1])*v;
    chart_eval(&c, ts1, ts2, t1, t2, &x[3*i]);
  }
}

void stellarator_geo(const i8 *mp_, const i8 *np_, const i8 *p_,
                     const double *x0, const double *D,
                     const double *uvs, const double *wts,
                     double *sxo, double *snxo, double *swo,
                     double *rtso, double *rpso)
{
  i8 mp = *mp_, np = *np_, p = *p_;
  i8 p2 = p*p, hdim = p*(p+1)/2;
  double ts1 = M_PI/(double)np, ts2 = M_PI/(double)mp;
  i8 nq, ipan, itri, i, k, m, a, b, d, r1, c1;
  chart_t *cl;
  double *xx1, *xx2, *sxk, *rts, *rps, *S;
  double *MtlA, *MtrA, *MtlB, *MtrB, *uvtmp;

  cl = (chart_t*) malloc(sizeof(chart_t)*(size_t)(2*mp*np));
  nq = build_charts(mp, np, p, x0, ts1, ts2, cl);

  MtlA  = (double*) malloc(sizeof(double)*(size_t)(p2*hdim));
  MtrA  = (double*) malloc(sizeof(double)*(size_t)(p2*hdim));
  MtlB  = (double*) malloc(sizeof(double)*(size_t)(p2*hdim));
  MtrB  = (double*) malloc(sizeof(double)*(size_t)(p2*hdim));
  uvtmp = (double*) malloc(sizeof(double)*(size_t)(2*hdim));

  for (i = 0; i < hdim; i++) {
    uvtmp[2*i]     = 2.0*uvs[2*i]     - 1.0;
    uvtmp[2*i + 1] = 2.0*uvs[2*i + 1] - 1.0;
  }
  build_Mt(uvtmp, hdim, x0, p, MtlA);
  for (i = 0; i < 2*hdim; i++) uvtmp[i] = -uvtmp[i];
  build_Mt(uvtmp, hdim, x0, p, MtrA);

  for (i = 0; i < hdim; i++) {
    double u = 2.0*uvs[2*i] - 1.0, v = 2.0*uvs[2*i + 1] - 1.0;
    uvtmp[2*i] = v;  uvtmp[2*i + 1] = -u;
  }
  build_Mt(uvtmp, hdim, x0, p, MtlB);
  for (i = 0; i < 2*hdim; i++) uvtmp[i] = -uvtmp[i];
  build_Mt(uvtmp, hdim, x0, p, MtrB);

  xx1 = (double*) malloc(sizeof(double)*(size_t)p2);
  xx2 = (double*) malloc(sizeof(double)*(size_t)p2);
  sxk = (double*) malloc(sizeof(double)*(size_t)(3*p2));
  rts = (double*) malloc(sizeof(double)*(size_t)(3*p2));
  rps = (double*) malloc(sizeof(double)*(size_t)(3*p2));
  S   = (double*) malloc(sizeof(double)*(size_t)(9*p2));

  for (b = 0; b < p; b++)
    for (a = 0; a < p; a++) { xx1[a + b*p] = x0[b];  xx2[a + b*p] = x0[a]; }

  itri = 0;
  for (ipan = 0; ipan < nq; ipan++) {
    double *Mtl, *Mtr, d1, d2;

    for (m = 0; m < p2; m++)
      chart_eval(&cl[ipan], ts1, ts2, xx1[m], xx2[m], &sxk[3*m]);

    for (c1 = 0; c1 < p; c1++)
      for (r1 = 0; r1 < p; r1++)
        for (d = 0; d < 3; d++) {
          double s = 0.0;
          for (k = 0; k < p; k++) s += D[r1 + k*p] * sxk[3*(k + c1*p) + d];
          rts[3*(r1 + c1*p) + d] = s;
        }
    for (r1 = 0; r1 < p; r1++)
      for (c1 = 0; c1 < p; c1++)
        for (d = 0; d < 3; d++) {
          double s = 0.0;
          for (k = 0; k < p; k++) s += D[c1 + k*p] * sxk[3*(r1 + k*p) + d];
          rps[3*(r1 + c1*p) + d] = s;
        }

    d1 = 0.0; d2 = 0.0;
    for (d = 0; d < 3; d++) {
      double u = sxk[3*0 + d]     - sxk[3*(p2-1) + d];
      double v = sxk[3*(p-1) + d] - sxk[3*(p2-p) + d];
      d1 += u*u;  d2 += v*v;
    }
    if (d1 + 1e-13 > d2) { Mtl = MtlA;  Mtr = MtrA; }
    else                 { Mtl = MtlB;  Mtr = MtrB; }

    for (m = 0; m < p2; m++)
      for (d = 0; d < 3; d++) {
        S[0 + d + m*9] = sxk[3*m + d];
        S[3 + d + m*9] = rts[3*m + d];
        S[6 + d + m*9] = rps[3*m + d];
      }

    for (k = 0; k < 2; k++) {
      const double *Mt = (k == 0) ? Mtl : Mtr;
      i8 off = (itri + k)*hdim;

      for (i = 0; i < hdim; i++) {
        double acc[9], xt[3], xp[3], nx[3], sp;
        for (d = 0; d < 9; d++) acc[d] = 0.0;
        for (m = 0; m < p2; m++) {
          double c = Mt[m + i*p2];
          for (d = 0; d < 9; d++) acc[d] += S[d + m*9]*c;
        }
        for (d = 0; d < 3; d++) { xt[d] = acc[3+d];  xp[d] = acc[6+d]; }
        nx[0] = xp[1]*xt[2] - xp[2]*xt[1];
        nx[1] = xp[2]*xt[0] - xp[0]*xt[2];
        nx[2] = xp[0]*xt[1] - xp[1]*xt[0];
        sp = sqrt(nx[0]*nx[0] + nx[1]*nx[1] + nx[2]*nx[2]);
        for (d = 0; d < 3; d++) {
          sxo [3*(off + i) + d] = acc[d];
          snxo[3*(off + i) + d] = nx[d]/sp;
          rtso[3*(off + i) + d] = 2.0*xt[d];
          rpso[3*(off + i) + d] = 2.0*xp[d];
        }
        swo[off + i] = sp*wts[i]*4.0;
      }
    }
    itri += 2;
  }

  free(S); free(rps); free(rts); free(sxk); free(xx2); free(xx1);
  free(uvtmp); free(MtrB); free(MtlB); free(MtrA); free(MtlA); free(cl);
}

#ifdef MATLAB_MEX_FILE
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  i8 mp, np, p, hdim, ntri;
  double *x0, *D, *uvs, *wts;

  if (nrhs != 8 && nrhs != 10)
    mexErrMsgIdAndTxt("stellarator_geo:nrhs", "8 inputs (geometry) or 10 (uv2x)");
  if (nlhs > 6)  mexErrMsgIdAndTxt("stellarator_geo:nlhs", "at most 6 outputs");

  mp = (i8) mxGetScalar(prhs[0]);
  np = (i8) mxGetScalar(prhs[1]);
  p  = (i8) mxGetScalar(prhs[2]);
  x0  = mxGetDoubles(prhs[3]);
  D   = mxGetDoubles(prhs[5]);
  uvs = mxGetDoubles(prhs[6]);
  wts = mxGetDoubles(prhs[7]);

  hdim = p*(p+1)/2;

  if (nrhs == 10) {
    i8 itri = (i8) mxGetScalar(prhs[8]);
    i8 nuv  = (i8) mxGetNumberOfElements(prhs[9]) / 2;
    plhs[0] = mxCreateDoubleMatrix(3, (mwSize)nuv, mxREAL);
    stellarator_geo_uv2x(&mp, &np, &p, x0, &itri, &nuv,
                         mxGetDoubles(prhs[9]), mxGetDoubles(plhs[0]));
    return;
  }

  if ((i8) mxGetNumberOfElements(prhs[6]) != 2*hdim)
    mexErrMsgIdAndTxt("stellarator_geo:uvs", "uvs must be 2-by-%d", (int)hdim);

  ntri = stellarator_geo_ntri(&mp, &np, &p, x0);

  plhs[0] = mxCreateDoubleMatrix(3, (mwSize)(ntri*hdim), mxREAL);
  plhs[1] = mxCreateDoubleMatrix(3, (mwSize)(ntri*hdim), mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1, (mwSize)(ntri*hdim), mxREAL);
  plhs[3] = mxCreateDoubleMatrix(3, (mwSize)(ntri*hdim), mxREAL);
  plhs[4] = mxCreateDoubleMatrix(3, (mwSize)(ntri*hdim), mxREAL);
  if (nlhs >= 6) plhs[5] = mxCreateDoubleScalar((double)ntri);

  stellarator_geo(&mp, &np, &p, x0, D, uvs, wts,
                  mxGetDoubles(plhs[0]), mxGetDoubles(plhs[1]),
                  mxGetDoubles(plhs[2]), mxGetDoubles(plhs[3]),
                  mxGetDoubles(plhs[4]));
}
#endif

