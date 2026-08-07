/*
 * stellarator_geo_mex.c
 *
 * One call for the whole geometry build: create_panels -> stellarator_pan_init
 * -> the get_high_order_quad / get_vioreanu_quadr loop -> the cellfun/horzcat.
 *
 * The compute core is plain C with a bind(C) signature, so one file serves the
 * MATLAB gateway and Fortran callers.  No BLAS or LAPACK: the Vandermonde
 * solve is p<=16 and the per-panel product is small, so both are local.  That
 * keeps libmwlapack out of the mex build, OpenBLAS out of the Fortran build,
 * and makes the two paths agree bit for bit.
 *
 * MATLAB:
 *
 *   mex -R2018a stellarator_geo_mex.c
 *
 *   [x0,w0,D] = gauss(order);
 *   [uvs,wts] = get_vioreanu_nodes(order-1);
 *   [sx,snx,sw,rts,rps,npan] = stellarator_geo_mex(mp,np,order,x0,w0,D,uvs,wts);
 *
 * Fortran, see test/stellarator_grf_debug.f90:
 *
 *   cc -O3 -c stellarator_geo_mex.c
 *
 *   ntri = stellarator_geo_ntri(mp, np, order, x0)
 *   call stellarator_geo(mp, np, order, x0, D, uvs, wts, sx, snx, sw, rts, rps)
 *
 * Outputs are concatenated over patches, two triangles per quad, in the order
 * pan{2*(k-1)+1}=sl, pan{2*(k-1)+2}=sr:
 *
 *   sx, snx, rts, rps   (3, ntri*hdim)      hdim = order*(order+1)/2
 *   sw                  (1, ntri*hdim)
 *
 * rts/rps already carry the factor 2 that the MATLAB applies as
 * sl.rts = 2*sl.xpt, so they drop straight into the source arrays.
 *
 * Second derivatives and curvature are not computed: nothing downstream of
 * this block reads .xptt/.xpss/.xpts/.hx/.cur/.kcur.
 *
 * get_high_order_quad is a no-op for side='e' -- it returns sl.x = pan.x and
 * everything else it computes is discarded by the caller -- so only the
 * side='i' index flip would need adding if an interior run is ever wanted.
 *
 * PSPLIT_MATCH_MATLAB reproduces a bug the .m used to have: its dist_p/dist_t
 * branch gave both halves the same chart (t(2,:)/2-1/2 twice), so a p-split
 * yielded two copies of the lower half and dropped the upper half.  Both
 * copies of stellarator_pan_init are fixed now, and the branch never fires for
 * the six standard cases (0 p-splits, t-splits only), so this is dead history
 * kept only for reproducing old meshes.
 */

#define PSPLIT_MATCH_MATLAB 0

#include <math.h>
#include <stdlib.h>

#define ASPECT 1.65

typedef long long i8;

/* the 7-term Fourier surface, utils/torus/stellaratorparam.m */
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

/* A boundary Fourier series in the VMEC convention,
     R(theta,phi) = sum Rmn cos(m*theta - n*nfp*phi)
     Z(theta,phi) = sum Zmn sin(m*theta - n*nfp*phi)
   with theta poloidal and phi toroidal.  nmode = 0 selects the built-in
   7-term test surface above, so every existing caller keeps its bits. */
typedef struct {
  i8 nfp, nmode;
  const double *m, *n, *rc, *zs;
} surf_t;

static const surf_t SURF_DEFAULT = { 1, 0, 0, 0, 0, 0 };

static void surf_eval(const surf_t *s, double t, double p, double *x)
{
  double R = 0.0, Z = 0.0, ang;
  i8 k;
  if (s == 0 || s->nmode <= 0) { stellaratorparam(t, p, x);  return; }
  for (k = 0; k < s->nmode; k++) {
    ang = s->m[k]*p - (double)s->nfp*s->n[k]*t;
    R += s->rc[k]*cos(ang);
    Z += s->zs[k]*sin(ang);
  }
  x[0] = R*cos(t);  x[1] = R*sin(t);  x[2] = Z;
}

typedef struct {
  double a1, b1, a2, b2;   /* t1 -> a1*t1+b1,  t2 -> a2*t2+b2 */
  double jj, ii;           /* panel indices, 1-based as in create_panels */
} chart_t;

static void chart_eval_s(const chart_t *c, double ts1, double ts2,
                         double t1, double t2, const surf_t *sf, double *x)
{
  double u = ((c->a1*t1 + c->b1) + 2.0*c->jj) * ts1;   /* orig = 1 */
  double v = ((c->a2*t2 + c->b2) + 2.0*c->ii) * ts2;
  surf_eval(sf, u, v, x);
}

static void chart_eval(const chart_t *c, double ts1, double ts2,
                       double t1, double t2, double *x)
{
  chart_eval_s(c, ts1, ts2, t1, t2, &SURF_DEFAULT, x);
}

/* pass 1: aspect test on three corners, fill the chart list, return the count.
   cl may be NULL to count only. */
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
      chart_eval(&c0, ts1, ts2, x0[0],   x0[0],   xll);   /* x(:,1) */
      chart_eval(&c0, ts1, ts2, x0[0],   x0[p-1], xul);   /* x(:,order) */
      chart_eval(&c0, ts1, ts2, x0[p-1], x0[0],   xlr);   /* x(:,end-order+1) */
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

/* interpmat_1d: the interpolant of interpmat_1d.m, evaluated in the second
   barycentric form instead of through its monomial Vandermonde.  Same
   polynomial, but the monomial solve at p=14 is conditioned around 1e5 and
   turns ulp into ~1e-13 in the node positions, which is the scale of the
   patch-to-patch seam mismatch that limits close evaluation near an edge.
   Set BARYCENTRIC to 0 to get the original monomial solve back.  L[i + j*q]. */
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

/* A X = B in place, A p-by-p, B p-by-q, column major, partial pivoting */
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

/* Mt[m + i*p2] = L1[i][a]*L2[i][b], m = b + a*p */
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

/* ================= adaptive chart list =================================
   The refinement w7x.m does, in C.  chart_t already represents any axis
   aligned affine sub-rectangle of a base cell, so a 4-way split is just
   a1 -> a1/2, b1 -> b1 +- a1/2 (and the same in 2): nothing new is stored,
   and the parameter box a chart covers follows from its own coefficients.

     1  iterate the aspect split until every panel is within ASPECT.  One
        pass is enough for the test torus; W7-X starts near 3:1.
     2  refine on curvature: expand the mean curvature H(t1,t2) of the panel
        in the tensor Legendre basis on its own p x p grid and split while the
        last two degrees carry more than restol of the whole spectrum.  H, not
        position: a panel's coordinates are dominated by where it sits and how
        big it is.  Panels that pass are marked and never rechecked.
     3  2:1 balance, so no two face neighbours differ by more than a factor of
        two along the face they share.

   restol <= 0 skips 2 and 3, which reproduces the old single-pass behaviour.
   ====================================================================== */

#define MAXASPECT_PASS 8
#define MAXCURV_LEVEL 12
#define MAXBAL_SWEEP  20

static void chart_box(const chart_t *c, double ts1, double ts2, double *bx)
{
  bx[0] = ((c->b1 - c->a1) + 2.0*c->jj)*ts1;
  bx[1] = ((c->b1 + c->a1) + 2.0*c->jj)*ts1;
  bx[2] = ((c->b2 - c->a2) + 2.0*c->ii)*ts2;
  bx[3] = ((c->b2 + c->a2) + 2.0*c->ii)*ts2;
}

/* the panel on the p x p tensor grid: node m = a + b*p is (t1,t2) = (x0[b],x0[a]) */
static void chart_grid(const chart_t *c, double ts1, double ts2, i8 p,
                       const double *x0, const surf_t *sf, double *sxk)
{
  i8 a, b;
  for (b = 0; b < p; b++)
    for (a = 0; a < p; a++)
      chart_eval_s(c, ts1, ts2, x0[b], x0[a], sf, &sxk[3*(a + b*p)]);
}

/* d/d(index a) of a 3 x p*p field, i.e. along t2 */
static void diff_a(const double *D, i8 p, const double *f, double *g)
{
  i8 a, b, k, d;
  for (b = 0; b < p; b++)
    for (a = 0; a < p; a++)
      for (d = 0; d < 3; d++) {
        double sm = 0.0;
        for (k = 0; k < p; k++) sm += D[a + k*p]*f[3*(k + b*p) + d];
        g[3*(a + b*p) + d] = sm;
      }
}

/* d/d(index b), i.e. along t1 */
static void diff_b(const double *D, i8 p, const double *f, double *g)
{
  i8 a, b, k, d;
  for (b = 0; b < p; b++)
    for (a = 0; a < p; a++)
      for (d = 0; d < 3; d++) {
        double sm = 0.0;
        for (k = 0; k < p; k++) sm += D[b + k*p]*f[3*(a + k*p) + d];
        g[3*(a + b*p) + d] = sm;
      }
}

static void legvander(i8 p, const double *x0, double *V)
{
  i8 i, k;
  for (i = 0; i < p; i++) {
    V[i + 0*p] = 1.0;
    if (p > 1) V[i + 1*p] = x0[i];
    for (k = 1; k < p-1; k++)
      V[i + (k+1)*p] = ((2.0*k + 1.0)*x0[i]*V[i + k*p] - k*V[i + (k-1)*p])/(k + 1.0);
  }
}

/* relative weight of the last two degrees in the Legendre spectrum of H */
static double panel_resolution(const chart_t *c, double ts1, double ts2, i8 p,
                               const double *x0, const double *D,
                               const surf_t *sf, const double *V, double *wk)
{
  i8 n2 = p*p, a, b, k, d;
  double *sxk = wk;                 /* 3*n2 */
  double *r1  = sxk + 3*n2, *r2 = r1 + 3*n2;
  double *r11 = r2 + 3*n2, *r12 = r11 + 3*n2, *r22 = r12 + 3*n2;
  double *H   = r22 + 3*n2;         /* n2   */
  double *C   = H + n2, *T = C + n2, *A = T + n2;   /* n2, n2, p*p */
  double tail = 0.0, tot = 0.0;

  chart_grid(c, ts1, ts2, p, x0, sf, sxk);
  diff_b(D, p, sxk, r1);   diff_a(D, p, sxk, r2);
  diff_b(D, p, r1,  r11);  diff_b(D, p, r2,  r12);  diff_a(D, p, r2, r22);

  for (k = 0; k < n2; k++) {
    double *u = &r1[3*k], *v = &r2[3*k];
    double nv[3], nn, E, F, G, L, M, N, den;
    nv[0] = u[1]*v[2] - u[2]*v[1];
    nv[1] = u[2]*v[0] - u[0]*v[2];
    nv[2] = u[0]*v[1] - u[1]*v[0];
    nn = sqrt(nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2]);
    if (nn <= 0.0) nn = 1.0;
    nv[0] /= nn;  nv[1] /= nn;  nv[2] /= nn;
    E = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    F = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
    G = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    L = nv[0]*r11[3*k] + nv[1]*r11[3*k+1] + nv[2]*r11[3*k+2];
    M = nv[0]*r12[3*k] + nv[1]*r12[3*k+1] + nv[2]*r12[3*k+2];
    N = nv[0]*r22[3*k] + nv[1]*r22[3*k+1] + nv[2]*r22[3*k+2];
    den = 2.0*(E*G - F*F);
    if (den == 0.0) den = 1.0;
    H[k] = (E*N - 2.0*F*M + G*L)/den;
  }

  /* C = V \ H / V.'  -- rows (t2) then columns (t1) */
  for (k = 0; k < n2; k++) C[k] = H[k];
  for (k = 0; k < p*p; k++) A[k] = V[k];
  lu_solve(A, p, C, p);
  for (a = 0; a < p; a++) for (b = 0; b < p; b++) T[b + a*p] = C[a + b*p];
  for (k = 0; k < p*p; k++) A[k] = V[k];
  lu_solve(A, p, T, p);
  for (a = 0; a < p; a++) for (b = 0; b < p; b++) C[a + b*p] = T[b + a*p];

  for (a = 0; a < p; a++)
    for (b = 0; b < p; b++) {
      double z = C[a + b*p];
      tot += z*z;
      if (a >= p-2 || b >= p-2) tail += z*z;
    }
  return sqrt(tail/(tot > 0.0 ? tot : 1.0));
}

static void split4(const chart_t *c, chart_t *out)
{
  i8 q;
  for (q = 0; q < 4; q++) {
    out[q]    = *c;
    out[q].a1 = c->a1/2.0;
    out[q].a2 = c->a2/2.0;
    out[q].b1 = c->b1 + ((q & 1) ? c->a1/2.0 : -c->a1/2.0);
    out[q].b2 = c->b2 + ((q & 2) ? c->a2/2.0 : -c->a2/2.0);
  }
}

/* Build the chart list.  Returns the count, or -1 if it would exceed cap.
   cl may be NULL to count only (the refinement then runs on a scratch list). */
static i8 build_charts_adaptive(i8 mp, i8 np, i8 p, const double *x0,
                                const double *D, double ts1, double ts2,
                                const surf_t *sf, double restol, i8 cap,
                                chart_t *cl)
{
  chart_t *L, kid[4];
  char *done;
  double *wk, *V, bx[4];
  i8 n, i, j, k, q, pass, nsplit, ret = -1;

  if (cap <= 0) cap = 200000;
  L    = (chart_t*) malloc(sizeof(chart_t)*(size_t)cap);
  done = (char*)    malloc((size_t)cap);
  V    = (double*)  malloc(sizeof(double)*(size_t)(p*p));
  wk   = (double*)  malloc(sizeof(double)*(size_t)(18*p*p + 3*p*p + p*p));
  if (!L || !done || !V || !wk) goto out;
  legvander(p, x0, V);

  n = 0;
  for (j = 1; j <= np; j++)
    for (i = 1; i <= mp; i++) {
      if (n >= cap) goto out;
      L[n].a1 = 1.0; L[n].b1 = 0.0; L[n].a2 = 1.0; L[n].b2 = 0.0;
      L[n].jj = (double)j;  L[n].ii = (double)i;  n++;
    }

  /* 1. aspect */
  for (pass = 0; pass < MAXASPECT_PASS; pass++) {
    nsplit = 0;
    for (k = 0; k < n; k++) {
      double xll[3], xul[3], xlr[3], dt, dp;
      chart_t c = L[k];
      chart_eval_s(&c, ts1, ts2, x0[0],   x0[0],   sf, xll);
      chart_eval_s(&c, ts1, ts2, x0[0],   x0[p-1], sf, xul);
      chart_eval_s(&c, ts1, ts2, x0[p-1], x0[0],   sf, xlr);
      dt = sqrt((xll[0]-xlr[0])*(xll[0]-xlr[0]) + (xll[1]-xlr[1])*(xll[1]-xlr[1])
              + (xll[2]-xlr[2])*(xll[2]-xlr[2]));
      dp = sqrt((xll[0]-xul[0])*(xll[0]-xul[0]) + (xll[1]-xul[1])*(xll[1]-xul[1])
              + (xll[2]-xul[2])*(xll[2]-xul[2]));
      if (dt/dp > ASPECT) {
        if (n >= cap) goto out;
        L[k].a1 = c.a1/2.0;  L[k].b1 = c.b1 - c.a1/2.0;
        L[n] = c;  L[n].a1 = c.a1/2.0;  L[n].b1 = c.b1 + c.a1/2.0;  n++;
        nsplit++;
      } else if (dp/dt > ASPECT) {
        if (n >= cap) goto out;
        L[k].a2 = c.a2/2.0;  L[k].b2 = c.b2 - c.a2/2.0;
        L[n] = c;  L[n].a2 = c.a2/2.0;  L[n].b2 = c.b2 + c.a2/2.0;  n++;
        nsplit++;
      }
    }
    if (!nsplit) break;
  }

  /* 2. curvature */
  if (restol > 0.0) {
    for (k = 0; k < n; k++) done[k] = 0;
    for (pass = 0; pass < MAXCURV_LEVEL; pass++) {
      i8 n0 = n;
      nsplit = 0;
      for (k = 0; k < n0; k++) {
        if (done[k]) continue;
        if (panel_resolution(&L[k], ts1, ts2, p, x0, D, sf, V, wk) <= restol) {
          done[k] = 1;  continue;
        }
        if (n + 3 > cap) goto out;
        split4(&L[k], kid);
        L[k] = kid[0];  done[k] = 0;
        for (q = 1; q < 4; q++) { L[n] = kid[q];  done[n] = 0;  n++; }
        nsplit++;
      }
      if (!nsplit) break;
    }

    /* 3. 2:1 balance */
    for (pass = 0; pass < MAXBAL_SWEEP; pass++) {
      i8 n0 = n;
      nsplit = 0;
      for (k = 0; k < n0; k++) done[k] = 0;
      for (k = 0; k < n0; k++) {
        double bk[4];
        chart_box(&L[k], ts1, ts2, bk);
        for (q = 0; q < n0; q++) {
          double bq[4], du, dv, ek, eq;
          if (q == k) continue;
          if (fabs(L[q].jj - L[k].jj) > 1.5 && fabs(fabs(L[q].jj - L[k].jj) - np) > 1.5) continue;
          if (fabs(L[q].ii - L[k].ii) > 1.5 && fabs(fabs(L[q].ii - L[k].ii) - mp) > 1.5) continue;
          chart_box(&L[q], ts1, ts2, bq);
          du = 1e-9*ts1;  dv = 1e-9*ts2;
          if ((fabs(bk[1]-bq[0]) < du || fabs(bk[0]-bq[1]) < du) &&
              fmin(bk[3],bq[3]) - fmax(bk[2],bq[2]) > dv) {
            ek = bk[3]-bk[2];  eq = bq[3]-bq[2];
            if (ek > 2.0*eq + dv) done[k] = 1;
          } else if ((fabs(bk[3]-bq[2]) < dv || fabs(bk[2]-bq[3]) < dv) &&
                     fmin(bk[1],bq[1]) - fmax(bk[0],bq[0]) > du) {
            ek = bk[1]-bk[0];  eq = bq[1]-bq[0];
            if (ek > 2.0*eq + du) done[k] = 1;
          }
        }
      }
      for (k = 0; k < n0; k++) {
        if (!done[k]) continue;
        if (n + 3 > cap) goto out;
        split4(&L[k], kid);
        L[k] = kid[0];
        for (q = 1; q < 4; q++) { L[n] = kid[q];  n++; }
        nsplit++;
      }
      if (!nsplit) break;
    }
  }

  if (cl) for (k = 0; k < n; k++) cl[k] = L[k];
  ret = n;

out:
  free(wk);  free(V);  free(done);  free(L);
  return ret;
}

/* ---- exported: triangle count, so a Fortran caller can allocate first ---- */
i8 stellarator_geo_ntri(const i8 *mp, const i8 *np, const i8 *p, const double *x0)
{
  return 2*build_charts(*mp, *np, *p, x0,
                        M_PI/(double)(*np), M_PI/(double)(*mp), NULL);
}

/* the (u,v) -> x map for one triangle of one chart; side comes from itri */
static void uv2x_core(const chart_t *c, i8 side, double ts1, double ts2, i8 p,
                      const double *x0, const surf_t *sf,
                      i8 nuv, const double *uv, double *x)
{
  double xll[3], xul[3], xlr[3], xur[3], d1 = 0.0, d2 = 0.0, V[3][2];
  i8 i, d;
  chart_eval_s(c, ts1, ts2, x0[0],   x0[0],   sf, xll);
  chart_eval_s(c, ts1, ts2, x0[p-1], x0[p-1], sf, xur);
  chart_eval_s(c, ts1, ts2, x0[0],   x0[p-1], sf, xul);
  chart_eval_s(c, ts1, ts2, x0[p-1], x0[0],   sf, xlr);
  for (d = 0; d < 3; d++) {
    d1 += (xll[d]-xur[d])*(xll[d]-xur[d]);
    d2 += (xul[d]-xlr[d])*(xul[d]-xlr[d]);
  }
  if (d1 + 1e-13 > d2) {
    if (side == 0) { V[0][0]=-1; V[0][1]=-1;  V[1][0]= 1; V[1][1]=-1;  V[2][0]=-1; V[2][1]= 1; }
    else                { V[0][0]= 1; V[0][1]= 1;  V[1][0]=-1; V[1][1]= 1;  V[2][0]= 1; V[2][1]=-1; }
  } else {
    if (side == 0) { V[0][0]=-1; V[0][1]= 1;  V[1][0]=-1; V[1][1]=-1;  V[2][0]= 1; V[2][1]= 1; }
    else                { V[0][0]= 1; V[0][1]=-1;  V[1][0]= 1; V[1][1]= 1;  V[2][0]=-1; V[2][1]=-1; }
  }
  for (i = 0; i < nuv; i++) {
    double u = uv[2*i], v = uv[2*i + 1];
    double t1 = V[0][0] + (V[1][0]-V[0][0])*u + (V[2][0]-V[0][0])*v;
    double t2 = V[0][1] + (V[1][1]-V[0][1])*u + (V[2][1]-V[0][1])*v;
    chart_eval_s(c, ts1, ts2, t1, t2, sf, &x[3*i]);
  }
}

/* ---- exported: reference-triangle (u,v) of one triangle -> chart point.
   itri is 1-based over the triangle list, so panel = (itri-1)/2 and
   side = (itri-1)%2 with 0 = sl, 1 = sr.  The triangle's vertices in quad
   parameters are exactly (+-1,+-1); which three, and in which order, is the
   same diagonal choice get_vioreanu_quadr makes.  Evaluating the chart here
   rather than a per-patch Koornwinder fit is what makes two patches agree on
   the edge they share: both reduce to the same floating point expression
   (t + 2*j)*ts, so the shared parameter line is bit identical. ---- */
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
  uv2x_core(&c, side, ts1, ts2, p, x0, &SURF_DEFAULT, nuv, uv, x);
  return;

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

/* ---- exported: the geometry ---- */
static void geo_core(const chart_t *cl, i8 nq, i8 mp, i8 np, i8 p,
                     const double *x0, const double *D,
                     const double *uvs, const double *wts, const surf_t *sf,
                     double *sxo, double *snxo, double *swo,
                     double *rtso, double *rpso)
{
  i8 p2 = p*p, hdim = p*(p+1)/2;
  double ts1 = M_PI/(double)np, ts2 = M_PI/(double)mp;
  i8 ipan, itri, i, k, m, a, b, d, r1, c1;
  double *xx1, *xx2, *sxk, *rts, *rps, *S;
  double *MtlA, *MtrA, *MtlB, *MtrB, *uvtmp;

  MtlA  = (double*) malloc(sizeof(double)*(size_t)(p2*hdim));
  MtrA  = (double*) malloc(sizeof(double)*(size_t)(p2*hdim));
  MtlB  = (double*) malloc(sizeof(double)*(size_t)(p2*hdim));
  MtrB  = (double*) malloc(sizeof(double)*(size_t)(p2*hdim));
  uvtmp = (double*) malloc(sizeof(double)*(size_t)(2*hdim));

  for (i = 0; i < hdim; i++) {                 /* uvsl = 2*uvs - 1 */
    uvtmp[2*i]     = 2.0*uvs[2*i]     - 1.0;
    uvtmp[2*i + 1] = 2.0*uvs[2*i + 1] - 1.0;
  }
  build_Mt(uvtmp, hdim, x0, p, MtlA);
  for (i = 0; i < 2*hdim; i++) uvtmp[i] = -uvtmp[i];       /* uvsr = -uvsl */
  build_Mt(uvtmp, hdim, x0, p, MtrA);

  for (i = 0; i < hdim; i++) {                 /* rotate by exp(-i*pi/2) */
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

  /* meshgrid(x0): node m = a + b*p has xx1 = x0[b], xx2 = x0[a] */
  for (b = 0; b < p; b++)
    for (a = 0; a < p; a++) { xx1[a + b*p] = x0[b];  xx2[a + b*p] = x0[a]; }

  itri = 0;
  for (ipan = 0; ipan < nq; ipan++) {
    double *Mtl, *Mtr, d1, d2;

    for (m = 0; m < p2; m++)
      chart_eval_s(&cl[ipan], ts1, ts2, xx1[m], xx2[m], sf, &sxk[3*m]);

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
      double u = sxk[3*0 + d]     - sxk[3*(p2-1) + d];    /* rll - rur */
      double v = sxk[3*(p-1) + d] - sxk[3*(p2-p) + d];    /* rul - rlr */
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
  free(uvtmp); free(MtrB); free(MtlB); free(MtrA); free(MtlA);
}

/* ---- exported: the geometry, unchanged interface and bits ---- */
void stellarator_geo(const i8 *mp_, const i8 *np_, const i8 *p_,
                     const double *x0, const double *D,
                     const double *uvs, const double *wts,
                     double *sxo, double *snxo, double *swo,
                     double *rtso, double *rpso)
{
  i8 mp = *mp_, np = *np_, p = *p_, nq;
  double ts1 = M_PI/(double)np, ts2 = M_PI/(double)mp;
  chart_t *cl = (chart_t*) malloc(sizeof(chart_t)*(size_t)(2*mp*np));
  nq = build_charts(mp, np, p, x0, ts1, ts2, cl);
  geo_core(cl, nq, mp, np, p, x0, D, uvs, wts, &SURF_DEFAULT,
           sxo, snxo, swo, rtso, rpso);
  free(cl);
}

/* ================= exported: the adaptive path =========================
   Build once, query many.  The chart list is the caller's buffer -- 6 doubles
   per chart, {a1,b1,a2,b2,jj,ii} -- so nothing is cached between calls and the
   routines stay safe to call from several threads.  The surface is passed in:
   nmode = 0 keeps the built-in test torus, otherwise mn is 2*nmode doubles
   holding (m,n) pairs and rc/zs are nmode each.
   ====================================================================== */

i8 stellarator_charts(const i8 *mp_, const i8 *np_, const i8 *p_,
                      const double *x0, const double *D,
                      const i8 *nfp_, const i8 *nmode_, const double *mn,
                      const double *rc, const double *zs,
                      const double *restol_, const i8 *cap_, double *cl)
{
  i8 mp = *mp_, np = *np_, p = *p_, nmode = *nmode_, cap = *cap_, n, k;
  double ts1 = M_PI/(double)np, ts2 = M_PI/(double)mp;
  double *mm = 0, *nn = 0;
  surf_t sf;
  chart_t *tmp;

  sf = SURF_DEFAULT;
  if (nmode > 0) {
    mm = (double*) malloc(sizeof(double)*(size_t)(2*nmode));
    if (!mm) return -1;
    nn = mm + nmode;
    for (k = 0; k < nmode; k++) { mm[k] = mn[2*k];  nn[k] = mn[2*k + 1]; }
    sf.nfp = *nfp_;  sf.nmode = nmode;  sf.m = mm;  sf.n = nn;
    sf.rc = rc;      sf.zs = zs;
  }

  tmp = cl ? (chart_t*) malloc(sizeof(chart_t)*(size_t)cap) : 0;
  if (cl && !tmp) { free(mm);  return -1; }
  n = build_charts_adaptive(mp, np, p, x0, D, ts1, ts2, &sf, *restol_, cap, tmp);
  if (cl && n > 0)
    for (k = 0; k < n; k++) {
      cl[6*k+0] = tmp[k].a1;  cl[6*k+1] = tmp[k].b1;
      cl[6*k+2] = tmp[k].a2;  cl[6*k+3] = tmp[k].b2;
      cl[6*k+4] = tmp[k].jj;  cl[6*k+5] = tmp[k].ii;
    }
  free(tmp);  free(mm);
  return n;
}

static void charts_unpack(const double *cl, i8 n, chart_t *out)
{
  i8 k;
  for (k = 0; k < n; k++) {
    out[k].a1 = cl[6*k+0];  out[k].b1 = cl[6*k+1];
    out[k].a2 = cl[6*k+2];  out[k].b2 = cl[6*k+3];
    out[k].jj = cl[6*k+4];  out[k].ii = cl[6*k+5];
  }
}

static void surf_pack(const i8 *nfp_, i8 nmode, const double *mn,
                      const double *rc, const double *zs,
                      surf_t *sf, double **scratch)
{
  i8 k;
  *scratch = 0;
  *sf = SURF_DEFAULT;
  if (nmode > 0) {
    double *mm = (double*) malloc(sizeof(double)*(size_t)(2*nmode));
    if (!mm) return;
    for (k = 0; k < nmode; k++) { mm[k] = mn[2*k];  mm[nmode + k] = mn[2*k + 1]; }
    sf->nfp = *nfp_;  sf->nmode = nmode;  sf->m = mm;  sf->n = mm + nmode;
    sf->rc = rc;      sf->zs = zs;
    *scratch = mm;
  }
}

void stellarator_geo_charts(const double *cl, const i8 *nchart,
                            const i8 *mp_, const i8 *np_, const i8 *p_,
                            const double *x0, const double *D,
                            const double *uvs, const double *wts,
                            const i8 *nfp_, const i8 *nmode_, const double *mn,
                            const double *rc, const double *zs,
                            double *sxo, double *snxo, double *swo,
                            double *rtso, double *rpso)
{
  i8 n = *nchart;
  chart_t *cc = (chart_t*) malloc(sizeof(chart_t)*(size_t)n);
  double *scratch;
  surf_t sf;
  if (!cc) return;
  charts_unpack(cl, n, cc);
  surf_pack(nfp_, *nmode_, mn, rc, zs, &sf, &scratch);
  geo_core(cc, n, *mp_, *np_, *p_, x0, D, uvs, wts, &sf,
           sxo, snxo, swo, rtso, rpso);
  free(scratch);  free(cc);
}

void stellarator_uv2x_charts(const double *cl, const i8 *nchart,
                             const i8 *mp_, const i8 *np_, const i8 *p_,
                             const double *x0, const i8 *itri_, const i8 *nuv_,
                             const double *uv,
                             const i8 *nfp_, const i8 *nmode_, const double *mn,
                             const double *rc, const double *zs, double *x)
{
  i8 n = *nchart, itri = *itri_, ipan = (itri - 1)/2, side = (itri - 1) % 2;
  double ts1 = M_PI/(double)(*np_), ts2 = M_PI/(double)(*mp_);
  double *scratch;
  surf_t sf;
  chart_t c;
  if (ipan < 0 || ipan >= n) return;
  c.a1 = cl[6*ipan+0];  c.b1 = cl[6*ipan+1];
  c.a2 = cl[6*ipan+2];  c.b2 = cl[6*ipan+3];
  c.jj = cl[6*ipan+4];  c.ii = cl[6*ipan+5];
  surf_pack(nfp_, *nmode_, mn, rc, zs, &sf, &scratch);
  uv2x_core(&c, side, ts1, ts2, *p_, x0, &sf, *nuv_, uv, x);
  free(scratch);
}

/* ================= MATLAB gateway ================= */
#ifdef MATLAB_MEX_FILE
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  i8 mp, np, p, hdim, ntri;
  double *x0, *D, *uvs, *wts;

  if (nrhs != 8 && nrhs != 10 && nrhs != 11)
    mexErrMsgIdAndTxt("stellarator_geo:nrhs",
                      "8 inputs (geometry), 10 (uv2x), or 11 (adaptive geometry)");
  if (nlhs > 6)  mexErrMsgIdAndTxt("stellarator_geo:nlhs", "at most 6 outputs");

  mp = (i8) mxGetScalar(prhs[0]);
  np = (i8) mxGetScalar(prhs[1]);
  p  = (i8) mxGetScalar(prhs[2]);
  x0  = mxGetDoubles(prhs[3]);
  D   = mxGetDoubles(prhs[5]);
  uvs = mxGetDoubles(prhs[6]);
  wts = mxGetDoubles(prhs[7]);

  hdim = p*(p+1)/2;

  if (nrhs == 10) {                     /* x = f(mp,np,order,x0,w0,D,uvs,wts,itri,uv) */
    i8 itri = (i8) mxGetScalar(prhs[8]);
    i8 nuv  = (i8) mxGetNumberOfElements(prhs[9]) / 2;
    plhs[0] = mxCreateDoubleMatrix(3, (mwSize)nuv, mxREAL);
    stellarator_geo_uv2x(&mp, &np, &p, x0, &itri, &nuv,
                         mxGetDoubles(prhs[9]), mxGetDoubles(plhs[0]));
    return;
  }

  if ((i8) mxGetNumberOfElements(prhs[6]) != 2*hdim)
    mexErrMsgIdAndTxt("stellarator_geo:uvs", "uvs must be 2-by-%d", (int)hdim);

  if (nrhs == 11) {   /* [...] = f(mp,np,order,x0,w0,D,uvs,wts, surf, restol, cap)
                         surf is a 4-by-nmode matrix, rows m, n, Rmn, Zmn, or []
                         for the built-in test torus. */
    const mxArray *S = prhs[8];
    i8 nmode = (i8) (mxGetNumberOfElements(S)/4), nfp = 5, cap, nchart, k;
    double restol = mxGetScalar(prhs[9]);
    double *sd = mxGetDoubles(S), *mn = 0, *rc = 0, *zs = 0, *cl;
    mxArray *clm;

    cap = (i8) mxGetScalar(prhs[10]);
    if (mxGetM(S) != 4 && nmode > 0)
      mexErrMsgIdAndTxt("stellarator_geo:surf", "surf must be 4-by-nmode: m; n; Rmn; Zmn");
    if (nmode > 0) {
      mn = (double*) mxMalloc(sizeof(double)*(size_t)(2*nmode));
      rc = (double*) mxMalloc(sizeof(double)*(size_t)nmode);
      zs = (double*) mxMalloc(sizeof(double)*(size_t)nmode);
      for (k = 0; k < nmode; k++) {
        mn[2*k] = sd[4*k];  mn[2*k+1] = sd[4*k+1];
        rc[k]   = sd[4*k+2];  zs[k]    = sd[4*k+3];
      }
    }
    clm = mxCreateDoubleMatrix(6, (mwSize)cap, mxREAL);
    cl  = mxGetDoubles(clm);
    nchart = stellarator_charts(&mp, &np, &p, x0, D, &nfp, &nmode, mn, rc, zs,
                                &restol, &cap, cl);
    if (nchart <= 0) {
      mxDestroyArray(clm);  mxFree(zs); mxFree(rc); mxFree(mn);
      mexErrMsgIdAndTxt("stellarator_geo:charts",
                        "chart build failed or exceeded cap = %d", (int)cap);
    }
    ntri = 2*nchart;
    plhs[0] = mxCreateDoubleMatrix(3, (mwSize)(ntri*hdim), mxREAL);
    plhs[1] = mxCreateDoubleMatrix(3, (mwSize)(ntri*hdim), mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, (mwSize)(ntri*hdim), mxREAL);
    plhs[3] = mxCreateDoubleMatrix(3, (mwSize)(ntri*hdim), mxREAL);
    plhs[4] = mxCreateDoubleMatrix(3, (mwSize)(ntri*hdim), mxREAL);
    if (nlhs >= 6) plhs[5] = mxCreateDoubleScalar((double)ntri);
    stellarator_geo_charts(cl, &nchart, &mp, &np, &p, x0, D, uvs, wts,
                           &nfp, &nmode, mn, rc, zs,
                           mxGetDoubles(plhs[0]), mxGetDoubles(plhs[1]),
                           mxGetDoubles(plhs[2]), mxGetDoubles(plhs[3]),
                           mxGetDoubles(plhs[4]));
    mxDestroyArray(clm);  mxFree(zs); mxFree(rc); mxFree(mn);
    return;
  }

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
