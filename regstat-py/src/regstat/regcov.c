/* regcov.c : shared pure-C back-end for the exact covariance-change test.
 * Bound by both R (.C) and Python (ctypes). Standard C and libm only. */
#include "regcov.h"
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <complex.h>

/* ---- self-contained RNG: splitmix64 seed -> xoshiro256**, Box-Muller normals ---- */
static uint64_t sm_x;
static uint64_t splitmix64(void) {
  uint64_t z = (sm_x += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
static uint64_t xs[4];
static inline uint64_t rotl(uint64_t x, int k){ return (x << k) | (x >> (64 - k)); }
static uint64_t xoshiro(void) {
  uint64_t r = rotl(xs[1] * 5, 7) * 9;
  uint64_t t = xs[1] << 17;
  xs[2] ^= xs[0]; xs[3] ^= xs[1]; xs[1] ^= xs[2]; xs[0] ^= xs[3]; xs[2] ^= t; xs[3] = rotl(xs[3], 45);
  return r;
}
static void rng_seed(unsigned int seed) {
  sm_x = (uint64_t) seed + 0x1234567ULL;
  for (int i = 0; i < 4; i++) xs[i] = splitmix64();
}
static double runif01(void) {              /* uniform (0,1) */
  return ((xoshiro() >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}
static double have_spare = 0.0; static int spare_ok = 0;
static double rnorm1(void) {               /* standard normal (Box-Muller) */
  if (spare_ok) { spare_ok = 0; return have_spare; }
  double u1 = runif01(), u2 = runif01();
  double r = sqrt(-2.0 * log(u1)), th = 6.283185307179586 * u2;
  have_spare = r * sin(th); spare_ok = 1;
  return r * cos(th);
}

/* ---- linear algebra (small p; column-major p-by-p matrices) ---- */
/* unbiased covariance of an n-by-p column-major data matrix into S (p-by-p). */
static void covariance(const double *X, int n, int p, double *S) {
  double *mu = (double*) calloc(p, sizeof(double));
  for (int j = 0; j < p; j++) { double s = 0; for (int i = 0; i < n; i++) s += X[i + (size_t)j*n]; mu[j] = s / n; }
  for (int a = 0; a < p; a++) for (int b = a; b < p; b++) {
    double s = 0;
    for (int i = 0; i < n; i++) s += (X[i + (size_t)a*n] - mu[a]) * (X[i + (size_t)b*n] - mu[b]);
    s /= (n - 1); S[a + (size_t)b*p] = s; S[b + (size_t)a*p] = s;
  }
  free(mu);
}
/* log-determinant via Cholesky of a p-by-p SPD matrix (column-major). Returns 0 on success,
 * 1 if not positive definite. Does not modify S. */
static int logdet_chol(const double *S, int p, double *out) {
  double *L = (double*) calloc((size_t)p*p, sizeof(double));
  int fail = 0;
  for (int j = 0; j < p && !fail; j++) {
    double d = S[j + (size_t)j*p];
    for (int k = 0; k < j; k++) d -= L[j + (size_t)k*p] * L[j + (size_t)k*p];
    if (d <= 0) { fail = 1; break; }
    double ljj = sqrt(d); L[j + (size_t)j*p] = ljj;
    for (int i = j+1; i < p; i++) {
      double s = S[i + (size_t)j*p];
      for (int k = 0; k < j; k++) s -= L[i + (size_t)k*p] * L[j + (size_t)k*p];
      L[i + (size_t)j*p] = s / ljj;
    }
  }
  double ld = 0; if (!fail) for (int j = 0; j < p; j++) ld += 2.0 * log(L[j + (size_t)j*p]);
  free(L); *out = ld; return fail;
}
/* Box's M from two within-group scatter/covariance inputs given as covariances SA, SB. */
static double mstat_from_cov(const double *SA, const double *SB, int nA, int nB, int p) {
  double *Sp = (double*) malloc((size_t)p*p*sizeof(double));
  double a = nA - 1, b = nB - 1, nu = a + b;
  for (int k = 0; k < p*p; k++) Sp[k] = (a*SA[k] + b*SB[k]) / nu;
  double ldp, lda, ldb; int f1 = logdet_chol(Sp, p, &ldp), f2 = logdet_chol(SA, p, &lda), f3 = logdet_chol(SB, p, &ldb);
  free(Sp);
  if (f1 || f2 || f3) return NAN;
  return nu*ldp - (a*lda + b*ldb);
}

void reg_mstat(const double *XA, const double *XB, const int *nA, const int *nB, const int *p, double *out) {
  int P = *p; double *SA = malloc((size_t)P*P*sizeof(double)), *SB = malloc((size_t)P*P*sizeof(double));
  covariance(XA, *nA, P, SA); covariance(XB, *nB, P, SB);
  *out = mstat_from_cov(SA, SB, *nA, *nB, P);
  free(SA); free(SB);
}

void reg_drawM(const int *p, const int *nA, const int *nB, const int *B, const int *seed, double *out) {
  int P = *p, na = *nA, nb = *nB, nrep = *B;
  rng_seed((unsigned int) *seed); spare_ok = 0;
  double *XA = malloc((size_t)na*P*sizeof(double)), *XB = malloc((size_t)nb*P*sizeof(double));
  double *SA = malloc((size_t)P*P*sizeof(double)), *SB = malloc((size_t)P*P*sizeof(double));
  for (int r = 0; r < nrep; r++) {
    for (size_t k = 0; k < (size_t)na*P; k++) XA[k] = rnorm1();
    for (size_t k = 0; k < (size_t)nb*P; k++) XB[k] = rnorm1();
    covariance(XA, na, P, SA); covariance(XB, nb, P, SB);
    out[r] = mstat_from_cov(SA, SB, na, nb, P);
  }
  free(XA); free(XB); free(SA); free(SB);
}

void reg_clx(const double *XA, const double *XB, const int *nA, const int *nB, const int *p, double *out) {
  int P = *p, na = *nA, nb = *nB;
  double *muA = calloc(P, sizeof(double)), *muB = calloc(P, sizeof(double));
  for (int j = 0; j < P; j++) { double sa=0, sb=0;
    for (int i = 0; i < na; i++) sa += XA[i + (size_t)j*na];
    for (int i = 0; i < nb; i++) sb += XB[i + (size_t)j*nb];
    muA[j] = sa/na; muB[j] = sb/nb; }
  double best = 0;
  for (int a = 0; a < P; a++) for (int b = a; b < P; b++) {
    double s1=0, s2=0, q1=0, q2=0;   /* MLE covariance (/n) and mean of squared products */
    for (int i = 0; i < na; i++) { double pr=(XA[i+(size_t)a*na]-muA[a])*(XA[i+(size_t)b*na]-muA[b]); s1+=pr; q1+=pr*pr; }
    for (int i = 0; i < nb; i++) { double pr=(XB[i+(size_t)a*nb]-muB[a])*(XB[i+(size_t)b*nb]-muB[b]); s2+=pr; q2+=pr*pr; }
    s1/=na; s2/=nb; double th1=q1/na - s1*s1, th2=q2/nb - s2*s2;
    double den = th1/na + th2/nb, num = (s1-s2)*(s1-s2);
    double v = den > 0 ? num/den : 0;
    if (v > best) best = v;
  }
  free(muA); free(muB); *out = best;
}

/* ---- exact null of Box's M by characteristic-function inversion ---- */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/* complex log-gamma: recurse to large real part, then Stirling series */
static double complex clgamma_c(double complex z) {
  double complex r = 0;
  while (creal(z) < 14.0) { r -= clog(z); z += 1.0; }
  double complex zi = 1.0 / z;
  return r + (z - 0.5) * clog(z) - z + 0.5 * log(2.0 * M_PI)
           + zi/12.0 - zi*zi*zi/360.0 + zi*zi*zi*zi*zi/1260.0;
}
/* characteristic function of M at real argument t (Box's exact null, product of gamma ratios) */
static double complex cf_M_c(double t, double nu1, double nu2, int p) {
  double nu = nu1 + nu2;
  double complex h = -2.0 * I * t;
  double logC = (nu*p/2.0)*log(nu) - (nu1*p/2.0)*log(nu1) - (nu2*p/2.0)*log(nu2);
  double complex s = h * logC;
  for (int j = 1; j <= p; j++) {
    double a = (j - 1) / 2.0;
    s += clgamma_c(nu/2.0 - a)          - clgamma_c(nu*(1.0+h)/2.0 - a)
       + clgamma_c(nu1*(1.0+h)/2.0 - a) - clgamma_c(nu1/2.0 - a)
       + clgamma_c(nu2*(1.0+h)/2.0 - a) - clgamma_c(nu2/2.0 - a);
  }
  return cexp(s);
}

void reg_pexact(const double *m, const double *nu1, const double *nu2, const int *p, double *out) {
  double M = *m, n1 = *nu1, n2 = *nu2; int P = *p;
  double EM = cimag(cf_M_c(1e-5, n1, n2, P)) / 1e-5;   /* E[M] from Phi'(0), for the t->0 endpoint */
  double T = 250.0; int N = 250000; double dt = T / N;
  double g0 = EM - M;                                  /* integrand limit as t -> 0 */
  double sum = 0.5 * g0, gl = 0.0;
  for (int k = 1; k <= N; k++) {
    double t = k * dt;
    double g = cimag(cexp(-I * t * M) * cf_M_c(t, n1, n2, P)) / t;
    if (k == N) gl = g; else sum += g;
  }
  sum += 0.5 * gl;                                      /* trapezoidal: half-weight the last node */
  double val = 0.5 + (1.0 / M_PI) * sum * dt;
  *out = val < 0.0 ? 0.0 : (val > 1.0 ? 1.0 : val);
}

/* --- Log-domain differential network (added for the functional-relationships paper) --- */

/* Symmetric eigendecomposition of a p-by-p matrix (column-major) by the cyclic Jacobi method.
 * On return, w holds the eigenvalues and V (p-by-p, column-major) their eigenvectors, so that
 * A = V diag(w) V'. A is copied internally and left unchanged. */
static void jacobi_sym(const double *A, int p, double *w, double *V) {
  double *a = (double*) malloc((size_t)p*p*sizeof(double));
  for (size_t i = 0; i < (size_t)p*p; i++) a[i] = A[i];
  for (int i = 0; i < p; i++) for (int j = 0; j < p; j++) V[i + (size_t)j*p] = (i == j) ? 1.0 : 0.0;
  for (int sweep = 0; sweep < 100; sweep++) {
    double off = 0.0;
    for (int q = 0; q < p; q++) for (int r = q+1; r < p; r++) off += a[q + (size_t)r*p]*a[q + (size_t)r*p];
    if (off < 1e-30) break;
    for (int q = 0; q < p; q++) for (int r = q+1; r < p; r++) {
      double apq = a[q + (size_t)r*p];
      if (fabs(apq) < 1e-300) continue;
      double app = a[q + (size_t)q*p], aqq = a[r + (size_t)r*p];
      double phi = 0.5 * atan2(2.0*apq, aqq - app);
      double c = cos(phi), s = sin(phi);
      for (int k = 0; k < p; k++) {
        double akq = a[k + (size_t)q*p], akr = a[k + (size_t)r*p];
        a[k + (size_t)q*p] = c*akq - s*akr;
        a[k + (size_t)r*p] = s*akq + c*akr;
      }
      for (int k = 0; k < p; k++) {
        double aqk = a[q + (size_t)k*p], ark = a[r + (size_t)k*p];
        a[q + (size_t)k*p] = c*aqk - s*ark;
        a[r + (size_t)k*p] = s*aqk + c*ark;
      }
      for (int k = 0; k < p; k++) {
        double vkq = V[k + (size_t)q*p], vkr = V[k + (size_t)r*p];
        V[k + (size_t)q*p] = c*vkq - s*vkr;
        V[k + (size_t)r*p] = s*vkq + c*vkr;
      }
    }
  }
  for (int i = 0; i < p; i++) w[i] = a[i + (size_t)i*p];
  free(a);
}

/* Matrix logarithm of a symmetric positive-definite p-by-p matrix (column-major) into out. */
static void logm_spd(const double *S, int p, double *out) {
  double *w = (double*) malloc((size_t)p*sizeof(double));
  double *V = (double*) malloc((size_t)p*p*sizeof(double));
  jacobi_sym(S, p, w, V);
  for (int i = 0; i < p; i++) for (int j = 0; j < p; j++) {
    double s = 0.0;
    for (int k = 0; k < p; k++) s += V[i + (size_t)k*p] * log(w[k]) * V[j + (size_t)k*p];
    out[i + (size_t)j*p] = s;
  }
  free(w); free(V);
}

/* Log-domain differential network D = log cov(XB) - log cov(XA) from two column-major data
 * matrices (XA is nA-by-p, XB is nB-by-p). Writes the p-by-p contrast D (column-major) to Dout.
 * D is inversion invariant: the same object results whether dependence is read through covariances
 * or precisions (log of a precision is minus log of the covariance). */
void reg_logdiff(const double *XA, const double *XB,
                 const int *nA, const int *nB, const int *p, double *Dout) {
  int P = *p;
  double *SA = (double*) malloc((size_t)P*P*sizeof(double));
  double *SB = (double*) malloc((size_t)P*P*sizeof(double));
  double *LA = (double*) malloc((size_t)P*P*sizeof(double));
  double *LB = (double*) malloc((size_t)P*P*sizeof(double));
  covariance(XA, *nA, P, SA); covariance(XB, *nB, P, SB);
  logm_spd(SA, P, LA); logm_spd(SB, P, LB);
  for (size_t i = 0; i < (size_t)P*P; i++) Dout[i] = LB[i] - LA[i];
  free(SA); free(SB); free(LA); free(LB);
}
