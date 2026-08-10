/* regcov.h : shared pure-C back-end for the exact covariance-change test.
 *
 * One C ABI, bound identically by R (.C) and Python (ctypes). Every entry point takes
 * only pointer arguments so the SAME compiled library serves both bindings. Data matrices
 * are column-major (an n-by-p matrix stored by columns, as R stores matrices and as numpy
 * stores order='F' arrays). The random-number generator is self-contained and seeded by an
 * argument, so no R runtime is required: the null of the statistic is covariance-free, so a
 * single simulation at the identity calibrates the test at every covariance.
 */
#ifndef REGCOV_H
#define REGCOV_H
#ifdef __cplusplus
extern "C" {
#endif

/* Box's M statistic for equality of two covariance matrices, from two data matrices.
 * XA is nA-by-p column-major, XB is nB-by-p column-major. Writes M to *out;
 * writes NaN if a within-group scatter matrix is not positive definite. */
void reg_mstat(const double *XA, const double *XB,
               const int *nA, const int *nB, const int *p, double *out);

/* B Monte-Carlo draws of M under the covariance-free Gaussian null (equal covariance),
 * generated with the self-contained RNG seeded by *seed. Writes out[0..B-1]. Because the
 * null is free of the common covariance, these draws (made at the identity) calibrate the
 * test at every covariance. */
void reg_drawM(const int *p, const int *nA, const int *nB,
               const int *B, const int *seed, double *out);

/* Cai--Liu--Xia (2013) max-type two-sample covariance statistic, from two data matrices
 * (column-major, as above). Writes the maximum standardised squared entry difference to *out. */
void reg_clx(const double *XA, const double *XB,
             const int *nA, const int *nB, const int *p, double *out);

/* Exact right-tail probability P(M > m) of Box's M under the covariance-free Gaussian null,
 * evaluated DETERMINISTICALLY by numerical inversion of the exact characteristic function (a product
 * of gamma ratios). nu1 = nA - 1, nu2 = nB - 1 are the within-group degrees of freedom. Writes the
 * tail probability to *out. This is the exact, simulation-free evaluation of the null. */
void reg_pexact(const double *m, const double *nu1, const double *nu2, const int *p, double *out);

/* Log-domain differential network D = log cov(XB) - log cov(XA) from two column-major data
 * matrices (XA is nA-by-p, XB is nB-by-p). Writes the p-by-p contrast D (column-major) to Dout.
 * Inversion invariant: identical whether dependence is read through covariances or precisions. */
void reg_logdiff(const double *XA, const double *XB,
                 const int *nA, const int *nB, const int *p, double *Dout);

#ifdef __cplusplus
}
#endif
#endif
