## regstat: an exact test for a change in covariance (dependence) structure.
## One shared pure-C back-end (compiled at install), also bound from Python (regstat on PyPI).
## The likelihood-ratio null is covariance-free (the real Jacobi ensemble), so a single simulation
## at the identity calibrates the test at every covariance, with no estimate of the nuisance covariance.
#' @useDynLib regstat, .registration = TRUE, .fixes = "C_"
#' @importFrom stats complete.cases
NULL

.as_mat <- function(X) {
  X <- as.matrix(X); storage.mode(X) <- "double"
  X[stats::complete.cases(X), , drop = FALSE]
}

#' Box's M statistic for equality of two covariance matrices (C back-end)
#'
#' @param XA numeric data matrix, group A (rows observations, columns variables).
#' @param XB numeric data matrix, group B, with the same number of columns.
#' @return the scalar likelihood-ratio statistic M; \code{NaN} if a group scatter is singular.
#' @export
cov_M <- function(XA, XB) {
  XA <- .as_mat(XA); XB <- .as_mat(XB)
  if (ncol(XA) != ncol(XB)) stop("XA and XB must have the same number of columns")
  .C(C_reg_mstat, XA = XA, XB = XB, nA = nrow(XA), nB = nrow(XB),
     p = ncol(XA), out = double(1))$out
}

#' Draws from the covariance-free null of M (C back-end)
#'
#' Simulates the exact Gaussian null of \code{\link{cov_M}} at the given dimensions. Because the null
#' is free of the common covariance, these draws (made at the identity) are the reference law at every
#' covariance.
#' @param p number of variables.
#' @param nA,nB per-group sample sizes.
#' @param B number of draws.
#' @param seed integer seed for the self-contained generator.
#' @return a numeric vector of \code{B} draws of M under equal covariance.
#' @export
cov_null <- function(p, nA, nB, B = 8000, seed = 1L) {
  .C(C_reg_drawM, p = as.integer(p), nA = as.integer(nA), nB = as.integer(nB),
     B = as.integer(B), seed = as.integer(seed), out = double(B))$out
}

#' Exact right-tail probability of Box's M (C back-end)
#'
#' Evaluates \eqn{P(M > m)} under the covariance-free null deterministically, by numerical inversion of
#' the exact characteristic function -- no simulation. Because the null is covariance-free the result
#' depends only on the dimensions.
#' @param m observed statistic value.
#' @param nu1,nu2 within-group degrees of freedom (sample size minus one).
#' @param p number of variables.
#' @return the exact tail probability.
#' @export
cov_pexact <- function(m, nu1, nu2, p)
  .C(C_reg_pexact, m = as.double(m), nu1 = as.double(nu1), nu2 = as.double(nu2),
     p = as.integer(p), out = double(1))$out

#' Exact test for a change in covariance structure
#'
#' Tests H0 that two groups share a covariance matrix, using the likelihood-ratio statistic M and its
#' covariance-free null. By default the null is evaluated exactly and deterministically by inverting the
#' characteristic function (\code{method = "exact"}); \code{method = "mc"} uses a Monte Carlo null
#' instead. Either way the test needs no estimate of the common covariance, its advantage over
#' high-dimensional competitors near the dimension barrier where that covariance is hardest to estimate.
#' @param XA,XB numeric data matrices with the same number of columns.
#' @param method \code{"exact"} (CF inversion) or \code{"mc"} (Monte Carlo).
#' @param B number of null draws when \code{method = "mc"}.
#' @param seed integer seed when \code{method = "mc"}.
#' @return an object of class \code{"htest"}.
#' @export
cov_test <- function(XA, XB, method = c("exact", "mc"), B = 8000, seed = 1L) {
  method <- match.arg(method)
  XA <- .as_mat(XA); XB <- .as_mat(XB)
  Mo <- cov_M(XA, XB)
  if (method == "exact") {
    p <- cov_pexact(Mo, nrow(XA) - 1, nrow(XB) - 1, ncol(XA))
    meth <- "Exact covariance-change test (covariance-free Jacobi null, CF inversion)"
  } else {
    nul <- cov_null(ncol(XA), nrow(XA), nrow(XB), B = B, seed = seed)
    p <- (1 + sum(nul >= Mo)) / (B + 1)
    meth <- "Exact covariance-change test (covariance-free Jacobi null, Monte Carlo)"
  }
  structure(list(statistic = c(M = Mo), p.value = p, method = meth,
                 data.name = paste(deparse(substitute(XA)), "and", deparse(substitute(XB)))),
            class = "htest")
}

#' Cai--Liu--Xia max-type covariance statistic (C back-end)
#'
#' @param XA,XB numeric data matrices with the same number of columns.
#' @return the maximum standardised squared difference of covariance entries.
#' @export
clx_stat <- function(XA, XB) {
  XA <- .as_mat(XA); XB <- .as_mat(XB)
  if (ncol(XA) != ncol(XB)) stop("XA and XB must have the same number of columns")
  .C(C_reg_clx, XA = XA, XB = XB, nA = nrow(XA), nB = nrow(XB),
     p = ncol(XA), out = double(1))$out
}

#' Max-type two-sample covariance test with the extreme-value calibration
#'
#' The Cai--Liu--Xia (2013) test as used out of the box. Its extreme-value calibration is accurate at
#' large dimension but inflates near the dimension barrier; see \code{\link{cov_test}} for the exact,
#' covariance-free alternative.
#' @param XA,XB numeric data matrices with the same number of columns.
#' @param alpha nominal level.
#' @return an object of class \code{"htest"}.
#' @export
clx_test <- function(XA, XB, alpha = 0.05) {
  XA <- .as_mat(XA); XB <- .as_mat(XB); p <- ncol(XA)
  M <- clx_stat(XA, XB)
  q <- 4 * log(p) - log(log(p)) - 2 * log(-sqrt(8 * pi) * log(1 - alpha))
  pval <- 1 - exp(-(1 / sqrt(8 * pi)) * exp(-(M - 4 * log(p) + log(log(p))) / 2))
  structure(list(statistic = c(Mmax = M), p.value = pval,
                 method = "Max-type two-sample covariance test (Cai-Liu-Xia, extreme-value null)",
                 data.name = paste(deparse(substitute(XA)), "and", deparse(substitute(XB)))),
            class = "htest")
}
