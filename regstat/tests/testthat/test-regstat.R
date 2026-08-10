## The exact covariance-change test. The properties worth pinning down are the ones the paper's
## argument rests on: Box's M vanishes when nothing differs and is invariant to a common change
## of basis, the exact null agrees with simulation, and the log-difference network is
## antisymmetric and symmetric-matrix-valued.

make_pair <- function(p = 3, nA = 60, nB = 50, scale = c(1, 1.6, 0.7), seed = 1) {
  set.seed(seed)
  list(XA = matrix(rnorm(nA * p), nA, p),
       XB = matrix(rnorm(nB * p), nB, p) %*% diag(scale[seq_len(p)]))
}

test_that("Box's M vanishes when the two samples are the same", {
  d <- make_pair()
  expect_equal(cov_M(d$XA, d$XA), 0)
})

test_that("Box's M is invariant to a common change of basis", {
  ## the statistic compares covariance STRUCTURE, so applying one nonsingular map to both
  ## samples must leave it alone. This is what lets the null be calibrated at the identity.
  d <- make_pair()
  m0 <- cov_M(d$XA, d$XB)
  set.seed(99)
  for (i in 1:3) {
    A <- matrix(rnorm(9), 3, 3)
    if (abs(det(A)) < 1e-6) next
    expect_lt(abs(cov_M(d$XA %*% A, d$XB %*% A) - m0), 1e-8 * max(1, abs(m0)))
  }
})

test_that("Box's M grows as the covariances are pulled apart", {
  m <- vapply(c(1, 1.2, 1.8, 3), function(s) {
    d <- make_pair(scale = c(1, s, 1)); cov_M(d$XA, d$XB)
  }, 0)
  expect_true(all(diff(m) > 0))
})

test_that("the exact null is a distribution function", {
  ## cov_pexact evaluates the null tail by characteristic-function inversion rather than by
  ## simulation, which is the paper's contribution; it must behave like a tail probability.
  pv <- vapply(c(1, 5, 15, 35, 60), function(m) cov_pexact(m, 59, 49, 3), 0)
  expect_true(all(pv >= 0 & pv <= 1))
  expect_false(is.unsorted(rev(pv)))                    # decreasing in the statistic
})

test_that("the exact null agrees with the simulated one", {
  ## the two routes are independent: inversion of the exact characteristic function against
  ## covariance-free Monte Carlo draws. They must agree in the body of the distribution.
  skip_on_cran()
  set.seed(3)
  p <- 3; nA <- 40; nB <- 40
  null <- cov_null(p, nA, nB, B = 20000, seed = 5L)
  for (q in c(0.5, 0.75, 0.9)) {
    m <- unname(quantile(null, q))
    expect_lt(abs(cov_pexact(m, nA - 1, nB - 1, p) - (1 - q)), 0.02)
  }
})

test_that("the test returns an htest with a probability in it", {
  d <- make_pair()
  for (meth in c("exact", "mc")) {
    ct <- cov_test(d$XA, d$XB, method = meth, B = 4000)
    expect_true(all(c("statistic", "p.value", "method") %in% names(ct)))
    expect_gte(ct$p.value, 0)
    expect_lte(ct$p.value, 1)
  }
  ## a genuine scale change in one coordinate must be detected
  expect_lt(cov_test(d$XA, d$XB, method = "exact")$p.value, 0.01)
  ## and two samples from the same law must not be
  set.seed(11)
  YA <- matrix(rnorm(60 * 3), 60, 3); YB <- matrix(rnorm(60 * 3), 60, 3)
  expect_gt(cov_test(YA, YB, method = "exact")$p.value, 0.01)
})

test_that("the log-difference network is antisymmetric and symmetric-matrix-valued", {
  ## D = log cov(XB) - log cov(XA). Swapping the samples must negate it exactly, and each D
  ## is a difference of symmetric matrix logarithms so it is symmetric.
  d <- make_pair()
  D1 <- cov_logdiff(d$XA, d$XB)
  D2 <- cov_logdiff(d$XB, d$XA)
  expect_lt(max(abs(D1 + D2)), 1e-12)
  expect_lt(max(abs(D1 - t(D1))), 1e-12)
  expect_lt(max(abs(cov_logdiff(d$XA, d$XA))), 1e-12)
})

test_that("the max-type statistic and its test behave", {
  d <- make_pair()
  expect_true(is.finite(clx_stat(d$XA, d$XB)))
  expect_gte(clx_stat(d$XA, d$XB), 0)
  ct <- clx_test(d$XA, d$XB)
  expect_true(is.list(ct) || is.numeric(ct))
})

test_that("the differential-network test returns a probability", {
  skip_on_cran()
  d <- make_pair()
  dt <- diffnet_test(d$XA, d$XB, B = 200, seed = 2L)
  pv <- if (is.list(dt)) dt$p.value else dt
  expect_gte(pv, 0)
  expect_lte(pv, 1)
})
