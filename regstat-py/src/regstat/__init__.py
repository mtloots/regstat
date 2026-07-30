"""regstat: an exact test for a change in covariance (dependence) structure.

The likelihood-ratio null is covariance-free (the real Jacobi ensemble), so a single Monte-Carlo
calibration at the identity serves every covariance with no estimate of the nuisance covariance --
the property that survives the dimension barrier. Built on the same pure-C back-end as the R package.
"""
from ._core import cov_M, cov_null, cov_pexact, cov_test, clx_stat, clx_test

__all__ = ["cov_M", "cov_null", "cov_pexact", "cov_test", "clx_stat", "clx_test"]
__version__ = "0.1.0"
