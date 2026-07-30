# regstat

An exact test for a change in covariance (dependence) structure — the omnibus form of the
differential-network question. Under the Gaussian null the likelihood-ratio statistic has a
distribution given by the real Jacobi ensemble that is **free of the unknown common covariance**, so a
single Monte-Carlo calibration at the identity serves every covariance with no estimate of the nuisance
covariance. This is the property that survives the dimension barrier, where estimating the covariance is
hardest. The Cai–Liu–Xia (2013) max-type high-dimensional test is included for comparison.

One shared pure-C back-end does the numerics and also backs the Python package `regstat`.

```python
import numpy as np; from regstat import cov_test, clx_test
XA = np.random.randn(30, 8); XB = np.random.randn(30, 8)
cov_test(XA, XB)   # dict with statistic and p_value
clx_test(XA, XB)
```
