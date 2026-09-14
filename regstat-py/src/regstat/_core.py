"""regstat: Python (ctypes) binding to the shared covariance-change C back-end. Pure standard library.

The same C source backs the R package 'regstat'; a parity harness in the repository requires every
exported quantity to agree between the two front ends. The exact test's null is covariance-free, so a
single simulation at the identity calibrates it at every covariance.

LOADING, and why it has two paths. An INSTALLED package carries a library that pip built at install
time from setup.py -- no compiler is needed here and nothing is written to disk. A SOURCE CHECKOUT
has the .c file but no built library, so it is compiled on demand into a cache directory to keep the
edit-a-C-file-and-re-import development loop working. The previous version compiled on every fresh
import AND wrote into its own site-packages directory, which fails on read-only installs and needs a
compiler at import rather than install.
"""
import ctypes as _ct
import glob as _glob
import math as _math
import os as _os
import sys as _sys

_here = _os.path.dirname(_os.path.abspath(__file__))

def _find_built():
    """A library placed beside this file by pip at install time."""
    for pat in ("_libregstat*.so", "_libregstat*.dylib", "_libregstat*.pyd", "_libregstat*.dll"):
        hits = sorted(_glob.glob(_os.path.join(_here, pat)))
        if hits:
            return hits[0]
    return None

def _cache_dir():
    """Somewhere writable, preferring the package directory only if it actually is."""
    if _os.access(_here, _os.W_OK):
        return _here
    base = _os.environ.get("XDG_CACHE_HOME") or _os.path.join(_os.path.expanduser("~"), ".cache")
    d = _os.path.join(base, "regstat")
    _os.makedirs(d, exist_ok=True)
    return d

def _lib_ext():
    return "dll" if _sys.platform == "win32" else ("dylib" if _sys.platform == "darwin" else "so")

def _compile_from_source():
    """Development fallback: build the bundled source. Never reached in an installed package."""
    import subprocess as _sp
    srcs = sorted(f for f in _glob.glob(_os.path.join(_here, "*.c"))
                  if _os.path.basename(f) != "init.c")
    if not srcs:
        return None
    out = _os.path.join(_cache_dir(), "libregstat." + _lib_ext())
    cc = _os.environ.get("CC") or "cc"
    try:
        _sp.check_call([cc, "-O2", "-fPIC", "-shared", "-I", _here] + srcs + ["-lm", "-o", out])
    except (OSError, _sp.CalledProcessError) as exc:
        raise ImportError(
            "regstat could not build its C back-end from the bundled source (%s). In an installed "
            "copy this should never happen -- pip builds the library at install time. In a source "
            "checkout, a C compiler is required; set CC if it is not 'cc'." % exc
        )
    return out

_libpath = _find_built()
if _libpath is None:
    _cached = _os.path.join(_cache_dir(), "libregstat." + _lib_ext())
    _libpath = _cached if _os.path.exists(_cached) else _compile_from_source()
if _libpath is None:
    raise ImportError("regstat: no compiled back-end and no C source to build one from.")

_lib = _ct.CDLL(_libpath)
_dp, _ip = _ct.POINTER(_ct.c_double), _ct.POINTER(_ct.c_int)
_lib.reg_mstat.argtypes  = [_dp, _dp, _ip, _ip, _ip, _dp]
_lib.reg_drawM.argtypes  = [_ip, _ip, _ip, _ip, _ip, _dp]
_lib.reg_clx.argtypes    = [_dp, _dp, _ip, _ip, _ip, _dp]
_lib.reg_pexact.argtypes = [_dp, _dp, _dp, _ip, _dp]
_lib.reg_pdet.argtypes   = [_dp, _dp, _dp, _ip, _dp]
_lib.reg_logdiff.argtypes = [_dp, _dp, _ip, _ip, _ip, _dp]

def _i(v):
    a = (_ct.c_int * 1)(); a[0] = int(v); return a
def _d1(v):
    a = (_ct.c_double * 1)(); a[0] = float(v); return a
def _mat(X):
    """Column-major ctypes double array from a 2-D sequence (numpy array or list of rows)."""
    rows = [list(map(float, r)) for r in X]
    n = len(rows); p = len(rows[0]) if n else 0
    arr = (_ct.c_double * (n * p))()
    for j in range(p):
        base = j * n
        for i in range(n):
            arr[base + i] = rows[i][j]
    return arr, n, p

def cov_M(XA, XB):
    """Box's M statistic for equality of the two groups' covariance matrices."""
    a, nA, p = _mat(XA); b, nB, pb = _mat(XB)
    if p != pb: raise ValueError("XA and XB must have the same number of columns")
    out = _d1(0.0); _lib.reg_mstat(a, b, _i(nA), _i(nB), _i(p), out); return out[0]

def cov_null(p, nA, nB, B=8000, seed=1):
    """Draws from the covariance-free null of M at the given dimensions."""
    out = (_ct.c_double * int(B))()
    _lib.reg_drawM(_i(p), _i(nA), _i(nB), _i(B), _i(seed), out); return list(out)

def cov_pexact(m, nu1, nu2, p):
    """Exact right-tail probability P(M > m) by CF inversion (deterministic, no simulation)."""
    out = _d1(0.0); _lib.reg_pexact(_d1(m), _d1(nu1), _d1(nu2), _i(p), out); return out[0]

def cov_pdet(x, nu1, nu2, p):
    """Exact upper-tail probability of D = log(det(A2)/det(A1)) under the covariance-free null,
    by the same safeguarded Gil-Pelaez inversion as cov_pexact. Same bytes as the R front."""
    out = _d1(0.0); _lib.reg_pdet(_d1(x), _d1(nu1), _d1(nu2), _i(p), out); return out[0]

def cov_powdet(lam, crit, nu1, nu2, p):
    """Exact power of the determinant chart against the proportional alternative
    Sigma2 = lam * Sigma1: the statistic shifts by exactly p*log(lam)."""
    return [cov_pdet(crit - p * _math.log(l), nu1, nu2, p) for l in lam]

def cov_test(XA, XB, method="exact", B=8000, seed=1):
    """Exact covariance-change test: statistic M against its covariance-free null. By default the null
    is evaluated deterministically by characteristic-function inversion (method='exact'); method='mc'
    uses a Monte Carlo null. Returns a dict with the statistic, p-value, and method."""
    a, nA, p = _mat(XA); b, nB, pb = _mat(XB)
    if p != pb: raise ValueError("XA and XB must have the same number of columns")
    Mo = cov_M(XA, XB)
    if method == "exact":
        pv = cov_pexact(Mo, nA - 1, nB - 1, p)
        meth = "Exact covariance-change test (covariance-free Jacobi null, CF inversion)"
    else:
        nul = cov_null(p, nA, nB, B, seed)
        pv = (1 + sum(1 for v in nul if v >= Mo)) / (B + 1)
        meth = "Exact covariance-change test (covariance-free Jacobi null, Monte Carlo)"
    return {"statistic": Mo, "p_value": pv, "method": meth}

def cov_logdiff(XA, XB):
    """Log-domain differential network D = log cov(XB) - log cov(XA), returned as a column-major
    list of p lists (columns). Inversion invariant: identical whether dependence is read through
    covariances or precisions, so it dissolves the covariance-versus-precision choice."""
    a, nA, p = _mat(XA); b, nB, pb = _mat(XB)
    if p != pb: raise ValueError("XA and XB must have the same number of columns")
    out = (_ct.c_double * (p * p))()
    _lib.reg_logdiff(a, b, _i(nA), _i(nB), _i(p), out)
    return [[out[i + j * p] for i in range(p)] for j in range(p)]

def clx_stat(XA, XB):
    """Cai-Liu-Xia max-type two-sample covariance statistic."""
    a, nA, p = _mat(XA); b, nB, pb = _mat(XB)
    if p != pb: raise ValueError("XA and XB must have the same number of columns")
    out = _d1(0.0); _lib.reg_clx(a, b, _i(nA), _i(nB), _i(p), out); return out[0]

def clx_test(XA, XB, alpha=0.05):
    """Max-type (Cai-Liu-Xia) covariance test with the extreme-value calibration."""
    _, _, p = _mat(XA); M = clx_stat(XA, XB)
    pv = 1.0 - _math.exp(-(1.0 / _math.sqrt(8 * _math.pi)) *
                         _math.exp(-(M - 4 * _math.log(p) + _math.log(_math.log(p))) / 2))
    return {"statistic": M, "p_value": pv,
            "method": "Max-type two-sample covariance test (Cai-Liu-Xia, extreme-value null)"}
