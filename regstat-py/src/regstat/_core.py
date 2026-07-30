"""regstat: Python (ctypes) binding to the shared covariance-change C back-end. Pure standard library.
The C source is compiled on first import; the same source backs the R package 'regstat'. The exact
test's null is covariance-free, so a single simulation at the identity calibrates it at every covariance."""
import ctypes as _ct, os as _os, math as _math
_here = _os.path.dirname(_os.path.abspath(__file__))
_ext = "dylib" if _os.uname().sysname == "Darwin" else "so"
_libpath = _os.path.join(_here, "libregstat." + _ext)
if not _os.path.exists(_libpath):
    import subprocess as _sp
    _sp.check_call(["cc", "-O2", "-fPIC", "-shared", "-I", _here,
                    _os.path.join(_here, "regcov.c"), "-lm", "-o", _libpath])
_lib = _ct.CDLL(_libpath)
_dp, _ip = _ct.POINTER(_ct.c_double), _ct.POINTER(_ct.c_int)
_lib.reg_mstat.argtypes  = [_dp, _dp, _ip, _ip, _ip, _dp]
_lib.reg_drawM.argtypes  = [_ip, _ip, _ip, _ip, _ip, _dp]
_lib.reg_clx.argtypes    = [_dp, _dp, _ip, _ip, _ip, _dp]
_lib.reg_pexact.argtypes = [_dp, _dp, _dp, _ip, _dp]

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
