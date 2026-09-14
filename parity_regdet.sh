#!/bin/sh
# Parity harness for the exact tail evaluations (regcov.c): cov_pexact and cov_pdet on a grid
# of arguments, byte-identical between the R and Python fronts. Both fronts compile the same
# C source; the guard below refuses to run against stale copies.
set -e
OUT=${TMPDIR:-/tmp}/parity_regdet
mkdir -p "$OUT"
PKG="$(cd "$(dirname "$0")" && pwd)"
sh "$PKG/parity_sync_guard.sh" "$PKG/regstat/src" "$PKG/regstat-py/src/regstat" || exit 1
Rscript -e '
suppressMessages(library(regstat))
vals <- c(
  vapply(c(5, 15, 35), function(m) cov_pexact(m, 59, 49, 3), 0),
  vapply(c(-2, 0, 1, 3), function(x) cov_pdet(x, 15, 15, 12), 0),
  vapply(c(-1, 0.5), function(x) cov_pdet(x, 59, 49, 1), 0),
  cov_powdet(c(0.5, 1, 1.5, 2), 1.0, 15, 15, 12))
cat(sprintf("%.17g", vals), sep="\n")
' > "$OUT/r.txt"
PYTHONPATH="$PKG/regstat-py/src" python3 - <<PY > "$OUT/py.txt"
import math
from regstat import cov_pexact, cov_pdet, cov_powdet
vals = ([cov_pexact(m, 59, 49, 3) for m in (5, 15, 35)]
        + [cov_pdet(x, 15, 15, 12) for x in (-2, 0, 1, 3)]
        + [cov_pdet(x, 59, 49, 1) for x in (-1, 0.5)]
        + cov_powdet([0.5, 1, 1.5, 2], 1.0, 15, 15, 12))
print("\n".join(("NaN" if math.isnan(v) else "%.17g" % v) for v in vals))
PY
if cmp -s "$OUT/r.txt" "$OUT/py.txt"; then
  echo "VERDICT: PARITY -- $(wc -l < "$OUT/r.txt" | tr -d ' ') values byte-identical"
else
  echo "VERDICT: MISMATCH"; diff "$OUT/r.txt" "$OUT/py.txt" | head; exit 1
fi
