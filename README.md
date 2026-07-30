# regstat

An exact, covariance-free test for a change in covariance (dependence) structure — the omnibus form of
the differential-network question. The likelihood-ratio null (the real Jacobi ensemble) is free of the
unknown common covariance, so one calibration serves every covariance with no estimate of the nuisance
covariance; the null is evaluated exactly and deterministically by characteristic-function inversion,
verified against large-scale simulation to the 0.1 per cent tail.

One shared pure-C back-end (`regcov.c`/`.h`) serves two front-ends:

	regstat/      R package    (.C bindings; R CMD check: Status OK)
	regstat-py/   Python package (ctypes; compiles the same C on first import)

Companion software to the manuscript "An exact, covariance-free test for a change in dependence
structure, evaluated by characteristic-function inversion" (M. T. Loots, submitted).
