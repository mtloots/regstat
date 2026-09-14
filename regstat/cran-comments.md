## Submission

This is a first submission of 'regstat'.

The package provides an exact finite-sample test for a change in covariance structure. Under the
Gaussian null the likelihood-ratio statistic follows a real Jacobi-ensemble law that does not depend
on the unknown common covariance, so one Monte-Carlo calibration at the identity serves every
covariance.

Two points raised by the CRAN team on my first submission of 'arcstat' in August 2026 were taken as
standing requirements here rather than waited for.

**`\value` tags.** Every `.Rd` documenting an exported function states what the function returns and
what it means. The package-level page `regstat-package.Rd` has no `\value`, which is correct for an
overview page. A check of the built tarball reports no undocumented code objects and no
code/documentation mismatches.

**References in the Description.** The comparison method is cited in the form CRAN asks for:
Cai, Liu and Xia (2013) <doi:10.1080/01621459.2012.758041>. The DOI was confirmed against Crossref
rather than copied from a bibliography.

The exact test itself is described in a manuscript under review at Metrika. It has no DOI yet, so it
is not cited in the Description; I would rather cite nothing than cite something a reader cannot
reach, and I will add it as soon as there is a DOI to give. This is the same position I took for
'arcstat', and the CRAN cookbook notes that a reference in the Description is optional.

Software and package names are single-quoted in the Description ('Python', 'regstat'), which the
cookbook asks for.

## Test environments

* local macOS 26 (arm64), R 4.6.1, `R CMD check --as-cran` on the built tarball

Anything beyond that line will be added only once it has actually been run. The 'arcstat'
submission file records why: an earlier version of it claimed a GitHub Actions matrix covering
R-devel, oldrel, macOS and Windows, when the matrix job carries
`if: github.event_name == 'workflow_dispatch'` and does not run on an ordinary push. Before this is
submitted I will either dispatch that workflow and wait for it, or run win-builder on this exact
tarball, and then say which.

## R CMD check results

0 errors | 0 warnings | 1 note

The note is "New submission", which is expected for a first submission.

## Notes for the reviewer

* The package contains compiled C. It uses only libm and R's own headers, needs no system
  requirement beyond a C compiler, and registers its native routines with `R_registerRoutines`
  and `R_useDynamicSymbols(dll, FALSE)`.
* The same C sources back a 'Python' package of the same name. The two front ends are compared
  value by value by a parity harness in the repository, and the tests here encode the mathematical
  identities the routines must satisfy rather than recording their current output.
* No example, test or vignette uses more than two cores.
