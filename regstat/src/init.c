/* Native-routine registration for regstat (CRAN-compliant .C interface). */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "regcov.h"

static const R_CMethodDef CEntries[] = {
  {"reg_mstat",  (DL_FUNC) &reg_mstat,  6},
  {"reg_drawM",  (DL_FUNC) &reg_drawM,  6},
  {"reg_clx",    (DL_FUNC) &reg_clx,    6},
  {"reg_pdet",     (DL_FUNC) &reg_pdet,     5},
  {"reg_pexact", (DL_FUNC) &reg_pexact, 5},
  {"reg_logdiff", (DL_FUNC) &reg_logdiff, 6},
  {NULL, NULL, 0}
};

void R_init_regstat(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
