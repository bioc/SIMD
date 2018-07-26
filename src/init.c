#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/
  
  /* .C calls */
extern void pvalueclassify(void *, void *, void *, void *, void *, void *, 
                           void *, void *, void *, void *, void *, void *);

static const R_CMethodDef cMethods[] = {
  {"pvalueclassify", (DL_FUNC) & pvalueclassify, 12},
  {NULL}
};

void R_init_SIMD(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
