#include <R.h>
#include <R_ext/Rdynload.h>

#include "flexclust.h"

#define CDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CMethodDef cMethods[] = {
    CDEF(kmeans, 12),
    CDEF(hardcl, 14),
    CDEF(neuralgas, 12),
    CDEF(predstr_calcSums, 5),
    CDEF(countPairs, 4),
    {NULL, NULL, 0}
};

void R_init_flexclust(DllInfo *dll)
{
    R_useDynamicSymbols(dll, FALSE);
    R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
}
