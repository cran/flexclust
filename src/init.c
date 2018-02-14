#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void countPairs(int *c1, int *c2, int *n, double *res);
extern void hardcl(int *xrows, int *xcols, double *x, int *ncenters,
	       double *centers, int *cluster,
	       int *itermax, int *iter, 
	       int *clustersize, int *verbose, int *dist, int *methrate,
	       double *par, double *weights);
extern void kmeans(int *xrows, int *xcols, double *x, int *ncenters,
	           double *centers, int *cluster,
	           int *itermax, int *iter, int *changes,
	           int *clustersize, int *verbose, int *dist);
extern void neuralgas(int *xrows, int *xcols, double *x, int *ncenters,
	              double *centers, int *cluster,int *itermax, int *iter,
	              int *clustersize, int *verbose, int *dist,double *par);

static const R_CMethodDef CEntries[] = {
    {"countPairs", (DL_FUNC) &countPairs,  4},
    {"hardcl",     (DL_FUNC) &hardcl,     14},
    {"kmeans",     (DL_FUNC) &kmeans,     12},
    {"neuralgas",  (DL_FUNC) &neuralgas,  12},
    {NULL, NULL, 0}
};

void R_init_flexclust(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
