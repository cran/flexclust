#ifndef _FLEXCLUST_H
#define _FLEXCLUST_H

void kmeans(int *xrows, int *xcols, double *x, int *ncenters,
 	    double *centers, int *cluster,
	    int *itermax, int *iter, int *changes,
	    int *clustersize, int *verbose, int *dist);
void hardcl(int *xrows, int *xcols, double *x, int *ncenters,
	    double *centers, int *cluster,
	    int *itermax, int *iter, 
	    int *clustersize, int *verbose, int *dist, int *methrate,
	    double *par, double *weights);
void neuralgas(int *xrows, int *xcols, double *x, int *ncenters,
	       double *centers, int *cluster, int *itermax, int *iter,
	       int *clustersize, int *verbose, int *dist, double *par);
void predstr_calcSums(int pred[], int tc[], int *length, int *levels, int res[]);
void countPairs(int *c1, int *c2, int *n, double *res);

#endif
