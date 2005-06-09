/*****************************************************************/
/*
 *  Copyright (C) 1998-2004 Evgenia Dimitriadou, Friedrich Leisch
 *                          and Andreas Weingessel
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <R_ext/Utils.h>

/**********************************************************/

double median(double *x, int n)
{
  double xmed;
  int n2;

  R_rsort (x, n);
  n2 = n / 2;
  if (n2 << 1 == n) {
	  xmed = (x[n2] + x[n2 + 1]) * .5;
  } else {
	  xmed = x[n2 + 1];
  }
  return xmed;
}


/**********************************************************/


int assign(int *xrows, int *xcols, double *x, int *ncenters,
	   double *centers, int *cluster, int *clustersize,
	   int *dist)
{
  int k, m, n;
  double dista, mindist;

  for(k=0; k<*xrows; k++){
    mindist=10e99;

    for(m=0; m< *ncenters; m++){
      dista=0;
      for(n=0;n<*xcols;n++){
	if(*dist == 0){
	  dista += (x[k+(*xrows)*n] - centers[m + (*ncenters)*n])*(x[k+(*xrows)*n] - centers[m + (*ncenters)*n]);
	}
	else if(*dist ==1){
	  dista += fabs(x[k+(*xrows)*n] - centers[m + (*ncenters)*n]);
	}
	
      }

      if(dista<mindist){
	cluster[k] = m;
	mindist = dista;
      }
    }
  }

  
  for(k=0;k<*ncenters;k++){
    clustersize[k]=0;
  } 
  for(k=0; k<*xrows; k++){
    clustersize[cluster[k]]++;
  }

  return 0;
}

/**********************************************************/
		 
int reloc(int *xrows, int *xcols, double *x, int *ncenters,
	  double *centers, int *cluster, int *clustersize,
	  int *dist)
{
  int k,l,m; 
  

  /* Initialize cluster size and centers with 0 */
  for(k=0;k<*ncenters;k++){
/*      clustersize[k]=0;*/
      for(m=0;m<*xcols;m++){
	  centers[k+(*ncenters)*m] = 0;
      }
  }
  
/*      for(k=0; k<*xrows; k++)
	clustersize[cluster[k]]++; */
  
  if(*dist==0){
    /* Euclidean distance */
    for(k=0; k<*xrows; k++){
      for(m=0;m<*xcols;m++){
	centers[cluster[k]+(*ncenters)*m] += x[k+(*xrows)*m];
      }
    }
    
    for(k=0; k<*ncenters; k++){
      for(m=0;m<*xcols;m++){
	centers[k+(*ncenters)*m] /= clustersize[k];
      }
    }
  }
  else if(*dist == 1){
    for(k=0; k<*ncenters; k++){
	double *xk;
	int i;
	xk = (double *) R_alloc(clustersize[k], sizeof(double));
	
	for (m=0; m<*xcols; m++){
	    i=0;
	    for (l=0; l<*xrows; l++)
		if(cluster[l]==k){
		    xk[i] = x[l+(*xrows)*m];
		    i++;
		}
	    centers[k+(*ncenters)*m] = median(xk, clustersize[k]);
	}
    }
  }
  return 0;
}

/**********************************************************/

int kmeans(int *xrows, int *xcols, double *x, int *ncenters,
	   double *centers, int *cluster,
	   int *itermax, int *iter, int *changes,
	   int *clustersize, int *verbose, int *dist)
{
  int m;
  int change;
  int *clustnew;

  clustnew = (int *) R_alloc(*xrows, sizeof(int));
  change = 1;
  *iter=0;
  while(change && ((*iter)++ < *itermax)){
    assign(xrows, xcols, x, ncenters, centers, clustnew, clustersize, dist);
    reloc(xrows, xcols, x, ncenters, centers, clustnew, clustersize, dist);
    change = 0;
    for(m=0; m<*xrows; m++){
      if(cluster[m] != clustnew[m]){
	change++;
	cluster[m] = clustnew[m];
      }
    }
    if(*verbose){
      Rprintf("Iteration: %3d    Changes: %13d \n", *iter, change);
    }
    changes[(*iter)-1] = change;
  }
  return 0;
}
     
/**********************************************************/

int  oncent(int *xrows, int *xcols, double *x, int *ncenters,
	   double *centers, int *cluster, int *clustersize,
	   int *dist,int *iter,int *itermax,int *methrate,
	   double *par,int *t, int *verbose, double *weights)
{
  int k, m, n;
  double dista, mindist, e, e2, serror;
  double ermin,par1,par2;

  ermin = 0.0;
  serror = 0.0;
  e = e2 = 0.0;

  for(k=0; k<*xrows; k++){
    mindist=10e99;
    
    for(m=0; m<*ncenters; m++){
      
      
      dista=0;
      for(n=0;n<*xcols;n++){
	if(*dist == 0){
	    dista += pow((x[k+(*xrows)*n] - centers[m + (*ncenters)*n]), 2);
	}
	else if(*dist ==1){
          dista += fabs(x[k+(*xrows)*n] - centers[m + (*ncenters)*n]);
	}
      }    
      
      if(dista<mindist){
	cluster[k] = m;
	mindist = dista;
      }
    }
    /*NEW CENTERS*/
    /*----------------------------*/
    /*k-means learning rate */
    if (*methrate==0){
	par1=par[0];
	par2=par[1];
	t[cluster[k]]=t[cluster[k]]+1;
	e = pow(t[cluster[k]], -par1);
    }
    /*expon. decaying l.r*/
    else if(*methrate==1){
	par1=par[0];
	par2=par[1];
	e2=(double)*iter/(*itermax);
	e=par1 * pow(par2/par1,e2);  
    }
    for (n=0;n<*xcols;n++){
	centers[cluster[k]+(*ncenters)*n] +=
	    e*weights[k]*(x[k+(*xrows)*n]-centers[cluster[k]+(*ncenters)*n]);
    }
  }
  /*ERROR MINIMIZATION*/
  for (m=0;m<*ncenters;m++){
      for (k=0;k<*xrows;k++){
	  if (cluster[k]==m){
	      for(n=0;n<*xcols;n++){
		  if(*dist == 0){
		      serror +=
			  pow(x[k+(*xrows)*n] - centers[m+(*ncenters)*n], 2);
		  }
		  else if(*dist ==1){
		      serror +=
			  fabs(x[k+(*xrows)*n] - centers[m + (*ncenters)*n]);
		  }
	      }
	  }
      }
  }
  ermin=(serror)/(*xrows);
  if (*verbose){
      Rprintf("Iteration: %3d    Error:   %13.10f\n",*iter,ermin);
  }
  return 0; 
}

/**********************************************************/

int hardcl(int *xrows, int *xcols, double *x, int *ncenters,
	   double *centers, int *cluster,
	   int *itermax, int *iter, 
	   int *clustersize, int *verbose, int *dist, int *methrate,
	   double *par, double *weights)
{
  int m,k;
  int *t;
  

  t = (int *) R_alloc(*ncenters, sizeof(int));
  
  
  *iter=0;
  for (m=0;m<*ncenters;m++){
      t[m]=0;
  }
   
  while(((*iter)++ < *itermax)){
      oncent(xrows, xcols, x, ncenters, centers, cluster, clustersize,
	     dist,iter,itermax,methrate,par, t, verbose, weights);
  }
  for(k=0;k<*ncenters;k++){
      clustersize[k]=0;
  } 
  for(k=0; k<*xrows; k++){
      clustersize[cluster[k]]++;
  }
  return 0;
}
     
/**********************************************************/


int  oncentb(int *xrows, int *xcols, double *x, int *ncenters,
	   double *centers, int *cluster, int *clustersize,
	   int *dist,int *iter,int *itermax, double *par, int *verbose)
{
  int k, m, n, chang, a ,seira, minn;
  double e, h, l, aa, i,ermin,serror, mindist;
  double *dista;
  int  *ordd;

  dista = (double *) R_alloc(*ncenters, sizeof(double));
  ordd = (int *) R_alloc(*ncenters, sizeof(int));
  
  ermin=0.0;
  serror=0.0;


  for(k=0; k<*xrows; k++){
    
    for(m=0; m<*ncenters; m++){
         dista[m]=0.0;
    }
    
     for(m=0; m<*ncenters; m++){
      for(n=0;n<*xcols;n++){
	if(*dist == 0){
	  dista[m] += (x[k+(*xrows)*n] - centers[m +(*ncenters)*n])*(x[k+(*xrows)*n] - centers[m + (*ncenters)*n]); 
	}
	else if(*dist ==1){
          dista[m] += fabs(x[k+(*xrows)*n] - centers[m + (*ncenters)*n]);
	}
      }
     }
     
      /*ORDERING ACCORDING TO THE DISTANCE*/   
      for (m=0;m<*ncenters;m++){
	 ordd[m]=m; 
      }
      chang=1;
       while(chang!=0){
	     chang=0;
	for (m=0;m<(*ncenters-1);m++){
	
	  if (dista[m]> dista[m+1]){
	     aa=dista[m];
	     dista[m]=dista[m+1];
	     dista[m+1]=aa;
	    a=ordd[m];
	    ordd[m]=ordd[m+1];
	    ordd[m+1]=a;
	    chang=chang+1;
	  }
         
	}
       }
       
   
      /*NEW CENTERS*/
      for (m=0;m<*ncenters;m++){
	 seira=ordd[m];
	 /*printf("m:%d\n....ord:%d\n",m,seira  );*/
	  i=(double)(((*iter)-1)*(*xrows)+k)/((*itermax)*(*xrows));
          e=par[0]*pow(par[1]/par[0],i);
          l=par[2]*pow(par[3]/par[2],i);
          h=exp(-(double)m/l);
	  /*	  printf("par: %f,  %f, i:%f\n", par[0], par[2], i);
	  printf("m: %i, seira: %i\n", m,seira); */
	  for (n=0; n<*xcols;n++){
	    centers[seira+(*ncenters)*n]+=e*h*(x[k+(*xrows)*n]-centers[seira+(*ncenters)*n]); 
	  }
	  /*	  printf("\n");*/
      }
  }

  for (k=0;k<*xrows;k++){
      mindist=0.0;/*just to avoid compiling warnings*/
      minn=0; /*the same reason*/
      for (m=0;m<*ncenters;m++){
	  dista[m] = 0.0;
	  for(n=0;n<*xcols;n++){
	      if(*dist == 0){
		  dista[m] += (x[k+(*xrows)*n] -
			       centers[m+(*ncenters)*n])*
		      (x[k+(*xrows)*n] - centers[m + (*ncenters)*n]); 
	      }
	      else if(*dist ==1){
		  dista[m] += fabs(x[k+(*xrows)*n] -
				   centers[m + (*ncenters)*n]);
	      }
	  }
	  if (m == 0)
	  {
	      mindist = dista[0]; minn = 0;
	  }
	  else
	  {
	      if (dista[m] < mindist)
	      {
		  mindist = dista[m];
		  minn = m;
	      }
	  }
      }
      cluster[k] = minn;
  }
  
  
   /*ERROR MINIMIZATION*/
  for (m=0;m<*ncenters;m++){
    for (k=0;k<*xrows;k++){
      if (cluster[k]==m){
	 for(n=0;n<*xcols;n++){
       	if(*dist == 0){
	  serror += (x[k+(*xrows)*n] - centers[m
					      +(*ncenters)*n])*(x[k+(*xrows)*n] - centers[m +(*ncenters)*n]);                                        
	}
	else if(*dist ==1){
          serror += fabs(x[k+(*xrows)*n] - centers[m + (*ncenters)*n]);
	}
	 }
      }
    }
  }
  ermin=(serror)/(*xrows);
  
  /*if (*iter==1 | *iter==*itermax){*/
  if (*verbose){
      Rprintf("Iteration: %3d    Error:   %13.10f\n",*iter,ermin);
       }
  /*}*/
  /* printf("Iter:%d\n...",*iter);*/
 

  
    return 0;
}

/**********************************************************/

int neuralgas(int *xrows, int *xcols, double *x, int *ncenters,
	      double *centers, int *cluster,int *itermax, int *iter,
	      int *clustersize, int *verbose, int *dist,double *par)
{
    int k;
  
  *iter=0;
   
  while(((*iter)++ < *itermax)){

    oncentb(xrows, xcols, x, ncenters, centers, cluster, clustersize,
	   dist,iter,itermax,par,verbose);

     
  }
    for(k=0;k<*ncenters;k++){ 
     clustersize[k]=0;
      	 } 
  for(k=0; k<*xrows; k++){ 
  clustersize[cluster[k]]++; 
  /*printf("size%d\n...cluster%d\n",clustersize[cluster[k]],cluster[k]);*/
  }
  return 0;
}

