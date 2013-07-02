#include <R.h>

void predstr_calcSums(int pred[], int tc[], int *length, int *levels, int res[]){

  int mat = 0;
  int i, j, k;
  int co[*length];
  
  for(i = 0; i < *length; i++){
    for(j = i + 1; j < *length; j++){
      for(k = 1; k <= *levels; k++){
        if(pred[i] == pred[j])
          mat = 2;
        else 
          mat = 0;
        if(tc[j] == k && tc[i] == k)
          res[k-1] = res[k-1] + mat;
      }
    }
  }
}

void countPairs(int *c1, int *c2, int *n, double *res)
{
    int i, j;
    double n01=0, n10=0, n11=0;

    for(i=0; i<*n; i++){
	for(j=i+1; j<*n; j++){
	    if(c1[i]==c1[j]){
		if(c2[i]==c2[j])
		    n11++;
		else
		    n10++;
	    } else {
		if(c2[i]==c2[j]) n01++;
	    }
	}
    }
    res[0]=(*n)*(*n-1)/2 - n11 - n10 - n01;
    res[1]=n10;
    res[2]=n01;
    res[3]=n11;
}
		
