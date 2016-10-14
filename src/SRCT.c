//#include "stdafx.h"
#include "math.h"
#include "stdlib.h" 
#include "R.h"
#include "Rmath.h"
#include "stdio.h"


// double **alloc_matrix(int r, int c)
// {
    // /* allocate a matrix with r rows and c columns */
    // int i;
    // double **matrix;
    // matrix = (double **) calloc(r, sizeof(double *));
    // for (i = 0; i < r; i++)
    // matrix[i] = (double *) calloc(c, sizeof(double));
    // return matrix;
// }


// void free_matrix(double **matrix, int r, int c)
// {
    // /* free a matrix with r rows and c columns */
    // int i;
    // for (i = 0; i < r; i++){
		// free(matrix[i]);	
	// }
    // free(matrix);
// }


// void vector2matrix(double *x, double **y, int N, int d, int isroworder) {
    // /* copy a d-variate sample into a matrix, N samples in rows */
    // int i, k;
    // if (isroworder == TRUE) {
        // for (k=0; k<d; k++)
            // for (i=0; i<N; i++)
                // y[i][k] = (*(x+i*d+k));
        // }
    // else {
        // for (k=0; k<N; k++)
            // for (i=0; i<d; i++)
                // y[i][k] = (*(x+k*N+i));
        // }
    // return;
// }

// void Euclidean_distance(double *x, double **Dx, int n, int d)
// {
    // /*
        // interpret x as an n by d matrix, in row order (n vectors in R^d)
        // compute the Euclidean distance matrix Dx
    // */
    // int i, j, k, p, q;
    // double dsum, dif;
    // for (i=1; i<n; i++) {
        // Dx[i][i] = 0.0;
        // p = i*d;
        // for (j=0; j<i; j++) {
            // dsum = 0.0;
            // q = j*d;
            // for (k=0; k<d; k++) {
                // dif = *(x+p+k) - *(x+q+k);
                // dsum += dif*dif;
            // }
            // Dx[i][j] = Dx[j][i] = sqrt(dsum);
        // }
    // }
// }


// void distance(double *x, double *Dx, int *n, int *d)
// {
    /*
        interpret x as an n by d matrix, in row order (n vectors in R^d)
        compute the Euclidean distance matrix Dx
    */
    // int i, j, k, p, q;
    // double dsum, dif;
    // for (i=1; i<n[0]; i++) {
        // p = i*d[0];
        // for (j=0; j<i; j++) {
            // dsum = 0.0;
            // q = j*d[0];
            // for (k=0; k<d[0]; k++) {
                // dif = *(x+p+k) - *(x+q+k);
                // dsum += dif*dif;
            // }
            // Dx[i*n[0]+j] = Dx[j*n[0]+i] = sqrt(dsum);
        // }
    // }
// }

void SRCT(double *x, double *t, double *delta, double *Sc, int *n, double *RCTV) 
{
	/*  computes RCT(x,y)  */	

    int    i, j;
    double jp=0, p1 = 0, p2 = 0;
    
    *RCTV = 0;

	/*  I need to change  */	
    for(i=0;i<(*n-1);i++){
		for(j=0;j<(*n);j++){
            if((x[j]>x[i]) && (t[j]>t[i]))
			    jp += 1;
            if(x[j]>x[i])
                p1 += 1;				
            if(t[j]>t[i])
                p2 += 1;
		}
		
		*RCTV += pow(jp/(*n)-p1*p2/((*n)*(*n)),2)*delta[i]/pow(Sc[i],3);
		jp = 0;
		p1 = 0;
		p2 = 0;
	}
    *RCTV = *RCTV/(1.0*(*n));	
      
  
	return;
}	 



