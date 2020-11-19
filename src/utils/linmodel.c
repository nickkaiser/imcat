/*
 * linmodel.c
 *
 * functions to fit a stream of data as a linear combination of parameters
 *
 * see linmodel.h for info
 */
 
#include	<math.h>
#include	"error.h"
#include	"lu.h"
 
 static	double 	**A, *B;
 static	int	np, *indx;
 
 void	linmodelinit(int NP)
 {
 	int	par;
	
 	np = NP;
	indx = (int *) calloc(np, sizeof(int));
	B = (double *) calloc(np, sizeof(double));
	A = (double **) calloc(np, sizeof(double *));
	for (par = 0; par < np; par++)
		A[par] = (double *) calloc(np, sizeof(double));
 }
 
 
 
 
 void	linmodelincrement(double f, double *p)
 {
 	int	par1, par2;
	
	for (par1 = 0; par1 < np; par1++) {
		B[par1] += f * p[par1];
		for (par2 = 0; par2 < np; par2++)
			A[par1][par2] += p[par1] * p[par2];
	}
 }
 
 
 
 
 void	linmodelsolve(double *F)
 {
 	int	par;
	double	det;
	
	myludcmp(A, np, indx, &det);
	mylubksb(A, np, indx, B);
	
	for (par = 0; par < np; par++)
		F[par] = B[par];
		
 	free(B);
	free(indx);
	for (par = 0; par < np; par++)
		free(A[par]);
	free(A);
 }
 

