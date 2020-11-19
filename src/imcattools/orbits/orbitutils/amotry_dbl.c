/*
 * amotry_dbl.c
 *
 * original code from NR converted to double and zero-based indexing by Nick Kaiser
 */


#define NRANSI
#include "nrutil.h"

double amotry(double **p, double y[], double psum[], int ndim,
	double (*funk)(double []), int ihi, double fac)
{
	int j;
	double fac1,fac2,ytry,*ptry;

	ptry=(double *) calloc(ndim, sizeof(double));
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=0;j<ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	ytry=(*funk)(ptry);
	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (j=0;j<ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	free(ptry);
	return ytry;
}
#undef NRANSI
