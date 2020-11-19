#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "amoeba_dbl.h"
#include "amotry_dbl.h"
#include "../../utils/nrutil.h"

double	myfunc(double *x);

#define NDIM 3

main(int argc, char * argv[])
{
	double	*x, **p, *y, ftol;
	int	nfunk, i, j;

	p = (double **) calloc(NDIM + 1, sizeof(double *));
	for (i = 0; i < NDIM + 1; i++) {
		p[i] = (double *) calloc(NDIM, sizeof(double));
	}
	for (i = 0; i < NDIM + 1; i++) {
		for (j = 0; j < NDIM; j++) {
			p[i][j] = (double) drand48();
		}
	}
	y = (double *) calloc(NDIM + 1, sizeof(double));
	for (i = 0; i < NDIM + 1; i++) {
		y[i] = myfunc(p[i]);
	}
	ftol = 1.e-20;
	amoeba(p, y, NDIM, ftol, myfunc, &nfunk);
	fprintf(stdout, "%d function calls\n", nfunk);
	for (i = 0; i < NDIM + 1; i++) {
		for (j = 0; j < NDIM; j++) {
			fprintf(stdout, "%13g ", p[i][j]);
		}
		fprintf(stdout, "\n");
	}
}


double	myfunc(double *x)
{
	double	res = 0.0;
	int	dim;

	for (dim = 0; dim < NDIM; dim++) {
		res += sqrt((x[dim] - 1) * (x[dim] - 1));
	}
	return (res);
}
