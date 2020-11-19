#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../utils/nrutil.h"

void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret,
	double (*func)(double []));
double	myfunc(double *x);

#define NDIM 3

main(int argc, char * argv[])
{
	double	**xi, *p, *y, ftol, fret;
	int	iter, i;

	p = (double *) calloc(NDIM, sizeof(double));
	xi = (double **) calloc(NDIM, sizeof(double *));
	for (i = 0; i < NDIM; i++) {
		xi[i] = (double *) calloc(NDIM, sizeof(double));
		xi[i][i] = 1.0;
	}
	ftol = 1.e-20;
	powell(p, xi, NDIM, ftol, &iter, &fret, myfunc);
	for (i = 0; i < NDIM; i++) {
		fprintf(stdout, "%13g ", p[i]);
	}
	fprintf(stdout, "\n");
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
