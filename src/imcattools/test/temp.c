#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "amoeba_dbl.h"
#include "amotry_dbl.h"
#include "../../utils/nrutil.h"

float	myfunc(float *x);

#define NDIM 3

main(int argc, char * argv[])
{
	float	*x, **p, *y, ftol;
	int	nfunk, i, j;

	p = matrix(1, NDIM + 1, 1, NDIM);
	for (i = 1; i <= NDIM + 1; i++) {
		for (j = 1; j <= NDIM; j++) {
			p[i][j] = (float) drand48();
		}
	}
	y = vector(1, NDIM + 1);
	for (i = 1; i <= NDIM + 1; i++) {
		y[i] = myfunc(p[i]);
	}
	ftol = 1.e-20;
	amoeba(p, y, NDIM, ftol, myfunc, &nfunk);
	fprintf(stdout, "%d function calls\n", nfunk);
	for (i = 1; i <= NDIM + 1; i++) {
		for (j = 1; j <= NDIM; j++) {
			fprintf(stdout, "%13g ", p[i][j]);
		}
		fprintf(stdout, "\n");
	}
}


float	myfunc(float *x)
{
	float	res = 0.0;
	int	dim;

	for (dim = 1; dim <= NDIM; dim++) {
		res += sqrt((x[dim] - 1) * (x[dim] - 1));
	}
	return (res);
}
