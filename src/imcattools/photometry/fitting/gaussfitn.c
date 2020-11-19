/**
 ** function for fitting a Gaussian ellipsiod
 **/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../../utils/arrays.h"
#include "../../../utils/stats_stuff.h"
#include "../../../utils/frprmn.h"
#include "../../../utils/fitstatus.h"
#include "../../../utils/powell.h"
#include "../../../utils/nrutil.h"
#include "gaussfitn.h"



static int			n1, n2, nn;
static float 			**f, **fmodel;

int     gaussfitn(float **ff, int N1, int N2, float *f0, float *x, float *y, float *a, float *b, float *phi, int n, float **fmod)
{
	float 	*p, ftol, fret;
	int	i, iter;
	float	**xi;
	
	/* tolerance */
	ftol = 1.e-20;

	/* set globals */
	f = ff;	
	n1 = N1; 
	n2 = N2;
	nn = n;
	fmodel = fmod;
	
	/* create and initialise p-vector */
	p = vector(1, 6 * n);
	for (i = 0; i < n; i++) {
		p[6 * i + 1] = f0[i];
		p[6 * i + 2] = x[i];
		p[6 * i + 3] = y[i];
		p[6 * i + 4] = a[i];
		p[6 * i + 5] = b[i];
		p[6 * i + 6] = phi[i];
	}

	/* starting directions... */
	xi = matrix(1, 6 * n, 1, 6 * n);
	for (i = 1; i <= 6 * n; i++) {
		xi[i][i] = 1.0;
	}

	powell(p, xi, 6 * n, ftol, &iter, &fret, func);

	for (i = 0; i < n; i++) {	
		f0[i] 	= p[6 * i + 1];
		x[i] 	= p[6 * i + 2];
		y[i] 	= p[6 * i + 3];
		a[i] 	= p[6 * i + 4];
		b[i] 	= p[6 * i + 5];
		phi[i] 	= p[6 * i + 6];
	}

	/* one more call to make sure we have solution fmodel */
	func(p);

	/* clean up */
	free_vector(p, 1, 6 * n);
	free_matrix(xi, 1, 6 * n, 1, 6 * n);
}



float	func(float *p)
{
	int 	i, ix, iy;
	float 	chi2, f0, x, y, a, b, phi, xp, yp, c, s;
	
	/* compute the model */
	for (iy = 0; iy < n2; iy++) {	
		for (ix = 0; ix < n1; ix++) {
			fmodel[iy][ix] = 0.0;
		}
	}
	for (i = 0; i < nn; i++) {	
		f0  = p[6 * i + 1];
		x   = p[6 * i + 2];
		y   = p[6 * i + 3];
		a   = p[6 * i + 4];
		b   = p[6 * i + 5];
		phi = p[6 * i + 6];
		c = cos(phi);
		s = sin(phi);
	 	for (iy = 0; iy < n2; iy++) {	
			for (ix = 0; ix < n1; ix++) {
				xp = (ix - x) * c + (iy - y) * s;
				yp = (iy - y) * c - (ix - x) * s;
				fmodel[iy][ix] += f0 * exp(-0.5 * (xp * xp / (a * a) + yp * yp / (b * b)));
			}
		}
	}
	chi2 = 0.0;
	for (iy = 0; iy < n2; iy++) {	
		for (ix = 0; ix < n1; ix++) {
			chi2 += (f[iy][ix] - fmodel[iy][ix]) * (f[iy][ix] - fmodel[iy][ix]);
		}
	}
	return (chi2);
}


