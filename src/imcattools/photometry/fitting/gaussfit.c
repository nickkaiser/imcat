/**
 ** function for fitting a Gaussian ellipsiod
 **/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../../utils/stats_stuff.h"
#include "../../../utils/frprmn.h"
#include "../../../utils/fitstatus.h"
#include "../../../utils/powell.h"
#include "../../../utils/nrutil.h"
#include "gaussfit.h"



static int			n1, n2;
static float 			**f;

int     gaussfit(float **ff, int N1, int N2, float *f0, float *d, float *a, float *b, float *phi)
{
	float 	*p, ftol, fret;
	int	i, iter;
	float	**xi;
	
	ftol = 1.e-20;
	f = ff;				/* set global base address for f[][] */
	n1 = N1; n2 = N2;		/* set global image size */

	/* starting guess */
	p = vector(1, 6);
	p[1] = (*f0 != 0.0 ? *f0 : f[N2/2][N1/2]);
	p[2] = (d[0] != 0.0 ? d[0] : 0.5 * N1);
	p[3] = (d[1] != 0.0 ? d[1] : 0.5 * N2);
	p[4] = *a;
	p[5] = *b;
	p[6] = *phi;
	if (p[4] == p[6])
		p[6] = 1.03 * p[4];

	/* starting directions... */
	xi = matrix(1,6,1,6);
	for (i = 1; i <= 6; i++) {
		xi[i][i] = 1.0;
	}

	powell(p, xi, 6, ftol, &iter, &fret, func);

/*
	frprmn(p, 6, ftol, &iter, &fret, func, dfunc);
*/

	*f0 = p[1];
	d[0] = p[2];
	d[1] = p[3];
	*a = p[4];
	*b = p[5];
	*phi = p[6];

	free_vector(p, 1, 6);
	free_matrix(xi, 1, 6, 1, 6);
}



float	func(float *p)
{
	int 	ix, iy;
	float 	chi2 = 0.0, fmodel, f0, x, y, a, b, phi, xp, yp, c, s;

	f0 = p[1];
	x = p[2];
	y = p[3];
	a = p[4];
	b = p[5];
	phi = p[6];

	c = cos(phi);
	s = sin(phi);
	for (iy = 0; iy < n2; iy++) {	
		for (ix = 0; ix < n1; ix++) {
			xp = (ix - x) * c + (iy - y) * s;
			yp = (iy - y) * c - (ix - x) * s;
			fmodel = f0 * exp(-0.5 * (xp * xp / (a * a) + yp * yp / (b * b)));
			chi2 += (f[iy][ix] - fmodel) * (f[iy][ix] - fmodel);
		}
	}
	return (chi2);
}


void	makemodel(float **ff, int N1, int N2, float f0, float *d, float a, float b, float phi)
{
        int     ix, iy;
        float   x, y, xp, yp, c, s;

        x = d[0];
        y = d[1];

	c = cos(phi);
	s = sin(phi);
       	for (iy = 0; iy < n2; iy++) {
                for (ix = 0; ix < n1; ix++) {
 			xp = (ix - x) * c + (iy - y) * s;
			yp = (iy - y) * c - (ix - x) * s;
			ff[iy][ix] = f0 * exp(-0.5 * (xp * xp / (a * a) + yp * yp / (b * b)));
                 }
        }

}
