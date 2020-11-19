/**
 ** function for fitting a generalised Gaussian ellipsiod
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
#include "gengaussfit.h"



static float			galpha, gbeta;
static int			n1, n2;
static float 			**f;

/* 
 * gaussfit()
 * 	fits a generalised Gaussian ellipsoid 
 *		f0 exp(-[q_ij(r-d)_i(r-d)_j]^alpha)
 * 	to N1 * N2 image f by minimising
 *		sum (f - f_model)^beta
 *	use f0, a, b to generate starting values (if zero then we use sensible guess)
 */
int	gaussfit(float **ff, int N1, int N2, float *d, float **q, float *f0,
	float alpha, float beta)
{
	float 	*p, ftol, fret;
	int	i, iter;
	float	Cxx, Cyy, Cxy, Cp, Cm, lambda, ab, D;
	float	**xi;
	
	ftol = 1.e-20;
	f = ff;				/* set global base address for f[][] */
	n1 = N1; n2 = N2;		/* set global image size */

	galpha = alpha;
	gbeta = beta;

	/* starting guess */
	p = vector(1, 6);
	p[1] = (*f0 != 0.0 ? *f0 : f[N2/2][N1/2]);
	p[2] = (d[0] != 0.0 ? d[0] : 0.5 * N1);
	p[3] = (d[1] != 0.0 ? d[1] : 0.5 * N2);
	p[4] = (q[0][0] != 0.0 ? q[0][0] : 0.5);
	p[5] = (q[0][1] != 0.0 ? q[0][1] : 0.0);
	p[6] = (q[1][1] != 0.0 ? q[1][1] : 0.51);
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
	q[0][0] = p[4];
	q[0][1] = q[1][0] = p[5];
	q[1][1] = p[6];

	free_vector(p, 1, 6);
	free_matrix(xi, 1, 6, 1, 6);
}



float	func(float *p)
{
	int 	ix, iy;
	float 	value = 0, f0, x, y, Cxx, Cxy, Cyy, expfunc, quad;

	f0 = p[1];
	x = p[2];
	y = p[3];
	Cxx = p[4];
	Cxy = p[5];
	Cyy = p[6];

	for (iy = 0; iy < n2; iy++) {	
		for (ix = 0; ix < n1; ix++) {
			quad = ( 
                                Cxx * (ix - x) * (ix - x)
                                + 2 * Cxy * (ix - x) * (iy - y)
                                + Cyy * (iy - y) * (iy - y)
                        );
			expfunc = exp(-0.5 * pow(fabs(quad), galpha));
			value += pow(fabs(f[iy][ix] - f0 * expfunc), gbeta);
		}
	}
	return (value);
}


void    makemodel(float **ff, int N1, int N2, float *d, float **q, float f0)
{
        int     ix, iy;
        float   x, y, Cxx, Cxy, Cyy, expfunc, quad;

        x = d[0];
        y = d[1];
        Cxx = q[0][0];
        Cxy = q[0][1];
        Cyy = q[1][1];

        for (iy = 0; iy < n2; iy++) {
                for (ix = 0; ix < n1; ix++) {
                        quad = (
                                Cxx * (ix - x) * (ix - x)
                                + 2 * Cxy * (ix - x) * (iy - y)
                                + Cyy * (iy - y) * (iy - y)
                        );
                        expfunc = exp(-0.5 * pow(fabs(quad), galpha));
                        ff[iy][ix] = f0 * expfunc;
                }
        }

}
