/**
 ** function for fitting a beta model
 **/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../../utils/stats_stuff.h"
#include "../../../utils/powell.h"
#include "../../../utils/nrutil.h"
#include "betafit.h"
#include "amoeba.h"

static int			n1, n2;
static float 			**gn;
static float			gN;
static float			gbeta;


int     fitall(float **n, int N1, int N2, float *f0, float *x0, float *y0, float *rc2, float beta, float *loglhood, float N)
{
	float 	*p, **xi, ftol, fret;
	int	i, j, iter;
	
	ftol = 1.e-20;
	gn = n;				/* set global base address for n[][] */
	n1 = N1; n2 = N2;		/* set global image size */
	gN = N;				/* set global number of micropix per pixels */
	gbeta = beta;

	/* set up the starting simplex */
	p = vector(1,4);
	p[1] = log(*f0);
	p[2] = *x0;
	p[3] = *y0;
	p[4] = log(*rc2);
	
	/* and the directions */
	xi = matrix(1,4,1,4);
	for (i = 1; i <= 4; i++) {
		for (j = 1; j <= 4; j++) {
			xi[i][j] = 0.0;
		}
		xi[i][i] = 1;
	}


	powell(p, xi, 4, ftol, &iter, &fret, func);

	*f0   = exp(p[1]);
	*x0   = p[2];
	*y0   = p[3];
	*rc2  = exp(p[4]);
	*loglhood = func(p);

	free_vector(p, 1, 4);
	free_matrix(xi, 1, 4, 1, 4);
}





float	func(float *p)
{
	int 	ix, iy;
	float 	lhood, fmodel, x0, y0,  f0, rc2, beta;

	f0   = exp(p[1]);
	x0   = p[2];
	y0   = p[3];
	rc2  = exp(p[4]);

/* fprintf(stdout, "%13.8g %13.8g %13.8g %13.8g %13.8g\n", f0, x0, y0, rc2, gbeta); */
	lhood = 0; 
	for (iy = 0; iy < n2; iy++) {	
		for (ix = 0; ix < n1; ix++) {
			fmodel = f0 * pow(1 + ((ix - x0) * (ix - x0) + (iy - y0) * (iy - y0)) / rc2, -gbeta);
			lhood += gN * fmodel;
			if (gn[iy][ix] > 0.0) {
				lhood -= log(fmodel) * gn[iy][ix];
			}
		}
	}
	return (lhood);
}

void	makebetamodel(float **nmodel, int N1, int N2, float f0, float x0, float y0, float rc2, float beta, float N)
{
        int     ix, iy;

 
        for (iy = 0; iy < N2; iy++) {
                for (ix = 0; ix < N1; ix++) {
			nmodel[iy][ix] = N * f0 * pow(1 + ((ix - x0) * (ix - x0) + (iy - y0) * (iy - y0)) / rc2, -beta);
                 }
        }

}
