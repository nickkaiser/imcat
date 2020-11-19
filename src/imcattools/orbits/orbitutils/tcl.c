#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"
#include "planetdata.h"
#include "tcl.h"

static	int	dojupiter = 0;
static	double	jupiterphase, jupiteromega, jupitermass, *jupiterr;

void	tcl(double dt, int nsteps, double*r, double *v)
{
	int	it, i;
	double	*g;

	/* allocate space for ge */
	g = (double *) calloc(3, sizeof(double));

	for (it = 0; it < nsteps; it++) {
		/* lazy time centered algorithm */
		/* compute gravity */
		copy(g, r);
		scale(g, -1.0 / pow(length(r), 3.0));
		if (dojupiter) {
			addjupitergravity(g, r);
		}
		/* update velocities by half a step */
		for (i = 0; i < 3; i++) {
			v[i] += 0.5 * dt * g[i];
		}
		/* update positions by a full step */
		for (i = 0; i < 3; i++) {
			r[i] += dt * v[i];
		}
		/* compute gravity */
		copy(g, r);
		scale(g, -1.0 / pow(length(r), 3.0));
		if (dojupiter) {
			addjupitergravity(g, r);
			jupiterphase += jupiteromega * dt;
		}
		/* update velocities by second half step */
		for (i = 0; i < 3; i++) {
			v[i] += 0.5 * dt * g[i];
		}
	}
	free(g);
}

void	jupiteron(double phase)
{
	dojupiter = 1;
	jupiterphase = phase;
	jupiteromega = pow(A_JUPITER, -1.5);
	jupiterr = (double *) calloc(3, sizeof(double));
	jupitermass = M_JUPITER / M_SUN;
}

void	addjupitergravity(double *g, double *r)
{
	int	i;

	assign(jupiterr, A_JUPITER * cos(jupiterphase) - r[0], A_JUPITER * sin(jupiterphase) - r[1], -r[2]);
	scale(jupiterr, jupitermass / pow(length(jupiterr), 3.0));
	for (i = 0; i < 3; i++) {
		g[i] += jupiterr[i];
	}
}
