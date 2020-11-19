#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"
#include "Ffunc.h"

#define USAGE "usage : orbit\n"

#define	DTINY	1.0e-20
#define FTOL	1.0e-8
#define NITERMAX	100


main (int argc, char *argv[])
{
	double	*re0, *ve0, *ra0, *va0, *dr, *n, *re, *ra, *va, *ge, *ga, **nobs, *tmp, *reperp;
	double	*dotn, *ddotn, *ddotnperp, *dotnhat;
	double	mu, rre, D, dt, *t, d, dd, len1, len2, ltmp, Fobs, ndot, sigma;
	int	nt, it, i, niter;

	/* observational error */
	sigma = 0.0;

	/* set up observation times and directions */
	nt 	= 3;
	dt 	= 0.02;
	t  	= (double *) calloc(nt, sizeof(double));
	nobs  	= (double **) calloc(nt, sizeof(double *));
	for (it = 0; it < nt; it++) {
		t[it] = dt * (it - 1);
		nobs[it] = (double *) calloc(3, sizeof(double));
	}
	
	/* allocate initial position and velocities */
	re0 = (double *) calloc(3, sizeof(double));
	ve0 = (double *) calloc(3, sizeof(double));
	ra0 = (double *) calloc(3, sizeof(double));
	va0 = (double *) calloc(3, sizeof(double));
	/* dr = ra - re */
	dr = (double *) calloc(3, sizeof(double));
	/* n = |d| */
	n = (double *) calloc(3, sizeof(double));
	/* actual positions */
	re = (double *) calloc(3, sizeof(double));
	ra = (double *) calloc(3, sizeof(double));
	va = (double *) calloc(3, sizeof(double));
	/* gravity */
	ge = (double *) calloc(3, sizeof(double));
	ga = (double *) calloc(3, sizeof(double));
	/* direction derivatives */
	dotn = (double *) calloc(3, sizeof(double));
	ddotn = (double *) calloc(3, sizeof(double));
	ddotnperp = (double *) calloc(3, sizeof(double));
	dotnhat = (double *) calloc(3, sizeof(double));
	/* a temporary vector */
	tmp = (double *) calloc(3, sizeof(double));
	/* re perpendicular */
	reperp = (double *) calloc(3, sizeof(double));
	

	/* assign initial positions and velocities */
	assign(re0, 1.0, 0.0, 0.0);
	assign(ve0, 0.0, 1.0, 0.0);
	assign(ra0, 2.0, 1.0, 0.3);
	assign(va0, 0.0, 0.0, 0.5);

	/* print initial positions and velocities */
	printvec(re0, "re0");
	printvec(ve0, "ve0");
	printvec(ra0, "ra0");
	printvec(va0, "va0");

	/* compute gravity */
	copy(ge, re0);
	scale(ge, -1.0 / pow(length(re0), 3.0));
	copy(ga, ra0);
	scale(ga, -1.0 / pow(length(ra0), 3.0));
	printvec(ge, "ge");
	printvec(ga, "ga");

	/* compute observations */
	for (it = 0; it < nt; it++) {
		/* compute the position */
		for (i = 0; i < 3; i++) {
			re[i] = re0[i] + ve0[i] * t[it] + 0.5 * ge[i] * t[it] * t[it];
			ra[i] = ra0[i] + va0[i] * t[it] + 0.5 * ga[i] * t[it] * t[it];
		}
		fprintf(stdout, "t=%10.3f\n", t[it]);
		diff(dr, ra, re);
		copy(nobs[it], dr);
		scale(nobs[it], 1.0 / length(nobs[it]));
		/* add error */
		for (i = 0; i < 3; i++) {
			nobs[it][i] += (2 * drand48() - 1) * sigma;
		}
		scale(nobs[it], 1.0 / length(nobs[it]));
		printvec(nobs[it], "nobs");
	}
	fprintf(stdout, "\n");

	/* get estimates of n[], dotn[], ddotn[] */
	/* we should solve A N = B here, but for now we use the solution for nt=3 equal interval case */
	copy(n, nobs[1]);
	diff(dotn, nobs[2], nobs[0]);
	scale(dotn, 0.5 / dt);
	copy(dotnhat, dotn);
	scale(dotnhat, 1.0 / length(dotn));
	for (i = 0; i < 3; i++) {
		ddotn[i] = (nobs[2][i] - 2 * nobs[1][i] + nobs[0][i]) / (dt * dt);
	}
	getperp(ddotnperp, ddotn, n);
	printvec(n, "n");
	printvec(dotn, "dotn");
	printvec(dotnhat, "dotnhat");
	printvec(ddotn, "ddotn");
	printvec(ddotnperp, "ddotnperp");

	/* get mu */
	mu = dot(re0, n) / length(re0);
	fprintf(stdout, "mu = %14.8lf\n", mu);

	/* compute Fobs */
	getperp(tmp, ddotnperp, dotnhat);
	len1 = length(tmp);
	getperp(reperp, re0, n);
	printvec(reperp, "reperp");
	getperp(tmp, reperp, dotnhat);
	len2 = length(tmp);
	if (fabs(len2) > DTINY) {
		Fobs = len1 * pow(length(re0), 4.0) / len2;
		fprintf(stdout, "Fobs = %14.8lf\n", Fobs);
	}

	diff(dr, ra, re);
	D = length(dr);
	rre = length(re0);
	d = D / rre;
	fprintf(stdout, "d=%14.8lf\nF=%14.8lf\n", d, F(d, mu));

	/* solve for d */
	fprintf(stdout, "solving for d...\n");
	d *= 2.0;	/* initial d - make it twice the true value to exercise the solver */
	niter = 0;
	while (fabs(Fobs - F(d, mu)) > FTOL) { 
		d += (Fobs - F(d, mu)) / dFdd(d, mu);
		niter++;
		if (niter > NITERMAX) {
			fprintf(stderr, "orbit : niter > NITERMAX\n");
			exit(-1);
		}
	}
	fprintf(stdout, "d=%14.8lf\nF=%14.8lf\nniter = %d\n", d, F(d, mu), niter);
	
	/* compute the position and velocity */
	D = d * rre;
	ndot = length(dotn);
	for (i = 0; i < 3; i++) {
		ra[i] = re0[i] + D * n[i];
		tmp[i] = ddotnperp[i] - Fobs * reperp[i] / pow(rre, 4.0);
	}
	ltmp = length(tmp);
	for (i = 0; i < 3; i++) {
		va[i] = ve0[i] + D * dotn[i] - D * ltmp * n[i] / (2.0 * ndot);
	}
	printvec(ra, "ra");
	printvec(va, "va");

	exit(0);
}

