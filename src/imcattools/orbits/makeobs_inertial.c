#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"
#include "Ffunc.h"
#include "orbitutils/gaussdev.h"

#define usage "\nNAME\n\
	makeobs_inertial - make a set of observations for an inertial observatory\n\
\n\
SYNOPSIS\n\
	makeobs_inertial dt sigma\n\
\n\
DESCRIPTION\n\
	makeobs_inertial reads an lc format catalog from stdin containing at least\n\
	re[3], ve[3], ra[3], va[3], these being the position and velocity of\n\
	the earth and the asteroid at t=0. It then generates a catalog\n\
	containing direction of the asteroid at times t = -dt, 0, +dt\n\
	and those times.  The catalog also contains, the position of the\n\
	earth re[] and the observatory rho[] wrt the earth at those times.\n\
\n\
	Times are given in units of earths dynamical time (approx 58 days).\n\
	Distances are in AU.\n\
\n\
SEE ALSO\n\
	maketestpscoords.pl laplace3\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"


main (int argc, char *argv[])
{
	double	*re0, *ve0, *ra0, *va0, *dr, *re, *ra, *ge, *ga, *rho, *n;
	double	*dddotre, *dddotra;
	double	sigma, dt, *t;
	int	i, it, nt;
	FILE	*ipf, *opf;

	/* parse args */
	if (argc != 3) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[1], "%lf", &dt) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[2], "%lf", &sigma) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}

	/* create the observation times */
	nt = 3;
	t = (double *) calloc(nt, sizeof(double));
	for (it = 0; it < nt; it++) {
		t[it] = dt * (it - 1);
	}

	/* allocate space for re[3], ve[3], ra[3], va[3] */
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
	/* gravity */
	ge = (double *) calloc(3, sizeof(double));
	ga = (double *) calloc(3, sizeof(double));
	/* d^3 r / d t^3 */
	dddotre = (double *) calloc(3, sizeof(double));
	dddotra = (double *) calloc(3, sizeof(double));
	/* observatory position */
	rho = (double *) calloc(3, sizeof(double));

	/* read the phase-space coords */
	ipf = popen("lc -b -o re0 ve0 ra0 va0", "r");
	if (!ipf) {
		fprintf(stderr, "makeobs_inertial : failed to open lc-pipe for input\n");
		exit(-1);
	}
	fread(re0, sizeof(double), 3, ipf);
	fread(ve0, sizeof(double), 3, ipf);
	fread(ra0, sizeof(double), 3, ipf);
	fread(va0, sizeof(double), 3, ipf);
	pclose(ipf);

	/* compute gravity */
	copy(ge, re0);
	scale(ge, -1.0 / pow(length(re0), 3.0));
	copy(ga, ra0);
	scale(ga, -1.0 / pow(length(ra0), 3.0));

	/* compute d^3 r / d t^3 */
	makedddotr(dddotre, re0, ve0);
	makedddotr(dddotra, ra0, va0);

	/* open lc-pipe for output */
	opf = popen("lc -C -n t -N '1 3 re' -N '1 3 rho' -N '1 3 n' -n sigma", "w");
	if (!opf) {
		fprintf(stderr, "makeobs_inertial : failed to open lc-pipe for output\n");
		exit(-1);
	}

	/* compute observations */
	for (it = 0; it < nt; it++) {
		/* compute the position to quadratic order */
		for (i = 0; i < 3; i++) {
			re[i] = re0[i] + ve0[i] * t[it] + 0.5 * ge[i] * t[it] * t[it];
			ra[i] = ra0[i] + va0[i] * t[it] + 0.5 * ga[i] * t[it] * t[it];
		}
		/* add the cubic correction */
		for (i = 0; i < 3; i++) {
			re[i] += dddotre[i] * pow(t[it], 3.0) / 6.0;
			ra[i] += dddotra[i] * pow(t[it], 3.0) / 6.0;
		}
		diff(n, ra, re);
		scale(n, 1.0 / length(n));
		/* add error */
		for (i = 0; i < 3; i++) {
			n[i] += gaussdev() * sigma;
		}
		scale(n, 1.0 / length(n));
		fprintf(opf, "%14.8lg\n", t[it]);
		fprintvec(re, opf);
		fprintvec(rho, opf);
		fprintvec(n, opf);
		fprintf(opf, "%14.8lg\n", sigma);
	}
	pclose(opf);

	exit(0);
}

