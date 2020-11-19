#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"

#define usage "\nNAME\n\
	tcl_evolve\n\
\n\
SYNOPSIS\n\
	tcl_evolve dt nsteps\n\
\n\
DESCRIPTION\n\
	tcl_evolve reads an lc catalog containing positions re[3], ra[3]\n\
	and velocities ve[3], va[3] from stdin, evolves the positions and\n\
	velocities through n steps of length dt using time centered\n\
	leapfrog, and writes the results to stdout.\n\
\n\
	Times are given in units of earths dynamical time (approx 58 days).\n\
	Distances are in AU.\n\
\n\
SEE ALSO\n\
	maketestpparfile.pl makeobs_circ laplace3\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"


main (int argc, char *argv[])
{
	double	*re, *ra, *ve, *va, *ge, *ga;
	double	dt;
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
	if (sscanf(argv[2], "%d", &nt) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}

	/* allocate space for re[3], ve[3], ra[3], va[3], ge[3], ga[3] */
	re = (double *) calloc(3, sizeof(double));
	ve = (double *) calloc(3, sizeof(double));
	ra = (double *) calloc(3, sizeof(double));
	va = (double *) calloc(3, sizeof(double));
	ge = (double *) calloc(3, sizeof(double));
	ga = (double *) calloc(3, sizeof(double));

	/* open lc pipe for input */
	ipf = popen("lc -b -o re ve ra va", "r");
	if (!ipf) {
		fprintf(stderr, "tcl_evolve : failed to open lc-pipe for input\n");
		exit(-1);
	}

	/* open lc-pipe for output */
	opf = popen("lc -C -N '1 3 re' -N '1 3 ve' -N '1 3 ra' -N '1 3 va'", "w");
	if (!opf) {
		fprintf(stderr, "tcl_evolve : failed to open lc-pipe for output\n");
		exit(-1);
	}

	while (1) {
		/* read initial conditions */
		if (!fread(re, sizeof(double), 3, ipf)) {
			pclose(ipf);
			pclose(opf);
			exit(0);
		}
		fread(ve, sizeof(double), 3, ipf);
		fread(ra, sizeof(double), 3, ipf);
		fread(va, sizeof(double), 3, ipf);
		for (it = 0; it < nt; it++) {
			/* lazy time centered algorithm */
			/* compute gravity */
			copy(ge, re);
			scale(ge, -1.0 / pow(length(re), 3.0));
			copy(ga, ra);
			scale(ga, -1.0 / pow(length(ra), 3.0));
			/* update velocities by half a step */
			for (i = 0; i < 3; i++) {
				ve[i] += 0.5 * dt * ge[i];
				va[i] += 0.5 * dt * ga[i];
			}
			/* update positions by a full step */
			for (i = 0; i < 3; i++) {
				re[i] += dt * ve[i];
				ra[i] += dt * va[i];
			}
			/* compute gravity */
			copy(ge, re);
			scale(ge, -1.0 / pow(length(re), 3.0));
			copy(ga, ra);
			scale(ga, -1.0 / pow(length(ra), 3.0));
			/* update velocities by second half step */
			for (i = 0; i < 3; i++) {
				ve[i] += 0.5 * dt * ge[i];
				va[i] += 0.5 * dt * ga[i];
			}
		}
		fprintvec(re, opf);
		fprintvec(ve, opf);
		fprintvec(ra, opf);
		fprintvec(va, opf);
	}
}

