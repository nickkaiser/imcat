#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"

#define usage "\nNAME\n\
	r2n - get perpendicular components of a vector\n\
\n\
SYNOPSIS\n\
	r2n rname x y z\n\
\n\
DESCRIPTION\n\
	r2n reads an lc catalog containing positions rname[3]\n\
	from stdin, evolves the positions and outputs a two-vector\n\
	n[2] which are the components of the unit vector |rname|\n\
	perpendicular to the direction (x, y, z).\n\
\n\
	The reference vector (x, y, z) need not be normalized.\n\
\n\
SEE ALSO\n\
	maketestpparfile.pl makeobs_circ laplace3 tcl_evolve\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"


main (int argc, char *argv[])
{
	double	*r, *rref, *i0, *i1, n0, n1;
	double	x, y, z, theta, phi;
	char	*rname, lccom[1024];
	int	i, it, nt;
	FILE	*ipf, *opf;

	/* parse args */
	if (argc != 5) {
		fprintf(stderr, usage);
		exit(-1);
	}
	rname = argv[1];
	if ((sscanf(argv[2], "%lf", &x) != 1) || (sscanf(argv[3], "%lf", &y) != 1) || (sscanf(argv[4], "%lf", &z) != 1)) {
		fprintf(stderr, usage);
		exit(-1);
	}

	/* allocate space for re[3], rref[3] */
	r = (double *) calloc(3, sizeof(double));
	rref = (double *) calloc(3, sizeof(double));
	/* and basis vectors */
	i0 = (double *) calloc(3, sizeof(double));
	i1 = (double *) calloc(3, sizeof(double));

	/* construct normalized reference direction */
	assign(rref, x, y, z);
	scale(rref, 1.0 / length(rref));

	/* compute polar coords */
	theta = acos(rref[2]);
	phi = atan2(rref[1], rref[0]);

	/* construct the unit vectors on the sky */
	assign(i0, - sin(phi), cos(phi), 0.0);
	assign(i1, cos(theta) * cos(phi), cos(theta) * sin(phi), - sin(theta));

	/* open lc pipe for input */
	sprintf(lccom, "lc -b -o %s", rname);
	ipf = popen(lccom, "r");
	if (!ipf) {
		fprintf(stderr, "r2n : failed to open lc-pipe for input\n");
		exit(-1);
	}

	/* open lc-pipe for output */
	opf = popen("lc -C -N '1 2 n'", "w");
	if (!opf) {
		fprintf(stderr, "r2n : failed to open lc-pipe for output\n");
		exit(-1);
	}

	while (1) {
		/* read vector r[3] */
		if (!fread(r, sizeof(double), 3, ipf)) {
			pclose(ipf);
			pclose(opf);
			exit(0);
		}
		scale(r, 1.0 / length(r));
		n0 = dot(r, i0);
		n1 = dot(r, i1);
		fprintf(opf, "%14.8lg %14.8lg\n", n0, n1);
	}
}

