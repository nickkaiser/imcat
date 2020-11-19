#define usage "\n\n\n\
NAME\n\
	makevonkarmanS --- von Karman structure function\n\
SYNOPSIS\n\
	makevonkarmanS R r0 dz Nz [-N N] [-d dy] [-y ystar]\n\
\n\
DESCRIPTION\n\
	'makevonkarmanS' computes the function\n\
\n\
	S(z) = r0^(5/3) z^(-5/3) sum dy y (y^2 + (2 pi z / R)^2)^(-11/6)\n\
		(1 - J0(y) exp(-(y / ystar)^2))\n\
\n\
	for z = iz * dz and 0 <= iz < Nz.\n\
\n\
	Options are:\n\
		-N	Ny	# number of steps in y (10000)\n\
		-d	dy	# step size in y (0.01)\n\
		-y	ystar	# integration parameter (50)\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "utils/error.h"
#include "utils/args.h"
#include "imlib/fits.h"

int		main(int argc, char *argv[])	
{
	char		*flag, argstring[256];
	char		lcstr[256];
	int		N, Nz, iz, iy;
	double		r0, R, z, dz, y, dy, ystar, y0, yy, yy0, S, *b, gamma, fudgefac;

	/* defaults */
	N 	= 10000;
	dy	= 0.01;
	ystar	= 50.0;
	fudgefac = 6.88 / 0.98723873;

	/* parse args */
	argsinit(argc, argv, usage);
	if (nextargtype() == FLAG_ARG) {
		error_exit(usage);
	}
	R 	= getargd();
	r0 	= getargd();
	dz	= getargd();
	Nz	= getargi();
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'N':
				N = getargi();
				break;
			case 'd':
				dy = getargd();
				break;
			case 'y':
				ystar = getargd();
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}


	/* open the output stream */
	argsToString(argc, argv, argstring);
	sprintf(lcstr, "lc -C -x -a 'history: %s' -n z -n S -H 'R = %14.8lg' -H 'r0 = %14.8lg' < /dev/null", argstring, R, r0);
	system(lcstr);

	/* set up the bessel function */
	b = (double *) calloc(N, sizeof(double));
	for (iy = 0; iy < N; iy++) {
		y = iy * dy;
		b[iy] = (1 - j0(y) * exp(-(y / ystar) * (y / ystar)));
	}

	/* compute the integral */
	gamma = -11.0 / 6.0;
	fprintf(stdout, "  %14.8lg %14.8lg\n", 0.0, 0.0);
	for (iz = 1; iz < Nz; iz++) {
		z = iz * dz;
		y0 = 2 * M_PI * z / R;
		yy0 = y0 * y0;
		S = 0.0;
		for (iy = 0; iy < N; iy++) {
			y = dy * iy;
			yy = y * y;
			S += y * pow(yy + yy0, gamma) * b[iy];
		}
		S *= dy * pow(r0, - 5.0 / 3.0) * pow(z, 5.0 / 3.0) * fudgefac;
		fprintf(stdout, "  %14.8lg %14.8lg\n", z, S);
	}
	exit(0);
}

