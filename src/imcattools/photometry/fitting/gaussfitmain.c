#define usage "\n\n\n\
NAME\n\
	gaussfit --- fit image to gaussian ellipsoid object model\n\
\n\
SYNOPSIS\n\
	gaussfit [options....]\n\
		-u		# print this message\n\
		-f		# output the model as fits image\n\
		-n	n	# fit n gaussians\n\
\n\
DESCRIPTION\n\
	\"gaussfit\" reads a fits image f from stdin and fits\n\
	to this a simple gaussian ellipsoid model with central\n\
	value f0; semi-axes a,b; position angle phi; at x0, y0:\n\
\n\
		f(x,y) = f0 * exp(-((X / a)^2 + (Y / b)^2) / 2)\n\
\n\
	where\n\
\n\
		X = (x - x0) cos(phi) + (y - y0) sin(phi)\n\
		Y = (y - y0) cos(phi) - (x - x0) sin(phi)\n\
\n\
	Output is:\n\
		f0  x0  y0  a  b  phi\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@ifa.hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../../imlib/fits.h"
#include "../../../utils/error.h"
#include "gaussfitn.h"

#ifndef PI
#define PI M_PI
#endif

#define NMAX 10

int		main(int argc, char *argv[])	
{
	int		arg = 1, N1, N2, i, n;
	int		comc, pixtype;
	int		invertq, outputabphi, outputmodelimage;
	float		**f, **fmodel;
	float		f0[NMAX], x[NMAX], y[NMAX], a[NMAX], b[NMAX], phi[NMAX];
	fitsheader	*fits;

	/* defaults */
	outputmodelimage = 0;
	n = 1;

	/* parse args */	
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		}
		switch (argv[arg++][1]) {
			case 'f':
				outputmodelimage = 1;
				break;
			case 'n':
				if (arg == argc) {
					error_exit(usage);
				}
				if (1 != sscanf(argv[arg++], "%d", &n)) {
					error_exit(usage);
				}
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}

	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);
	allocFloatArray(&fmodel, N1, N2);

	/* initialise parameters */
	for (i = 0; i < n; i++) {	
		x[i] = N1 * drand48();
		y[i] = N2 * drand48();
		f0[i] = 1000.0;
		a[i] = b[i] = 2.0;
		phi[i] = 3.14159 * drand48();
	}

	/* do the fit */
	gaussfitn(f, N1, N2, f0, x, y, a, b, phi, n, fmodel);

	/* output results */
	if (outputmodelimage) {
		add_comment(argc, argv, fits);
		write2Dfloatimage(fmodel, fits);
		exit(0);
	}

	for (i = 0; i < n; i++) {
		fprintf(stdout, "%13g %13g %13g %13g %13g %13g\n", 
		f0[i], x[i], y[i], a[i], b[i], phi[i]);
	}
	exit(0);
}



