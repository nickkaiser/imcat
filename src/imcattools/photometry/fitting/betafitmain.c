#define usage "\n\n\n\
NAME\n\
	betafit --- fit image to beta model\n\
\n\
SYNOPSIS\n\
	betafit [options....]\n\
		-u		# print this message\n\
		-f		# output the model as fits image\n\
		-N N		# number of micropixels per pixel (1.e4)\n\
		-v		# verbose mode\n\
		-r rc		# initial core radius (1.0)\n\
		-b beta		# assumed beta\n\
\n\
DESCRIPTION\n\
	\"betafit\" reads a fits image of counts n[][] from stdin and fits\n\
	this to a simple beta model for mean counts with central\n\
	value f0; core radius rc; index beta centered on x0, y0:\n\
\n\
		n_model = n0 * (1 + r^2 / rc^2)^-beta\n\
\n\
	where n0 = f0 * f and r^2 = (x - x0)^2 + (y - y0)^2\n\
\n\
	by minimising\n\
		sum_pixels N * f - sum_pixels (n[iy][ix] * log(f))\n\
\n\
	Where N is the number of micropixels per real pixel.\n\
\n\
	output is\n\
\n\
		N * f0, x0, y0, rc, beta\n\
\n\
	We start with x0, y0 given by mean of x,y for input counts\n\
	and core radius unity by default and with initial\n\
	f = sum n / (N * rc2)\n\
\n\
BUGS\n\
	Couldn't get it to work stably with beta as free parameter.\n\
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
#include "betafit.h"

#ifndef PI
#define PI M_PI
#endif


int		main(int argc, char *argv[])	
{
	int		arg = 1, N1, N2;
	int		pixtype, ix, iy;
	int		outputmodelimage, verbose;
	float		**n, **nmodel;
	float		f0, x0, y0, rc, rc2, beta, loglhood, N;
	float		nsum, xsum, ysum, xxsum, yysum;
	fitsheader	*fits;

	/* defaults */
	outputmodelimage = 0;
	N = 1.e6;
	verbose = 0;
	beta = 1.0;
	rc = 1.0;

	/* parse args */	
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		}
		switch (argv[arg++][1]) {
			case 'f':
				outputmodelimage = 1;
				break;
			case 'v':
				verbose = 1;
				break;
			case 'N':
				if (1 != sscanf(argv[arg++], "%f", &N)) {
					error_exit(usage);
				}
				break;
			case 'r':
				if (1 != sscanf(argv[arg++], "%f", &rc)) {
					error_exit(usage);
				}
				break;
			case 'b':
				if (1 != sscanf(argv[arg++], "%f", &beta)) {
					error_exit(usage);
				}
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}

	read2Dfloatimage(&n, &N1, &N2, &fits, stdin);

	/* compute mean count, x, y etc */
	nsum = xsum = ysum = xxsum = yysum = 0.0;
	for (iy = 0; iy < N2; iy++) {
		for (ix = 0; ix < N2; ix++) {
			nsum += n[iy][ix];
			xsum += ix * n[iy][ix];
			ysum += iy * n[iy][ix];
			xxsum += ix * ix * n[iy][ix];
			yysum += iy * iy * n[iy][ix];
		}
	}

	/* initialise parameters */
	x0 = xsum / nsum;
	y0 = ysum / nsum;
	rc2 = rc * rc;
	f0 = nsum / (N * rc2);	

	if (verbose) {
		fprintf(stderr, "# starting with:\n# f0=%13.8g\n# x0=%13.8g\n# y0=%13.8g\n# rc2=%13.8g\n", 
			f0, x0, y0, rc2);
	}

	/* do the fit */
	fitall(n, N1, N2, &f0, &x0, &y0, &rc2, beta, &loglhood, N);

	/* output results */
	if (outputmodelimage) {
		add_comment(argc, argv, fits);
		allocFloatArray(&nmodel, N1, N2);
		makebetamodel(nmodel, N1, N2, f0, x0, y0, rc2, beta, N);
		fits->extpixtype = FLOAT_PIXTYPE; 
		write2Dfloatimage(nmodel, fits);
		exit(0);
	}

	if (verbose) {
		fprintf(stdout, "#          n0            x0            y0            rc      loglhood\n");
	}
	fprintf(stdout, "%13.8g %13.8g %13.8g %13.8g %13.8g\n", N * f0, x0, y0, sqrt(rc2), -loglhood);
	exit(0);
}



