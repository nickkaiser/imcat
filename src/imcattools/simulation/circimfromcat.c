/*
 * generate a circularly symmetric image from a 1-D catalogue
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>

#include "imlib/fits.h"
#include "utils/error.h"
#include "utils/arrays.h"


char *usage = "\n\
NAME\n\
	circimfromcat --- generate a circularly symmetric image\n\
\n\
SYNOPSIS\n\
	circimfromcat [options]\n\
\n\
DESCRIPTION\n\
	'circimfromcat' reads an lc-format catalogue from standard\n\
	input containing at least a single scalar (default name F)\n\
	which is interpreted as a set uf uniformly spaced samples\n\
	of a function F(r), and generates a square circularly\n\
	symmetric 2-D image f(x,y) = F(r = sqrt(x^2 + y^2)).\n\
	Options are:\n\
		-d dr		# spacing of input samples in pixels (1.0)\n\
		-n N		# size of output image (512)\n\
		-N nrmax	# input array size (10000)\n\
		-F Fname	# name for the input value\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n";


main(int argc, char *argv[])
{
	int		arg = 1, ir, nr, ix, iy, nrmax;
	int		N, pixtype;
	char 		lcstring[128], *defname = "F", *thename;
	double		dr, *F, x, y, r, *ipbuff[1];
	float		**f;
	FILE		*lcpipe;
	fitsheader      *fits;

	/* defaults */
	dr = 1.0;
	N = 512;
	nrmax = 10000;
	thename = defname;

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		}
		switch (argv[arg++][1]) {
			case 'd':
				if (1 != sscanf(argv[arg++], "%lf", &dr)) {
					error_exit(usage);
				}
				break;
			case 'n':
				if (1 != sscanf(argv[arg++], "%d", &N)) {
					error_exit(usage);
				}
				break;
			case 'N':
				if (1 != sscanf(argv[arg++], "%d", &nrmax)) {
					error_exit(usage);
				}
				break;
			case 'F':
				thename = argv[arg++];
				break;
			case 'u':
			default:
				error_exit(usage);
				break;
		}
	}

	/* allocate space for the cl array */
	F = (double *) calloc(nrmax, sizeof(double));
	
	/* read the cl data */
	ir = 0;
	sprintf(lcstring, "lc -o -b %s", thename);
	lcpipe = popen(lcstring, "r");
	if (!lcpipe) {
		error_exit("circimfromcat: failed to open lc-pipe for input\n");
	}
	while(fread(F + ir++, sizeof(double), 1, lcpipe)) {
		if (ir == nrmax) {
			error_exit("circimfromcat: use a bigger input array\n");
		}
	}
	pclose(lcpipe);
	nr = ir;

	/* generate the 2-D power spectrum image */
	allocFloatArray(&f, N, N);
	for (iy = 0; iy < N; iy++) {
		y = iy - N / 2;
		for (ix = 0; ix < N; ix++) {
			x = ix - N / 2;
			r = sqrt(x * x + y * y);
			ir = floor(0.5 + r / dr);
			if (ir >= 0 && ir < nr) {
				f[iy][ix] = (float) F[ir];
			}
		}
	} 

	/* and output */
        fits = new2Dfitsheader(N, N, FLOAT_PIXTYPE);
        add_comment(argc, argv, fits);
        write2Dfloatimage(f, fits);
	exit(0);
}


