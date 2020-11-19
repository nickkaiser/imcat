/*
 * makesheariage.c
 */

#define usage "\n\
NAME\n\
	makeshearimage --- generate a shear image from catalogue\n\
\n\
SYNOPSIS\n\
	makeshearimage [options...]\n\
		-N N1 N2	# size of image (128,128)\n\
		-x xname	# name for spatial coordinate ('x')\n\
		-X dx		# pixel size in x-coords (16.0)\n\
		-o x0 y0	# origin in x-space (0.0, 0.0)\n\
		-r x1 x2 y1 y2	# box for calculating nbar\n\
		-e ename	# name for polarisation ('e')\n\
		-p psh		# polarisability (1.0)\n\
		-n nbar		# density of points in x-space\n\
\n\
DESCRIPTION\n\
	'makeshearimage' reads a catalogue which must contain a spatial\n\
	2-vector coordinate ('x' by default) and a 2-vector\n\
	polarisation ('e' by default), and creates a N1 by (2 N2)\n\
	fits image which is a simple binned average of the input\n\
	shear values. The shear is defined as\n\
		gamma_i = e_i / psh,\n\
	and gamma_0 and gamma_1 are stored in the first and last N2 lines\n\
	of the output image respectively.\n\
	By default the density of points is estimated as the number\n\
	of objects divided by the area of the image (in x-coord units)\n\
	but you can supply an alternative x-rectange with the -r\n\
	option.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n"

#include <stdio.h>
#include <math.h>
#include "../../imlib/fits.h"
#include "../../utils/arrays.h"

main (int argc, char *argv[])
{
	int		arg = 1, N1, N2, ix, iy, pixtype, count;
	int		neednbar, needrange;
	float		**gamma0, **gamma1;
	double		*x, *e, psh, dx, x0[2], nbar, x1, x2, y1, y2;
	double		*lcrec;
	char		lcstring[128];
	char		*xname, defxname[2] = "x", *ename, defename[2] = "e";
	FILE		*lcpipe;
	fitsheader	*fits;

	/* defaults */
	N1 = 128;
	N2 = 128;
	psh = 1;
	dx = 16.0;
	x0[0] = x0[1] = 0.0;
	neednbar = 1;
	xname = defxname;
	ename = defename;
	needrange = 1;
	

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			fprintf(stderr, "%s", usage);
			exit(-1);
		}
		switch (argv[arg++][1]) {
			case 'N':
				sscanf(argv[arg++], "%d", &N1);
				sscanf(argv[arg++], "%d", &N2);
				break;
			case 'x':
				xname = argv[arg++];
				break;
			case 'X':
				sscanf(argv[arg++], "%lf", &dx);
				break;
			case 'o':
				sscanf(argv[arg++], "%lf", x0);
				sscanf(argv[arg++], "%lf", x0 + 1);
				break;
			case 'r':
				sscanf(argv[arg++], "%lf", &x1);
				sscanf(argv[arg++], "%lf", &x2);
				sscanf(argv[arg++], "%lf", &y1);
				sscanf(argv[arg++], "%lf", &y2);
				if ((x2 - x1) <= 0.0 || (y2 - y1) <= 0.0)
					error_exit("makeshearimage: args for -r option!\n");
				needrange = 0;
				break;
			case 'e':
				ename = argv[arg++];
				break;
			case 'p':
				sscanf(argv[arg++], "%lf", &psh);
				break;
			case 'n':
				sscanf(argv[arg++], "%lf", &nbar);
				neednbar = 0;
				break;
			default:
				fprintf(stderr, "%s", usage);
				exit(-1);
				break;
		}
	}

	/* create output image */
	allocFloatArray(&gamma0, N1, 2 * N2);
	gamma1 = gamma0 + N2;

	/* get rectangle for nbar */
	if (needrange) {
		x1 = y1 = 0.0;
		x2 = N1 * dx;
		y2 = N2 * dx;
	}

	/* create input record */
	lcrec = (double *) calloc(4, sizeof(double));
	x = lcrec;
	e = lcrec + 2;

	/* open the lc-pipe */
	sprintf(lcstring, "lc -b -o %s %s", xname, ename);
	lcpipe = popen(lcstring, "r");
	if (!lcpipe)
		error_exit("makeshearimage: unable to open lc-pipe for input\n");

	/* accumulate gamma arrays */
	count = 0;
	while (4 == fread(lcrec, sizeof(double), 4, lcpipe)) {
		ix = floor((x[0] - x0[0]) / dx);
		iy = floor((x[1] - x0[1]) / dx);
		if (ix < 0 || ix >= N1 || iy < 0 || iy >= N2)
			continue;
		gamma0[iy][ix] += e[0] / psh;
		gamma1[iy][ix] += e[1] / psh;
		if (x[0] >= x1 && x[0] < x2 && x[1] >= y1 && x[1] < y2)
			count++;
	}
	if (!count && neednbar)
		error_exit("makeshearimage: no objects for determining nbar\n");

	/* figure nbar */
	if (neednbar) {
		nbar = count / ((x2 - x1) * (y2 - y1));
	}

	/* normalise the gamma arrays */
	nbar = nbar * dx * dx;
	for (iy = 0; iy < N2; iy++) {
		for (ix = 0; ix < N1; ix++) {
			gamma0[iy][ix] /= nbar;
			gamma1[iy][ix] /= nbar;
		}
	}
				
	fits = new2Dfitsheader(N1, N2, FLOAT_PIXTYPE);
	add_comment(argc, argv, fits);
	fits->ndim = 3;
	fits->n[2] = 2;
	writefitsheader(fits);
	for (iy = 0; iy < 2 * N2; iy++) {
		writefitsline(gamma0[iy], fits);
	}
	writefitstail(fits);
	exit(0);
}





