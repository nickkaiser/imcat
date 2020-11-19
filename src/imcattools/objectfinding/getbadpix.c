#define	usage "\n\n\n\
NAME\n\
	getbadpix --- reject discrepant pixels\n\
\n\
SYNOPSIS\n\
	getbadpix nu referenceimage [option...] \n\
		-s	sigma	# rms noise value\n\
		-l		# output pixel values as well as positions\n\
		-f	d rg	# allow for difference in psf shape\n\
\n\
DESCRIPTION\n\
	By default \"getbadpix\" reads a fits image from stdin, extracts\n\
	the value of the image header item with keyword 'SIGMA', and\n\
	sends to stdout an 'lc' format catalogue containing the a list\n\
	of pixels for which the pixel value differs from the reference\n\
	image pixel value by more than nu * SIGMA.\n\
\n\
	By default, the catalogue contains sinply the pixel position 'x[2]',\n\
	but with the -l option it will also contain the input and reference\n\
	image pixel values.\n\
\n\
	If the '-s' option is given the rms noise is read from the following\n\
	command line argument rather than from the header.\n\
\n\
	The -f option is provided to allow for the fact the for bright\n\
	objects such as stars the difference between the source and\n\
	reference images may greatly exceed the statistical nu * sigma\n\
	limit, so instead we reject pixels if |f - fref| exceeds the\n\
	greater of:\n\
\n\
		nu  sigma\n\
		2 d fref\n\
		d fref rg^2 grad(fref)^2 / (fref^2 + sigma^2)\n\
\n\
	the last two expressions being the difference in the flux of an\n\
	approximately gaussian star of radius r in the centre and on the\n\
	edge where d is the assumed fractional change in the gaussian\n\
	scale length: d = delta rg / rg = d ln(rg).\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n\n\n"		

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../imlib/fits.h"
#include "../../utils/error.h"
#include "../../catlib/cat.h"

#define MAGIC FLOAT_MAGIC


main(int argc, char *argv[])	
{
	float 		**f, **fref;
	int		arg, N1, N2, N1ref, N2ref, ix, iy; 
	fitsheader	*fitsin, *fitsref;
	char		*refimagename;
	FILE		*refimagef, *lcpipe;
	cathead		*thecathead;
	object		*theobject;
	item		*xitem, *fitem, *frefitem;
	double		x[2], ff, ffref, nu, sigma, dlnrg, rg, df;
	double		dfmaxstat, dfmaxedge, dfmaxcentre;
	int		needheadersigma, outputpixvals, allowfractionalerror;
	double		dfx1, dfx2, dfy1, dfy2, gradfsq;
	
	/* defaults */
	needheadersigma = 1;
	outputpixvals = 0;
	allowfractionalerror = 0;

	/* parse args */
	if (argc < 3) {
		error_exit(usage);
	}
	arg = 1;
	if (1 != sscanf(argv[arg++], "%lf", &nu)) {
		error_exit(usage);
	}
	refimagename = argv[arg++];
	while (arg < argc) {
		if (*argv[arg] != '-') {
			error_exit(usage);
		}
		switch (*(argv[arg++]+1)) {
			case 's':
				if (1 != sscanf(argv[arg++], "%lf", &sigma))
					error_exit(usage);
				needheadersigma = 0;
				break;
			case 'l':
				outputpixvals = 1;
				break;
			case 'f':
				if (1 != sscanf(argv[arg++], "%lf", &dlnrg))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%lf", &rg))
					error_exit(usage);
				allowfractionalerror = 1;
				break;
			default:
				error_exit(usage);
				break;
		}
	}

	/* open reference image and read header */
	refimagef = fopen(refimagename, "r");
	if (!refimagef) {
		error_exit("getbadpix: can't open reference image fits file\n");
	}
	read2Dfloatimage(&fref, &N1ref, &N2ref, &fitsref, refimagef);

	/* read source image header from stdin and get sigma if nbecessary */
	read2Dfloatimage(&f, &N1, &N2, &fitsin, stdin);
	if (needheadersigma) {
		 sigma = getnumericvalue(getcommentbyname("SIGMA", fitsin));
	}
	
	/* check the images have the same size */
	if ((N1 != N1ref) || (N2 != N2ref)) {
		error_exit("getbadpix: input and reference images must have same dimensions\n");
	}

	/* create the cathead with 'lc */
	if (outputpixvals) {
		lcpipe = popen("lc -C -x -N '1 2 x' -n f -n fref < /dev/null", "r");
	} else {
		lcpipe = popen("lc -C -x -N '1 2 x' < /dev/null", "r");
	}
	setcatipf(lcpipe);
	thecathead = readcathead();
	pclose(lcpipe);

	/* now add the history */
	addargscomment(argc, argv, thecathead);

	/* and write the cathead */
	setcatopfiletype(BINARY_FILE_TYPE);
	writecathead(thecathead);

	/* now create the object */
	theobject = newobject(thecathead);
	connectcatheadtoobject(theobject);
	setaddress(theobject, getobjectitemindex("x", theobject), x);
	if (outputpixvals) {
		setaddress(theobject, getobjectitemindex("f", theobject), &ff);
		setaddress(theobject, getobjectitemindex("fref", theobject), &ffref);
	}

	/* compute the threshold value for statistical noise */
	dfmaxstat = nu * sigma;

	/* and now find hotpixels */
	for (iy = 0; iy < N2; iy++) {
		for (ix = 0; ix < N1; ix++) {
			if (f[iy][ix] == MAGIC || fref[iy][ix] == MAGIC) {
				continue;
			}
			ff = (double) f[iy][ix];
			ffref = (double) fref[iy][ix];
			df = fabs(ff - ffref);
			if (df > dfmaxstat) {
				if (allowfractionalerror) {
					dfmaxcentre = 2 * dlnrg * fabs(ffref);
					if (df < dfmaxcentre) {
						continue;
					}
					dfx1 = (ix == 0 ? 0 : (fref[iy][ix - 1] == MAGIC ? 0 : fref[iy][ix] - fref[iy][ix - 1]));
					dfx2 = (ix == (N1 - 1) ? 0 : (fref[iy][ix + 1] == MAGIC ? 0 : fref[iy][ix + 1] - fref[iy][ix]));
					dfy1 = (iy == 0 ? 0 : (fref[iy - 1][ix] == MAGIC ? 0 : fref[iy][ix] - fref[iy - 1][ix]));
					dfy2 = (iy == (N2 - 1) ? 0 : (fref[iy + 1][ix] == MAGIC ? 0 : fref[iy + 1][ix] - fref[iy][ix]));
					gradfsq = 0.5 * (dfx1 * dfx1 + dfx2 * dfx2 + dfy1 * dfy1 + dfy2 * dfy2);
					dfmaxedge = dlnrg * ffref * rg * rg * gradfsq / (ffref * ffref + sigma * sigma);
					if (df < dfmaxedge) {
						continue;
					}
				}
				x[0] = (double) ix + 0.5;
				x[1] = (double) iy + 0.5;
				ff = df;
				ffref = dfmaxedge;
				writeobject(theobject);
			}
		}
	}
	exit(0);
}










