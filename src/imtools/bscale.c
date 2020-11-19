#define usage "\n\n\n\
NAME\n\
	bscale - set or check BSCALE, BZERO header values in FITS file\n\
\n\
SYNOPSIS\n\
	bscale [fmin fmax | -g | -u ]\n\
\n\
DESCRIPTION\n\
	\"bscale fmin fmax\" reads an image (which must be an integer\n\
	format BITPIX = 8, 16 or 32) from stdin and modifies\n\
	or generates BZERO, BSCALE values such that the max and\n\
	min representable values are fmin, fmax.\n\
	This has the effect of rescaling the image values.\n\
\n\
	Use \"bscale -g\" to get the current fmin, fmax values.\n\
\n\
	Use -u option to print this message.\n\
\n\
	Beware that the minimum representable values for\n\
	BITPIX 16, 32 and the maximum representable values\n\
	for BITPIX 8 are 'magic' values used by imcat to flag bad\n\
	or missing data.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


#include "../imlib/fits.h"
#include "../utils/error.h"
#include "../utils/args.h"

#define SETRANGEMODE	0
#define GETRANGEMODE 	1


int		main(int argc, char *argv[])	
{
	fitsheader	*fits;
	int		opmode, ifmin, ifmax, dim, nlines, linesize;
	char		*flag;
	double		fmin, fmax;
	void		*f;

	/* defaults */
	opmode = SETRANGEMODE;

	/* parse args */
	argsinit(argc, argv, usage);
	while (nextargtype() == FLAG_ARG) {
		flag = getflag();
		switch(flag[0]) {
			case 'g':
				opmode = GETRANGEMODE;
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}
	if (opmode == SETRANGEMODE) {
		fmin = getargd();
		fmax = getargd();
	}

	/* read the fits header */
	fits = readfitsheader(stdin);

	/* set bscaling if necessary */
	if (!(fits->bscaling)) {
		fits->bscaling = 1;
		fits->bscale = 1.0;
		fits->bzero = 0.0;
	}

	/* get the unscaled range */
	switch (fits->extpixtype) {
		case UCHAR_PIXTYPE:
			ifmin = 0;
			ifmax = UCHAR_MAX;
			break;
		case SHORT_PIXTYPE:	
			ifmin = SHRT_MIN;
			ifmax = SHRT_MAX;
			break;
		case INT_PIXTYPE:
			ifmin = INT_MIN;
			ifmax = INT_MAX;
			break;
		case FLOAT_PIXTYPE:
		case DBL_PIXTYPE:
			error_exit("bscale: error: non-integer image\n");
	}

	if (opmode == GETRANGEMODE) {
		fmin = fits->bzero + fits->bscale * ifmin;
		fmax = fits->bzero + fits->bscale * ifmax;
		fprintf(stdout, "%14.8g %14.8g\n", fmin, fmax);
		exit(0);
	} else {
		fits->bscale = (fmax - fmin) / (ifmax - ifmin);
		fits->bzero = fmax - fits->bscale * ifmax;
	}

	/* write the fitsheader */
	add_comment(argc, argv, fits);
	writefitsheader(fits);

	/* and process the data */
	linesize = fits->n[0] * pixsize(fits->extpixtype);
	nlines = 1;
	for (dim = 1; dim < fits->ndim; dim++) {
		nlines *= fits->n[dim];
	}
	f = calloc(linesize, sizeof(char));
	while (nlines--) {
		fread(f, sizeof(char), linesize, stdin);
		fwrite(f, sizeof(char), linesize, stdout);
	}
	writefitstail(fits);
	exit(0);
}


