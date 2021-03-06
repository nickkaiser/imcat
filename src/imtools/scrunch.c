#define usage "\n\n\n\
NAME\n\
	scrunch - demagnify image by factor 2\n\
\n\
SYNOPSIS\n\
	scrunch [-m | -c]\n\
		-m	# take median (mean is default)\n\
		-c	# conservative mean\n\
\n\
DESCRIPTION\n\
	Scrunches a fits image down by a factor 2.\n\
	Works in stream processing mode.\n\
	With -c option output is MAGIC if any of source pixels are MAGIC.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../imlib/fits.h"
#include "../utils/error.h"
#include "../utils/args.h"
#include "../imlib/scrunch_stuff.h"

int		main(int argc, char *argv[])	
{
	int		i, nplanes, N1in, N1out, N2in, N2out, arg = 1, mode;
	fitsheader	*fitsin, *fitsout;
	char		*flag;

	/* defaults */
	mode = MEAN;

	/* parse args */
	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'm':
				mode = MEDIAN;
				break;
			case 'c':
				mode = CMEAN;
				break;
			default:
				error_exit(usage);
		}
	}	
	
	fitsin = readfitsheader(stdin);
	N1in = fitsin->n[0];
	N2in = fitsin->n[1];
	nplanes = 1;
	for (i = 2; i < fitsin->ndim; i++) {
		nplanes *= fitsin->n[i];
	}
	if ((2 * (N1in / 2) != N1in) || (2 * (N2in / 2) != N2in)) 
		error_exit("scrunch: image dimensions must be even");
	fitsout = copyfitsheader(fitsin);
	fitsout->n[0] = N1in / 2;
	fitsout->n[1] = N2in / 2;
	add_comment(argc, argv, fitsout);
	writefitsheader(fitsout);
	for (i = 0; i < nplanes; i++) {
		scrunch_stream(fitsin, fitsout, mode);
	}
	writefitstail(fitsout);
	exit(0);
}

