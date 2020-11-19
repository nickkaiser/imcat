#define usage "\n\n\n\
NAME\n\
	stackplanes --- concatenate a set of fits files to multidimensional image\n\
\n\
SYNOPSIS\n\
	stackplanes fits1 fits2...\n\
	stackplanes nplanes\n\
\n\
DESCRIPTION\n\
	In the first invocation form (with 2 or more arguments)\n\
	\"stackplanes\" reads a set of N-dimensional images from\n\
	the files given as arguments and sends to stdout a single\n\
	(N+1)-dimensional fits images whose planes are the input\n\
	images.\n\
\n\
	stackplanes uses the 'iostream' mechanism, so you may provide,\n\
	in place of FITS file names, FITS file generating commands\n\
	(temrinated by a pipe symbol '|').\n\
\n\
	In the second invocation form, \"stackplanes\" reads a\n\
	single N-dimensional fits image, whose slowest dimension\n\
	must be an integral multiple of nplanes, from stdin and\n\
	generates a N+1 dimensional image whose planes are simply\n\
	concatenated in the input image.\n\
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
#include "../utils/iostream.h"

int		main(int argc, char *argv[])	
{
	int		arg, nplanes, ndim1, *n1, pixtype1, dim, nlines, line, good, singleimagemode;
	fitsheader	*fitsin, *fitsout;
	iostream	*ipstream;
	void		*f;

	/* check args */
	if (argc < 2) {
		error_exit(usage);
	}
	if (argc == 2) {
		if (!strcmp(argv[1], "-u")) error_exit(usage);
		if (1 != sscanf(argv[1], "%d", &nplanes)) error_exit(usage);
		fitsin = readfitsheader(stdin);
		if (fmod(fitsin->n[fitsin->ndim - 1], nplanes) > 0.0) {
			error_exit("stackplanes: slowest dimension must be integral multiple of nplanes\n");
		}
		fitsout = copyfitsheader(fitsin);
		fitsout->n[fitsout->ndim - 1] /= nplanes;
		fitsout->ndim++;
		fitsout->n[fitsout->ndim - 1] = nplanes;
		writefitsheader(fitsout);
		nlines = 1;
		for (dim = 1; dim < fitsin->ndim; dim++) {
			nlines *= fitsin->n[dim];
		}
		f = calloc(fitsin->n[0] * pixsize(fitsin->extpixtype), sizeof(char));
		for (line = 0; line < nlines; line++) {
			fread(f, sizeof(char), fitsin->n[0] * pixsize(fitsin->extpixtype), fitsin->ipstream);
			fwrite(f, sizeof(char), fitsin->n[0] * pixsize(fitsin->extpixtype), fitsout->opstream);
		}
		writefitstail(fitsout);
		free(f);
		exit(0);			
	}

	/* process files */
	for (arg = 1; arg < argc; arg++) {
		ipstream = openiostream(argv[arg], "r");
		if (!ipstream) {
			error_exit("stackplane: failed to open input file\n");
		}
		fitsin = readfitsheader(ipstream->f);
		if (arg == 1) {
			pixtype1 = fitsin->extpixtype;
			ndim1 = fitsin->ndim;
			n1 = (int *) calloc(ndim1, sizeof(int));
			nlines = 1;
			for (dim = 0; dim < ndim1; dim++) {
				n1[dim] = fitsin->n[dim];
			}
			nlines = 1;
			for (dim = 1; dim < ndim1; dim++) {
				nlines *= fitsin->n[dim];
			}
			fitsout = copyfitsheader(fitsin);
			add_comment(argc, argv, fitsout);
			nplanes = argc - 1;
			fitsout->n[fitsout->ndim] = nplanes;
			fitsout->ndim++;
			writefitsheader(fitsout);
			f = calloc(fitsin->n[0] * pixsize(pixtype1), sizeof(char));
		} else {
			good = 1;
			if (fitsin->extpixtype != pixtype1) good = 0;
			if (fitsin->ndim != ndim1) good = 0;
			for (dim = 0; dim < ndim1; dim++) {
				if (fitsin->n[dim] != n1[dim]) good = 0;
			}
			if (!good) {
				error_exit("stackplane: images must have identical dimensions and pixtype\n");
			}
		}
		for (line = 0; line < nlines; line++) {
			fread(f, sizeof(char), fitsin->n[0] * pixsize(pixtype1), fitsin->ipstream);
			fwrite(f, sizeof(char), fitsin->n[0] * pixsize(pixtype1), fitsout->opstream);
		}
		closeiostream(ipstream);
	}
	writefitstail(fitsout);
	free(f);

	exit(0);
}
