#define usage "\n\n\n\
NAME\n\
	getplanes --- extracts a set of planes from multi-dimensional fits image\n\
\n\
SYNOPSIS\n\
	getplanes [-f srcfits] i1 i2...\n\
\n\
DESCRIPTION\n\
	\"getplanes\" reads from stdin a N-dimensional fits image\n\
		f_in[Nm]....[N2][N1]\n\
	and writes to stdout a N dimensional image with planes\n\
		f_out[i1]....[N2][N1]\n\
		f_out[i2]....[N2][N1]\n\
	etc.\n\
\n\
	With -f option we read from named file (and we use fseek() to\n\
	move around in the file for efficiency).\n\
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

int		main(int argc, char *argv[])	
{
	int		npin, npout, pixsizein, pixsizeout, *i, arg, planesize, j, readfromfile, ind;
	long		offset, baseoffset;
	fitsheader	*fitsin, *fitsout;
	char		*srcfitsname;
	void		*f;
	FILE		*ipf;

	/* defaults */
	readfromfile = 0;
	ipf = stdin;

	/* parse args */
	if (argc < 2) {
		error_exit(usage);
	}
	arg = 1;
	if (argv[arg][0] == '-') {
		if (argv[arg][1] == 'f') {
			if (argc > 3) {
				arg++;
				srcfitsname = argv[arg];
				readfromfile = 1;
				ipf = fopen(srcfitsname, "r");
				if (!ipf) {
					error_exit("getplanes: failed to open srcfits\n");
				}
				arg++;
			} else {
				error_exit(usage);
			}
		} else {
			error_exit(usage);
		}
	}

	/* read image header */
	fitsin = readfitsheader(ipf);
	npin = fitsin->n[fitsin->ndim - 1];

	/* create the array of indices */
	npout = argc - arg;
	i = (int *) calloc(npout, sizeof(int));
	for (; arg < argc; arg++) {
		ind = arg - argc + npout;
		if (1 != sscanf(argv[arg], "%d", i + ind)) {
			error_exit("getplanes: error parsing arg\n");
		}
		if ((i[ind] < 0) || (i[ind] >= npin)) {
			error_exit("getplanes: index out of range\n");
		}		
	}
	
	/* construct and write output image header */
	fitsout = copyfitsheader(fitsin);
	fitsout->n[fitsout->ndim - 1] = npout;
	add_comment(argc, argv, fitsout);
	writefitsheader(fitsout);
	
	/* compute planesize */
	planesize = pixsize(fitsin->extpixtype);
	for (j = 0; j < fitsin->ndim - 1; j++) {
		planesize *= fitsin->n[j];
	}	

	if (readfromfile) {
		/* allocate enough space for a single plane */
		f = (void *) calloc(planesize, sizeof(char));
		baseoffset = ftell(ipf);
		for (j = 0; j < npout; j++) {
			offset = baseoffset + i[j] * planesize;
			fseek(ipf, offset, SEEK_SET);
			fread(f, planesize, sizeof(char), ipf);
			fwrite(f, planesize, sizeof(char), stdout);
		}
	} else {
		/* we are reading from stdin, so we read the whole FITS file */
		/* allocate and read the planes */
		f = (void *) calloc(npin, sizeof(void *));
		for (j = 0; j < npin; j++) {
			((void **) f)[j] = (void *) calloc(planesize, sizeof(char));
			fread(((void **) f)[j], planesize, sizeof(char), ipf);
		}
		/* write the planes */
		for (j = 0; j < npout; j++) {
			fwrite(((void **) f)[i[j]], planesize, sizeof(char), stdout);
		}
	}

	writefitstail(fitsout);
	fclose(ipf);
	exit(0);
}
