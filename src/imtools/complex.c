#define usage "\n\n\n\
NAME\n\
	complex --- perform operations on complex images\n\
SYNOPSIS\n\
	complex op [files...]\n\
\n\
DESCRIPTION\n\
	Complex works with 2-D complex images stored in f[2][n2][n1] FITS format.\n\
	Examples\n\
\n\
		complex multiply f1 f2\t\n\
		complex conjugate f\t\n\
		complex ab2rphi f\t\n\
		complex rphi2ab f\t\n\
\n\
	where f, f1, f2 are filenames or '-' for standard input.\n\
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
#include "../utils/arrays.h"

#define		MULTIPLY 	0
#define		CONJUGATE	1
#define		AB2RPHI		2
#define		RPHI2AB		3

fitsheader *openfits(char *fn);

int		main(int argc, char *argv[])	
{
	int		x, y, mode, N1, N2;
	fitsheader	*fits;
	float		**ar, **ai, **br, **bi, **cr, **ci;
	char		*op;
	fitsheader	*fitsa, *fitsb, *fitsc;
	
	argsinit(argc, argv, usage);
	if (nextargtype() == FLAG_ARG || nextargtype() == NO_ARG) {
		error_exit(usage);
	}
	op = getargs();
	fitsa = openfits(getargs());
	N1 = fitsa->n[0];
	N2 = fitsa->n[1];
	allocFloatArray(&ar, N1, N2);
	allocFloatArray(&ai, N1, N2);
	readfitsplane((void **) ar, fitsa);
	readfitsplane((void **) ai, fitsa);
	allocFloatArray(&cr, N1, N2);
	allocFloatArray(&ci, N1, N2);

	switch (op[0]) {
		case 'm':
			mode = MULTIPLY;
			fitsb = openfits(getargs());
			if (fitsb->n[0] != N1 || fitsb->n[1] != N2) {
				error_exit("complex : input images must have same size\n");
			}
			allocFloatArray(&br, N1, N2);
			allocFloatArray(&bi, N1, N2);
			readfitsplane((void **) br, fitsb);
			readfitsplane((void **) bi, fitsb);
			break;
		case 'c':
			mode = CONJUGATE;
			break;
		case 'a':
			mode = AB2RPHI;
			break;
		case 'r':
			mode = RPHI2AB;
			break;
		default:
			error_exit(usage);			
	}

	for (y = 0; y < N2; y++) {
		for (x = 0; x < N1; x++) {
			switch(mode) {
				case MULTIPLY:
					cr[y][x] = ar[y][x] * br[y][x] - ai[y][x] * bi[y][x];
					ci[y][x] = ar[y][x] * bi[y][x] + ai[y][x] * br[y][x];
					break;
				case CONJUGATE:
					cr[y][x] = ar[y][x];
					ci[y][x] = -ai[y][x];
					break;
				case AB2RPHI:
					cr[y][x] = sqrt(ar[y][x] * ar[y][x] + ai[y][x] * ai[y][x]);
					ci[y][x] = atan2(ai[y][x], ar[y][x]);
					break;
				case RPHI2AB:
					cr[y][x] = ar[y][x] * cos(ai[y][x]);
					ci[y][x] = ar[y][x] * sin(ai[y][x]);
					break;
				default:
					error_exit("complex: illegal mode\n");
			}
		}
	}

	fitsc = fitsa;
	add_comment(argc, argv, fitsc);
	writefitsheader(fitsc);
	writefitsplane((void **) cr, fitsc);
	writefitsplane((void **) ci, fitsc);
	writefitstail(fitsc);
	exit(0);
}


fitsheader *openfits(char *fn)
{
	FILE		*f;
	fitsheader 	*fits;

	if (strcmp(fn, "-")) {
		f = fopen(fn, "r");
		if (!f) {
			error_exit("complex: faile to open input image\n");
		}
	} else {
		f = stdin;
	}
	fits = readfitsheader(f);
	if (fits->ndim != 3 || fits->n[2] != 2) {
		fprintf(stderr, "complex: %s doesn;t look like a complex image\n", fn);
		exit(1);
	}
	return(fits);
}

