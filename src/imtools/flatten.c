#define usage "\n\n\n\
NAME\n\
	flatten --- convert 3-D fits to 2_D\n\
\n\
SYNOPSIS\n\
	flatten [options...] \n\
\n\
DESCRIPTION\n\
	\"flatten\" reads a 3-dimensional fits file from standard\n\
	input and writes a 2-dimensional version to stdout.\n\
	By default the output image is simply the concatenation\n\
	of the input image planes.\n\
\n\
	Options are:\n\
		-c ncols	# write images in rows ncols wide (1)\n\
		-g		# flatten a 60 plane image a la GPC1\n\
		-G		# as above, but rotate RHS images by 180 deg\n\
		-b dx dy	# pad with a boundary of dx and dy\n\
\n\
BUGS\n\
	With -b option, we only set the boundary to MAGIC for short (signed 16 bit)\n\
	or float pixel types. Otherwise, boundary is set to zero internally, but will\n\
	be non-zero on output if BZERO is non-zero.\n\
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
	int		arg, i, j, pixsiz, nlines, nbytes, ncols, dx, dy, N1, N2, nrows, nplanes, nblanks;
	int		row, col, plane, offset, gpc1format = 0, gpc1rotaterhs = 0;
	fitsheader	*fits;
	void		*f;
	short		*fshort;
	float		*ffloat;

	/* defaults */
	ncols = 1;
	dx = dy = 0;

	/* parse args */
	arg = 1;
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		}
		switch (argv[arg++][1]) {
			case 'c':
				sscanf(argv[arg++], "%d", &ncols);
				break;
			case 'g':
				gpc1format = 1;
				break;
			case 'G':
				gpc1format = 1;
				gpc1rotaterhs = 1;
				break;
			case 'b':
				sscanf(argv[arg++], "%d", &dx);
				sscanf(argv[arg++], "%d", &dy);
				break;
			default:
				error_exit(usage);
		}
	}

	fits = readfitsheader(stdin);
	N1 = fits->n[0];
	N2 = fits->n[1];
	add_comment(argc, argv, fits);
	/* size of a pixel in bytes */
	pixsiz = pixsize(fits->extpixtype);
	/* number of bytes per line */
	nbytes = fits->n[0] * pixsiz;
	if (gpc1format) {
		nplanes = 1;
		for (i = 2; i < fits->ndim; i++) {
			nplanes *= fits->n[i];
		}
		if (nplanes != 60) {
			error_exit("flatten : error : expecting 60 plane image in with -g option");
		}
		ncols = nrows = 8;
	} else {
		nlines =  1;
		for (i = 1; i < fits->ndim; i++) {
			nlines *= fits->n[i];
		}
		nplanes = nlines / fits->n[1];
		nrows = (int) ceil((double) nplanes / ncols);
		nblanks = nrows * ncols - nplanes;
	}

	fits->ndim = 2;
	fits->n[0] = ncols * (N1 + 2 * dx);
	fits->n[1] = nrows * (N2 + 2 * dy);
	writefitsheader(fits);
	for (row = 0; row < nrows; row++) {
		f = calloc(ncols * (N2 + 2 * dy) * (nbytes + 2 * dx * pixsiz), sizeof(char));
		if (dx || dy) {
			/* set pixels to magic */
			switch (fits->extpixtype) {
				case SHORT_PIXTYPE:
					fshort = (short *) f;
					for (i = 0; i < ncols * (N1 + 2 * dx) * (N2 + 2 * dy); i++) {
						fshort[i] = SHORT_MAGIC;
					}
					if (fits->opbyteorder == NON_NATIVE_BYTE_ORDER) {
						byteswapline(fshort, ncols * (N1 + 2 * dx) * (N2 + 2 * dy), pixsiz);
					}
					break;
				case FLOAT_PIXTYPE:
					ffloat = (float *) f;
					for (i = 0; i < ncols * (N1 + 2 * dx) * (N2 + 2 * dy); i++) {
						ffloat[i] = FLOAT_MAGIC;
					}
					if (fits->opbyteorder == NON_NATIVE_BYTE_ORDER) {
						byteswapline(ffloat, ncols * (N1 + 2 * dx) * (N2 + 2 * dy), pixsiz);
					}
					break;
			}
		}
		for (col = 0; col < ncols; col++) {
			plane = ncols * row + col;
			for (i = 0; i < N2; i++) {
				offset = (ncols * (i + dy) + col) * (nbytes + 2 * dx * pixsiz) + dx * pixsiz;
				if (gpc1format) {
					if ((row % 7) || (col % 7)) {
						if ((col > 3) && gpc1rotaterhs) {
							offset = (ncols * (N2 - (i - dy) - 1) + col) * (nbytes + 2 * dx * pixsiz) + dx * pixsiz;
							for (j = N1 - 1; j >= 0; j--) {
								fread(((char *) f) + offset + j * pixsiz, pixsiz, sizeof(char), stdin);
							}
						} else {
							fread(((char *) f) + offset, nbytes, sizeof(char), stdin);
						}
					} 
				} else {
					if (plane < nplanes) {
						fread(((char *) f) + offset, nbytes, sizeof(char), stdin);
					}
				}
			}
		}
		fwrite(f, ncols * (N2 + 2 * dy) * (nbytes + 2 * dx * pixsiz), sizeof(char), stdout);
		free(f);		
	}
	writefitstail(fits);
	exit(0);
}

