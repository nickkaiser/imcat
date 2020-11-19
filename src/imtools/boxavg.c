#define usage "\n\n\n\
NAME\n\
	boxavg - boxcar smooth FITS image of arbitrary dimensionality\n\
\n\
SYNOPSIS\n\
	boxavg [-u] [-m] [-d d0 d1...]\n\
\n\
DESCRIPTION\n\
	\"boxavg\" reads a FITS image of arbitrary dimensionality from\n\
	stdin and applies a simple boxcar smoothing. By default the width\n\
	of the boxcar is 3 in each dimension, but you can specify an array\n\
	d[] with the -d option such that the width in the i'th direction\n\
	is 1 + 2 * d[i].\n\
\n\
	Options:\n\
		-u		# print this message\n\
		-d d0 d1...	# specify (width-1)/2 of box car (1,1,...)\n\
		-m		# take median\n\
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
#include "../utils/fmedian.h"

static float		*f0, *f1, fsum, *flist;
static int		*x0, *x1, *dx, *fsize, domedian, nsum;
static fitsheader	*fits;

void	getpix(int level);
void	accumulate(int level);

int		main(int argc, char *argv[])	
{
	char		*flag = NULL;
	int		i, npixinbox, nlines;

	/* defaults */
	domedian = 0;

	/* check for -u argument */
	argsinit(argc, argv, usage);
	if (nextargtype() == FLAG_ARG) {
		flag = getflag();
		if (flag[0] == 'u') {
			error_exit(usage);
		}
	}	

	/* read the fits header and allocate space for x0, x1, dx, fsize */
	fits = readfitsheader(stdin);
	x0 = (int *) calloc(fits->ndim, sizeof(int));
	x1 = (int *) calloc(fits->ndim, sizeof(int));
	dx = (int *) calloc(fits->ndim, sizeof(int));
	fsize = (int *) calloc(fits->ndim + 1, sizeof(int));

	/* initialise dx for 3 x 3 x 3 x ... smoothing */
	for (i = 0; i < fits->ndim; i++) {
		dx[i] = 1;
	}

	/* reinitialise the args machinery */
	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch(flag[0]) {
			case 'd':
				for (i = 0; i < fits->ndim; i++) {
					dx[i] = getargi();
				}
				break;
			case 'm':
				domedian = 1;
				break;
			default:
				error_exit(usage);
		}
	} 
	
	/* compute fsize[] = size (in floats) of pixel, line, cube.... */
	fsize[0] = 1;
	for (i = 1; i <= fits->ndim; i++) {
		fsize[i] = fsize[i-1] * fits->n[i-1];
	}

	/* allocate input and output arrays */
	f0 = (float *) calloc(fsize[fits->ndim], sizeof(float));
	f1 = (float *) calloc(fsize[fits->ndim], sizeof(float));
	if ((!f0) || (!f1)) {
		error_exit("boxavg: failed to allocate memory for image data\n");
	}

	/* read the data */
	nlines = 1;
        for (i = 1; i < fits->ndim; i++) {
                nlines *= fits->n[i];
        }
        for (i = 0; i < nlines; i++) {
                readfitsline((void *) (f0 + i * fits->n[0]), fits);
        }

	/* compute number of pixels in the smoothing box */
	npixinbox = 1;
	for (i = 0; i < fits->ndim; i++) {
		npixinbox *= (1 + 2 * dx[i]);
	}
	/* and allocate flist */
	flist = (float *) calloc(npixinbox, sizeof(float));

	/* recursively find all pixels and run accumulate() */
	getpix(fits->ndim);

	add_comment(argc, argv, fits);
	writefitsheader(fits);

	/* write the output image */
        for (i = 0; i < nlines; i++) {
                writefitsline((void *) (f1 + i * fits->n[0]), fits);
        }

	writefitstail(fits);
	exit(0);
}


/* recursive function to loop over all pixel indices x1[] */
void	getpix(int level)
{
	int	x1off, i;

	level--;
	for (x1[level] = 0; x1[level] < fits->n[level]; x1[level]++) {
		if (level) {
			getpix(level);
		} else {
			fsum = 0.0;
			nsum = 0;
			accumulate(fits->ndim);
			x1off = 0;
			for (i = 0; i < fits->ndim; i++) {
				x1off += fsize[i] * x1[i];
			}
			if (nsum) {
				if (domedian) {
					f1[x1off] = fmedian(flist, nsum);
				} else {
					f1[x1off] = fsum / nsum;
				}
			}
		}
	}
}

/* recursive function to loop over all x0[] within smoothing box and sum pixel values to make fout */
void	accumulate(int level)
{
	int	x0off, i;

	level--;
	for (x0[level] = x1[level] - dx[level]; x0[level] <= x1[level] + dx[level]; x0[level]++) {
		if (x0[level] >= 0 && x0[level] < fits->n[level]) {
			if (level) {
				accumulate(level);
			} else {
				/* compute indices for fin value */
				x0off = 0;
				for (i = 0; i < fits->ndim; i++) {
					x0off += fsize[i] * x0[i];
				}
				if (f0[x0off] != FLOAT_MAGIC) {
					if (domedian) {
						flist[nsum] = f0[x0off];
					} else {
						fsum += f0[x0off];
					}
					nsum += 1;
				}
			}
		}
	}
}

