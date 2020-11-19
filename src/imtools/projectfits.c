#define usage "\n\n\n\
NAME\n\
	projectfits - average over rows, cols etc of a FITS file\n\
\n\
SYNOPSIS\n\
	projectfits [-u] [-d M] d\n\
\n\
DESCRIPTION\n\
	\"projectfits\" reads a FITS image of arbitrary dimensionality D from\n\
	stdin and outputs a D-1 dimensional image to stdout which\n\
	contains the average along the d'th dimension, where\n\
	d=0 is the fastest direction (row average), d=1 is the next\n\
	fastest direction (column average) etc.\n\
\n\
	Options:\n\
		-u		# print this message\n\
		-d M		# deproject   \n\
\n\
	With the -d option we stretch out a N-dimensional image along\n\
	the d'th direction to make a N+1 dimensional image.  For example,\n\
	with a 3-D input image fin[Nz][Ny][Nx], the result of\n\
		projectfits -d 10 2\n\
	is to make a 4-dimensional image fout[Nz][10][Ny][Nx] with\n\
		fout[z][t][y][x] = fin[z][y][x].\n\
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

static	float	*fout, *fin;
static	int	*fsizein, *fsizeout, *xin, *xout, d, deproj, Ndeproj, nlines;
static	fitsheader	*fitsin, *fitsout;

int	project(int level);
int	deproject(int level);

int		main(int argc, char *argv[])	
{
	char		*flag;
	int		i;
	float		avfac;

	/* parse the args */
	deproj = 0;
	argsinit(argc, argv, usage);
	while (nextargtype() == FLAG_ARG) {
		flag = getflag();
		switch (flag[0]) {
			case 'd':
				deproj = 1;
				Ndeproj = getargi();
				break;
			default:
				error_exit(usage);
		}
	}
	d = getargi();	

	/* read the fits header */
	fitsin = readfitsheader(stdin);

	/* check the dimension is OK */
	if (!deproj) {
		if (d < 0 || d >= fitsin->ndim) {
			error_exit("projectfits: bad dimension d\n");
		}
	}

	/* create the output fits header */
	fitsout = copyfitsheader(fitsin);
	if (deproj) {
		fitsout->ndim++;
		for (i = 0; i < d; i++) {
			fitsout->n[i] = fitsin->n[i];
		}
		fitsout->n[d] = Ndeproj;
		for (i = d + 1; i < fitsout->ndim; i++) {
			fitsout->n[i] = fitsin->n[i-1];
		}
	} else {
		fitsout->ndim--;
		for (i = 0; i < d; i++) {
			fitsout->n[i] = fitsin->n[i];
		}
		for (i = d; i < fitsout->ndim; i++) {
			fitsout->n[i] = fitsin->n[i+1];
		}
	}
	add_comment(argc, argv, fitsout);
	writefitsheader(fitsout);

	/* compute the sizes of pixel, line, plane, cube etc of images */
	fsizein = (int *) calloc(fitsin->ndim + 1, sizeof(int));
	fsizein[0] = 1;
	for (i = 0; i < fitsin->ndim; i++) {
		fsizein[i+1] = fsizein[i] * fitsin->n[i];
	}
	fsizeout = (int *) calloc(fitsout->ndim + 1, sizeof(int));
	fsizeout[0] = 1;
	for (i = 0; i < fitsout->ndim; i++) {
		fsizeout[i+1] = fsizeout[i] * fitsout->n[i];
	}
	
	/* allocate the output image */
	fout = (float *) calloc(fsizeout[fitsout->ndim], sizeof(float));
	if (!fout) {
		error_exit("projectfits: failed to allocate memory for output image\n");
	}

	/* allocate the input image */
	fin = (float *) calloc(fsizein[fitsin->ndim], sizeof(float));
	if (!fin) {
		error_exit("projectfits: failed to allocate memory for input image\n");
	}
	
	/* allocate the pixel index arrays */
	xout = (int *) calloc(fitsout->ndim, sizeof(int));
	xin  = (int *) calloc(fitsin->ndim, sizeof(int));

	/* read the input image */
	nlines = 1;
	for (i = 1; i < fitsin->ndim; i++) {
		nlines *= fitsin->n[i];
	}
	for (i = 0; i < nlines; i++) {
		readfitsline((void *) (fin + i * fitsin->n[0]), fitsin);
	}

	if (deproj) {
		deproject(fitsin->ndim);
	} else {
		/* process the data recursively */
		project(fitsout->ndim);
		/* and normalise the average */
		avfac = 1.0 / fitsin->n[d];
		for (i = 0; i < fsizeout[fitsout->ndim]; i++) {
			fout[i] *= avfac;
		}
	}
	
	/* read the output image */
	nlines = 1;
	for (i = 1; i < fitsout->ndim; i++) {
		nlines *= fitsout->n[i];
	}
	for (i = 0; i < nlines; i++) {
		writefitsline((void *) (fout + i * fitsout->n[0]), fitsout);
	}

	writefitstail(fitsout);
	exit(0);
}


/* recursive function to loop over output pixels and project */
int	project(int level)
{
	int	i, offout, offin;
	float	ff;

	level--;
	for (xout[level] = 0; xout[level] < fitsout->n[level]; xout[level]++) {
		if (level < d) {
			xin[level] = xout[level];
		} else {
			xin[level + 1] = xout[level];
		}
		if (level) {
			project(level);
		} else {
			offout = 0;
			for (i = 0; i < fitsout->ndim; i++) {
				offout += xout[i] * fsizeout[i];
			} 
			offin = 0;
			for (i = 0; i < fitsin->ndim; i++) {
				if (i != d) {
					offin += xin[i] * fsizein[i];
				} 
			}
			for (i = 0; i < fitsin->n[d]; i++) {
				ff = fin[offin + i * fsizein[d]];
				if (ff != FLOAT_MAGIC) {
					fout[offout] += ff;
				}
			}
		}
	}
}	

/* recursive function to loop over input pixels and deproject */
int	deproject(int level)
{
	int	i, offout, offin;

	level--;
	for (xin[level] = 0; xin[level] < fitsin->n[level]; xin[level]++) {
		if (level < d) {
			xout[level] = xin[level];
		} else {
			xout[level + 1] = xin[level];
		}
		if (level) {
			deproject(level);
		} else {
			offin = 0;
			for (i = 0; i < fitsin->ndim; i++) {
				offin += xin[i] * fsizein[i];
			} 
			offout = 0;
			for (i = 0; i < fitsout->ndim; i++) {
				if (i != d) {
					offout += xout[i] * fsizeout[i];
				} 
			}
			for (i = 0; i < fitsout->n[d]; i++) {
				fout[offout + i * fsizeout[d]] = fin[offin];
			}
		}
	}
}	