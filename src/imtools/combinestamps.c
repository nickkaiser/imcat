#define usage "\n\n\n\
NAME\n\
	combinestamps --- take mean, median etc of a stack of images\n\
\n\
SYNOPSIS\n\
	combinestamps [options...]\n\
\n\
DESCRIPTION\n\
	'combinestamps' reads a set of Nim images in the form of a\n\
	Nim x N2 x N1 fits image (as created by 'makestamps' for\n\
	instance) and performs some kind of\n\
	average.  By default it takes the median. Options are:\n\
		-m	# take mean instead\n\
		-a clip	# take average clipping at +- clip * sigma   \n\
	Input images must have identical sizes.  Combinestamps is\n\
	functionally very similar to 'combineimages', but is unsuitable\n\
	for large images as it reads all of the images into memory.\n\
	However, it is therefore capable for reading very large\n\
	*numbers* of images unlike 'combinimages' which reads images\n\
	line by line and must therefore have a file open for each\n\
	input image.\n\
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
#include "combine_stuff.h"

#define	MEAN_MODE	0
#define MEDIAN_MODE	1
#define	ASG_MODE	2

int		main(int argc, char *argv[])	
{
	int		arg = 1, mode, first = 1, imno, x, y, nimages, pixcount;
	FILE		*ipf, *opf;
	int		N1, N2;
	fitsheader	*fits;
	float		***f, **fout, *pixlist, clip;

	/* defaults */
	mode = MEDIAN_MODE;

	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		}
		switch (argv[arg++][1]) {
			case 'm':
				mode = MEAN_MODE;
				break;
			case 'a':
				mode = ASG_MODE;
				sscanf(argv[arg++], "%f", &clip);
				break;
			default:
				error_exit(usage);
				break;
		}
	}
				
	/* read the header */
	fits = readfitsheader(stdin);
	if (fits->ndim != 3) {
		error_exit("combinestamp: non 3-dimensional fits image\n");
	}
	nimages = fits->n[2];
	N2 = fits->n[1];
	N1 = fits->n[0];

	/* allocate and read the image */
	f = (float ***) calloc(nimages, sizeof(float **));
	for (imno = 0; imno < nimages; imno++) {
		allocFloatArray(&(f[imno]), N1, N2);
		for (y = 0; y < N2; y++) {
			readfitsline(f[imno][y], fits);
		}
	}
	close(fits->ipstream);

	pixlist = (float *) calloc(nimages, sizeof(float));
	allocFloatArray(&fout, N1, N2);
	add_comment(argc, argv, fits);
	fits->ndim = 2;
	for (y = 0; y < N2; y++) {
		for (x = 0; x < N1; x++) {
			pixcount = 0;
			for (imno = 0; imno < nimages; imno++) {
				if (f[imno][y][x] != FLOAT_MAGIC) {
					pixlist[pixcount++] = f[imno][y][x];
				}
			}
			if (pixcount) {
				switch (mode) {
					case MEDIAN_MODE:
						fout[y][x] = median(pixlist, pixcount);
						break;
					case MEAN_MODE:
						fout[y][x] = mean(pixlist, pixcount);
						break;
					case ASG_MODE:
						fout[y][x] = avsigclip(pixlist, pixcount, clip);
						break;
					default:
						error_exit("combinestamps: bad mode\n");
				}
			} else {
				fout[x][y] = FLOAT_MAGIC;
			}
		}
	}
	write2Dfloatimage(fout, fits);	
	return (0);
}


