#define usage "\n\n\n\
NAME\n\
	colorize --- convert gray fits file to rgb colormapped version\n\
\n\
SYNOPSIS\n\
	colorize colormapindex [options...] \n\
\n\
DESCRIPTION\n\
	If \"colorize\" reads a 2-dimensional fits file f[N2][N1] from standard\n\
	input and writes a 3-dimensional rgb version f[3][N2][N1] with\n\
	f[0][][] = r[][], f[1][][] = g[][], f[2][][] = b[][].\n\
\n\
	If it reads a 3 dimensional image f[N3][N2][N1] it generates a four\n\
	dimensional version f[N3][3][N2][N1] -- i.e. a stream of 3-D rgb\n\
	images.\n\
\n\
	The rgb colors are computed by linear interpolation on a\n\
	color ramp.\n\
\n\
	Options are:\n\
		-f fmin fmax	# limits for input image (0, 255)\n\
		-a		# use autoscaling\n\
		-p bitpix	# output pixtype (8)   \n\
\n\
	There are currently three colormaps available (0,1,2).\n\
\n\
BUGS\n\
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
#include "../utils/colormaps.h"
#include "../utils/arrays.h"

#define FRAC_MAX 0.9999999999999

int		main(int argc, char *argv[])	
{
	int		bitpixout, N1, N2, autorange, x, y, plane, nplanes;
	fitsheader	*fitsin, *fitsout;
	float		fmin, fmax, **f, *fout;
	char		*flag;
	float		*lp, *rp, *gp, *bp, *c[3], contrast, bright, frac;
	int		level, nlevels, color, colormapindex;

	/* defaults */
	autorange = 0;
	bitpixout = 8;
	fmin = 0.0;
	fmax = 255.0;

	/* parse args */
	argsinit(argc, argv, usage);
	if (nextargtype() == FLAG_ARG) {
		error_exit(usage);
	}
	colormapindex = getargi();
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'a':
				autorange = 1;
				break;
			case 'f':
				fmin = getargf();
				fmax = getargf();
				autorange = 0;
				break;
			case 'p':
				bitpixout = getargi();
				break;
			default:
				error_exit(usage);
		}
	}

	/* read the color map */
	getcolormap(&lp, &rp, &gp, &bp, &contrast, &bright, &nlevels, colormapindex);
	c[0] = rp;
	c[1] = gp;
	c[2] = bp;

	/* read the source image header */
	fitsin = readfitsheader(stdin);

	/* create and write the output image header */
	fitsout = copyfitsheader(fitsin);
	switch (fitsin->ndim) {
		case 2:
			fitsout->ndim = 3;
			fitsout->n[2] = 3;
			nplanes = 1;
			break;
		case 3:
			fitsout->ndim = 4;
			nplanes = fitsout->n[3] = fitsout->n[2];
			fitsout->n[2] = 3;
			break;
		default:
			error_exit("colorize: image must be 2 or 3 dimensional\n");
	}
	N1 = fitsin->n[0];
	N2 = fitsin->n[1];
	fitsout->extpixtype = bitpixout;
	add_comment(argc, argv, fitsout);
	writefitsheader(fitsout);
	

	allocFloatArray(&f, fitsin->n[0], fitsin->n[1]);
	readfitsplane((void *) f, fitsin);

	/* find min and max if necessary */
	if (autorange) {
		fmin = fmax = f[0][0];
		for (y = 0; y < N2; y++) {
			for (x = 0; x < N1; x++) {
				fmin = (f[y][x] < fmin ? f[y][x] : fmin);
				fmax = (f[y][x] > fmax ? f[y][x] : fmax);
			}
		}
	}
	if (fmin == fmax) {
		/* something wrong */
		fmin = 0.0;
		fmax = 256.0;
	}

	/* allocate space for output line */ 
	fout = (float *) calloc(N1, sizeof(float));

	for (plane = 0; plane < nplanes; plane++) {
		if (plane) {
			readfitsplane((void *) f, fitsin);
		}
		/* convert fvalues to level.fraction values */
		for (y = 0; y < N2; y++) {
			for (x = 0; x < N1; x++) {
				/* scale to unit range and clamp */
				f[y][x] = (fmax - f[y][x]) / (fmax - fmin);	
				if (f[y][x] > FRAC_MAX) {
					f[y][x] = FRAC_MAX;
				}
				if (f[y][x] < 0.0) {
                               	 f[y][x] = 0.0;
                        	}
				/* find the level we lie above */
				level = 0;
				while (level < nlevels - 1) {
					if (f[y][x] < lp[level + 1]) {
						break;
					}
					level++;
				}
				if (level >= nlevels - 1) {
					level--;
				}
				/* compute fraction of level gap */
				frac = (f[y][x] - lp[level]) / (lp[level + 1] - lp[level]);
				f[y][x] = (float) level + (frac < FRAC_MAX ? frac : FRAC_MAX);
			}
		}
		/* process the data */
		for (color = 0; color < 3; color++) {
			for (y = 0; y < N2; y++) {
				for (x = 0; x < N1; x++) {
					level = (int) floor(f[y][x]);
					fout[x] = c[color][level];
					frac = f[y][x] - level;
					fout[x] += frac * (c[color][level + 1] - c[color][level]);
					if ((bitpixout == 8) || (bitpixout == 16) || (bitpixout == 32)) {
        	                        	fout[x] *= 255.0;
	                       	 	}
				}
				writefitsline(fout, fitsout);
			}
		}
	}
	writefitstail(fitsout);

	/* all done */
	exit(0);
}

