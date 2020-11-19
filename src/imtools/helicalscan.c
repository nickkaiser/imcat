#define usage "\n\n\n\
NAME\n\
	helicalscan --- generate a helical scan of a 2D FITS image\n\
SYNOPSIS\n\
	helicalscan M d phi [-N nframes] [-o xo yo]\n\
\n\
DESCRIPTION\n\
	'helicalscan' reads a 2 dimensional FITS image fin[y][x] from stdin and\n\
	writes to stdout a 3 dimensional FITS image fout[i][y][x] consisting\n\
	of a set of M x M subimages with origin placed at a sequence of\n\
	positions on the input image:\n\
\n\
		x0 = (xo + i * d * cos(phi)) % N\n\
		y0 = (yo + i * d * sin(phi)) % N\n\
\n\
	The angle phi at which the subimage moves across the source image\n\
	is given in degrees.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "utils/error.h"
#include "utils/args.h"
#include "utils/ipbuff.h"
#include "imlib/fits.h"

#define MIN(a,b) ((a)<(b)?(a):(b))

int		main(int argc, char *argv[])	
{
	double		d, phi, dx, dy, xo, yo, xmin, ymin, xoff, yoff;
	int		x0, y0, x, y, i, M, N1, N2, nframes;
	float		**fin, **fout;
	fitsheader	*fits;
	char		*flag;

	/* set defaults */
	xo = yo = 0.0;
	nframes = 999999;

	/* parse args */
	argsinit(argc, argv, usage);
	if (nextargtype() == FLAG_ARG) {
		error_exit(usage);
	}
	M = getargi();
	d = getargd();
	phi = getargd();
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'N':
				nframes = getargi();
				break;
			case 'o':
				xo = getargd();
				yo = getargd();
				break;
			default:
				error_exit(usage);
		}
	}
	
	phi *= M_PI / 180.0;
	dx = d * cos(phi);
	dy = d * sin(phi);

	read2Dfloatimage(&fin, &N1, &N2, &fits, stdin);

	fits->ndim = 3;
	fits->n[0] = M;
	fits->n[1] = M;
	fits->n[2] = nframes;
	add_comment(argc, argv, fits);
	writefitsheader(fits);

	allocFloatArray(&fout, M, M);

	/* compute xoff, yoff which are multiples of N1, N2 and which */
	/* ensure that 1st arg of % operator is always positive */
	xmin = MIN(xo, xo + nframes * dx);
	ymin = MIN(yo, yo + nframes * dy);
	xoff = (xmin < 0 ? N1 * (1 - floor(xmin / N1)) : 0);
	yoff = (ymin < 0 ? N2 * (1 - floor(ymin / N2)) : 0);
 	
	for (i = 0; i < nframes; i++) {
		x0 = (int) floor(xoff + xo + i * dx);
		y0 = (int) floor(yoff + yo + i * dy);
		for (y = 0; y < M; y++) {
			for (x = 0; x < M; x++) {
				fout[y][x] = fin[(y0 + y) % N2][(x0 + x) % N1];
			}
		}
		writefitsplane((void *) fout, fits);			
	}
	exit(0);
}

