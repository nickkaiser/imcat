/*
 * stackpsfs.c --- read a fits stream of PSFs and average after recentering
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils/error.h"
#include "utils/arrays.h"
#include "utils/args.h"
#include "imlib/fits.h"


#define usage "\n\
NAME\n\
	stackpsfs --- read a fits stream of PSFs and average after recentering\n\
\n\
SYNOPSIS\n\
	stackpsfs [-n nplanes]\n\
\n\
DESCRIPTION\n\
\n\
	stackpsfs reads a fits stream of PSFs and averages after recentering.\n\
	It generates a 4 x N2 x N1 image consiting of 4 image planes:\n\
		0	# straight average\n\
		1	# center on centroid\n\
		2	# center on g^2 wighted centoid\n\
		3	# center on peak\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n\n"

int	shiftandadd(int dx, int dy, float **fsrc, int N, float **fdst);


main (int argc, char *argv[])
{
	int		np, p, ix, iy, N, ix0, iy0, ix1, iy1, ixmax, iymax;
	float		**g, **gnatl, **gcent, **gpeak, **ggcent;
	double		x0, y0, x1, y1, gmax, gsum, ggsum;
	fitsheader	*fits;
	char		*flag;	

	/* defaults */
	np = 0;

	/* parse args */
	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'n':
				np = getargi();
				break;
			default:
				error_exit(usage);
		}
	}

	fits = readfitsheader(stdin);
	if (fits->ndim != 3 || fits->n[0] != fits->n[1]) {
		error_exit("stackpsfs: bad source image dimensions\n");
	}
	if (np) {
		if (np > fits->n[2]) {
			error_exit("stackpsfs: source file contains too few planes\n");
		}
	} else {
		np = fits->n[2];
	}
	N = fits->n[0];

	/* allocation of arrays */
	allocFloatArray(&g, N, N);
	allocFloatArray(&gnatl, N, N);
	allocFloatArray(&gcent, N, N);
	allocFloatArray(&gpeak, N, N);
	allocFloatArray(&ggcent, N, N);

	/* accumulate with shifts */
	for (p = 0; p < np; p++) {
		readfitsplane((void *) g, fits);
		x0 = y0 = x1 = y1 = 0.0;
		gmax = gsum = ggsum = 0.0;
		for (iy = 0; iy < N; iy++) {
			for (ix = 0; ix < N; ix++) {
				/* centroid */
				x0 += ix * g[iy][ix];
				y0 += iy * g[iy][ix];
				gsum += g[iy][ix];
				/* gg centroid */
				x1 += ix * g[iy][ix] * g[iy][ix];
				y1 += iy * g[iy][ix] * g[iy][ix];
				ggsum += g[iy][ix] * g[iy][ix];						
				/* peak tracking */
				if (g[iy][ix] > gmax) {
					gmax = g[iy][ix];
					ixmax = ix;
					iymax = iy;
				}
			}
		}
		/* natl */
		shiftandadd(0, 0, g, N, gnatl);
		/* centroid */
		ix0 = (int) floor(0.5 + x0 / gsum) - N / 2;
		iy0 = (int) floor(0.5 + y0 / gsum) - N / 2;
		shiftandadd(-ix0, -iy0, g, N, gcent); 
		/* super centroid */
		ix0 = (int) floor(0.5 + x1 / ggsum) - N / 2;
		iy0 = (int) floor(0.5 + y1 / ggsum) - N / 2;
		shiftandadd(-ix0, -iy0, g, N, ggcent); 
		/* shift and add */
		ixmax -= N / 2;
		iymax -= N / 2;
		shiftandadd(-ixmax, -iymax, g, N, gpeak); 	
	}
	
	/* normalise and output */
	for (iy = 0; iy < N; iy++) {
		for (ix = 0; ix < N; ix++) {
			gnatl[iy][ix] /= np;
			gcent[iy][ix] /= np;
			ggcent[iy][ix] /= np;
			gpeak[iy][ix] /= np;
		}
	}
	fits->n[2] = 4;
	add_comment(argc, argv, fits);
	writefitsheader(fits);
	writefitsplane((void *) gnatl, fits);
	writefitsplane((void *) gcent, fits);
	writefitsplane((void *) ggcent, fits);
	writefitsplane((void *) gpeak, fits);
	writefitstail(fits);
	exit(0);	
}




int	shiftandadd(int dx, int dy, float **fsrc, int N, float **fdst)
{
	int 	x, y, xp, yp;

	for (x = 0; x < N; x++) {
		xp = x + dx;
		if (xp >= 0 && xp < N) {
			for (y = 0; y < N; y++) {
				yp = y + dy;
				if (yp >= 0 && yp < N) {
					fdst[yp][xp] += fsrc[y][x];
				}
			}
		}
	}
}


