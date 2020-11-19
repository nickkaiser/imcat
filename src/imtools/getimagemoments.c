#define usage "\n\n\n\
NAME\n\
	getimagemoments --- compute flux weighted moments of a FITS image\n\
\n\
SYNOPSIS\n\
	getimagemoments [option...] p lmax\n\
		-u	# print man page\n\
		-c	# use (N1/2, N2/2) as the centroid\n\
\n\
DESCRIPTION\n\
	getimagemoments computes flux weighted spatial moments up to order lmax.\n\
	It first calculates the centroid weighted by f[iy][ix] raised to the p'th power:\n\
\n\
		xbar = sum (ix + 0.5) pow(f[iy][ix], p) / sum pow(f[iy][ix], p)\n\
		ybar = sum (iy + 0.5) pow(f[iy][ix], p) / sum pow(f[iy][ix], p)\n\
\n\
	and then outputs first sum f[iy][ix] and then, for l = lmin through lmax and\n\
	for m = 0 through l:\n\
\n\
		xybar = sum x^m y^(l-m) f[iy][ix] / sum f[iy][ix]\n\
\n\
	where x = ix + 0.5 - xbar; y = iy + 0.5 - ybar and lmin = 2 if p = 1 and\n\
	lmin = 1 otherwise.\n\
\n\
AUTHOR\n\
	Nick Kaiser\n\
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

int		main(int argc, char *argv[])	
{
	int		N1, N2, nplanes = 1, iplane, ix, iy, p, l, m, lmin, lmax, computecentroid = 1;
	fitsheader	*fits;
	float		**f;
	double		fsum, x, y, xy, fxysum, fxsum, fysum, xbar, ybar;
	char		*flag;

 	argsinit(argc, argv, usage);
        while (nextargtype() == FLAG_ARG) {
                flag = getflag();
		switch (flag[0]) {
			case 'c':
				computecentroid = 0;
				break;
			default :
				error_exit(usage);
		}
	}
	p = getargi();
	lmax = getargi();	

	/* read the fits header */
	fits = readfitsheader(stdin);
	if (fits->ndim > 2) {
		nplanes = fits->n[2];
	}
	N1 = fits->n[0];
	N2 = fits->n[1];
	allocFloatArray(&f, N1, N2);
	for (iplane = 0; iplane < nplanes; iplane++) {
	readfitsplane(f, fits);
	/* read2Dfloatimage(&f, &N1, &N2, &fits, stdin);*/

	if (computecentroid) {
		/* compute the f**p weighted centroid */
		fsum = fxsum = fysum = 0.0;
		for (iy = 0; iy < N2; iy++) {
			y = iy + 0.5;
			for (ix = 0; ix < N1; ix++) {
				x = ix + 0.5;
				fsum += pow(f[iy][ix], p);
				fxsum += x * pow(f[iy][ix], p);
				fysum += y * pow(f[iy][ix], p);
			}
		}
		xbar = fxsum / fsum;
		ybar = fysum / fsum;
	} else {
		xbar = N1 / 2;
		ybar = N2 / 2;
	}

	fsum = 0.0;
	lmin = (p == 1 ? 2 : 1);
	for (l = lmin; l <= lmax; l++) {
		for (m = 0; m <= l; m++) {
			fxysum = 0.0;
		        for (iy = 0; iy < N2; iy++) {
		                y = iy + 0.5 - ybar;
                		for (ix = 0; ix < N1; ix++) {
                		        x = ix + 0.5 - xbar;
					if ((l == lmin) && (m == 0)) {
                        			fsum += f[iy][ix];
					}
					xy = pow(x, (double) m) * pow(y, (double) (l - m));
                        		fxysum += xy * f[iy][ix];
				}
			}
			if ((l == lmin) && (m == 0)) {
				fprintf(stdout, "%13.8le ", fsum);
			}
			fprintf(stdout, "%13.8lf ", fxysum / fsum);
                }
        }

			
	fprintf(stdout, "\n");
	}
	exit(0);
}




