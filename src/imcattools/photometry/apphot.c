/*
 * apphot.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "../../imlib/fits.h"
#include "apphot.h"

#define GC_MAX 		32
#define R_AP_MIN 	1
#define R_SOFTEN	0.3
#define P_MIN		5

#define MAGIC FLOAT_MAGIC

int	apphot(double *flux, double *rh, double *rql, double *rqu, double *rp, 
		double *nbad, double *fmax,
		double *x, double rap, double *fb0, double *dfb, 
		float **f, int N1, int N2)
{
	pixel	*pix;
	int	rapint, npmax, p, np, pplus;
	int	ix, iy, ix0, iy0, dx, dy, i;
	double	rx, ry, r, ff, rval[3] = {0.0, 0.0, 0.0}, flev[3] = {0.25, 0.50, 0.75}, dr;
	int	usesky;

	usesky = (fb0 == NULL ? 0 : 1);

	/* create a sufficiently large array to hold all the pixels */
	rapint = (int) ceil(rap);
	npmax = (1 + 2 * rapint) * (1 + 2 * rapint);
	pix = (pixel *) calloc(npmax, sizeof(pixel));

	ix0 = (int) floor(x[0]);
	iy0 = (int) floor(x[1]);

	np = 0;
	*nbad = 0.0;
	*fmax = 0.0;
	for (dy = -rapint; dy <= rapint; dy++) {
		iy = iy0 + dy;
		if (iy < 0 || iy >= N2) {
			*nbad += 1.0;
			continue;
		}
		for (dx = -rapint; dx <= rapint; dx++) {
			ix = ix0 + dx;
			if (ix < 0 || ix >= N1 || iy < 0 || iy >= N2) {
				*nbad += 1.0;
				continue;
			}
			if (f[iy][ix] == MAGIC) {
				*nbad += 1.0;
				continue;
			}
			/* distance to centre of pixel */
			rx = 0.5 + ix - x[0];
			ry = 0.5 + iy - x[1];
			r = sqrt(rx * rx + ry * ry + R_SOFTEN * R_SOFTEN);
			if (r <= rap) {
				ff = (usesky ? f[iy][ix] - *fb0 - dx * dfb[0] - dy * dfb[1] : (double) f[iy][ix]);
				pix[np].r = r;
				pix[np].f = ff;
				np++;
			}
			*fmax = (ff > *fmax ? ff : *fmax);
		}
	}
	if (np < 1) {
		*flux = *rh = *rql =*rqu = *rp = *fmax = 0.0;		
		return (0);
	}
	qsort((void *) pix, np, sizeof(pixel), (int *) pixcmp);

	/* cumulate the flux */
	pix[0].fcum = pix[0].f;
	pix[0].nucum = pix[0].nu = pix[0].fcum / pix[0].r;
	for (p = 1; p < np; p++) {
		pix[p].fcum = pix[p-1].fcum + pix[p].f;
		pix[p].nu = pix[p].fcum / pix[p].r;
		pix[p].nucum = pix[p-1].nucum + pix[p].nu;
	}
	for (p = 1; p < np; p++) {
		pix[p].nus = (pix[p].nucum - pix[p/2].nucum) / (p - p/2);
	}


	/* estimate the petrosian radius as first peak in nus which is also > pix[2n].nus*/
	*rp = 0.0;
	for (p = P_MIN; p < np - 1; p++) {
		if (*rp != 0.0) {
			break;
		}
		pplus = (2 * p >= np ? np - 1 : 2 * p);
		if ((pix[p].nus > pix[p+1].nu) && (pix[p].nus > pix[pplus].nu)) {
				*rp = pix[p].r;
		}
	}
	if (*rp == 0.0) {
		*rp = pix[np-1].r;
	}
	

/*
for (p = 0; p < np; p++) {
	fprintf(stdout, "%lf %lf %lf %lf %d\n", pix[p].r, pix[p].fcum, pix[p].nu, pix[p].nus, pix[p].r < *rp);
}
exit(0);

*/

	*flux = pix[np - 1].fcum;
	for (i = 0; i < 3; i++) {
		flev[i] *= *flux;
	}
	for (p = 1; p < np; p++) {
		for (i = 0; i < 3; i++) {
			if ((pix[p].fcum > flev[i]) && (rval[i] == 0.0)) {
				rval[i] = pix[p-1].r;
				if ((pix[p].fcum - pix[p-1].fcum) > 0.0) {
					dr = (pix[p].r - pix[p-1].r) * (flev[i] - pix[p-1].fcum) / (pix[p].fcum - pix[p-1].fcum);
				}
				if (dr >= 0.0)
					rval[i] += dr;
			}
		}
	}
	*rql = rval[0];
	*rh  = rval[1];
	*rqu = rval[2];
		
	free(pix);
	return (1);
}


double	rpetrosian(double *x, float **f, int N1, int N2, double *fb0, double *dfb)
{
	int	usesky;
	double	sumf[GC_MAX], sum1[GC_MAX], ff, rr, rp;
	int	ix0, iy0, ix, iy, dx, dy, ir, rnumax, rmax;

	usesky = (fb0 == NULL ? 0 : 1);


	/* calculate the growth curves */
	for (ir = 0; ir < GC_MAX; ir++)
		sumf[ir] = sum1[ir] = 0.0;
	ix0 = (int) floor(x[0]);
	iy0 = (int) floor(x[1]);
	for (dy = -GC_MAX; dy <= GC_MAX; dy++) {
		iy = iy0 + dy;
		if (iy < 0 || iy >= N2)
			continue;
		for (dx = -GC_MAX; dx <= GC_MAX; dx++) {
			ix = ix0 + dx;
			if (ix < 0 || ix >= N1 || f[iy][ix] == MAGIC)
				continue;
			rr = dx * dx + dy * dy;
			ir = (rr > 0 ? (int) floor(0.5 + sqrt(rr)) : 0);
			if (ir >= GC_MAX)
				continue;
			ff = (usesky ? f[iy][ix] - *fb0 - dx * dfb[0] - dy * dfb[1] : (double) f[iy][ix]);
			sum1[ir] += 1.0;
			sumf[ir] += ff;
		}
	}

	/* get rnumax */
	rnumax = GC_MAX;	
	for (ir = 1; ir < GC_MAX; ir++) {				/* cumulate sums */		
		sum1[ir] += sum1[ir - 1];			/* and find max-nu radius */
		sumf[ir] += sumf[ir - 1];
		if (sumf[ir] * sumf[ir] * sum1[ir-1] < sumf[ir-1] * sumf[ir-1] * sum1[ir]
				&& rnumax == GC_MAX) {
			rnumax = ir - 1;			/* peak-nu radius*/
		}
	}

	return ((double) rnumax);
}


int	pixcmp(pixel *pix1, pixel *pix2)
{
	return (pix1->r > pix2->r ? (1) : (pix1->r == pix2->r ? 0 : -1));
}
