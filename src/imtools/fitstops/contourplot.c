/*
 * contourplot.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../imlib/fits.h"
#include "psutils.h"
#include "contourplot.h"

/* we break pixel into 4 triangles bounded by 2 adjacent pixel corners and the centre of the pixel */

void	contourplot(fitsheader *fits, int N1, int N2, int ncont, float *fcont, float width, float height)
{
	/* corners */
	static int	xc[2][4] = {{0,1, 1, 0}, {1, 1, 0, 0}};
	static int	yc[2][4] = {{0,0, 1, 1}, {0, 1, 1, 0}};
	int	ix, iy, ic, icont, leg, p;
	float	*ff[2], *ftemp, f[3];
	double	x[3], y[3], xp[2], yp[2];

	/* coordinates of centre never change */
	x[2] = y[2] = 0.5;

	/* allocate space for two lines */
	ff[0] = (float *) calloc(N1, sizeof(float));
	ff[1] = (float *) calloc(N1, sizeof(float));

	/* read the first line into ff[1] */
	readfitsline(ff[1], fits);

	for (iy = 1; iy < N2; iy++) {
		/* swap ff[1] -> ff[0] and read ff[1] */
		ftemp = ff[0];
		ff[0] = ff[1];
		ff[1] = ftemp;
		readfitsline(ff[1], fits);
		for (ix = 0; ix < N1 - 1; ix++) {
			for (icont = 0; icont < ncont; icont++) {
				/* compute the value of centre pixel */
				f[2] = 0.25 * (ff[0][ix] + ff[0][ix + 1] + ff[1][ix] + ff[1][ix + 1]) - fcont[icont];
				for (leg = 0; leg < 4; leg++) {
					/* compute values of corners */
					for (ic = 0; ic < 2; ic++) {
						x[ic] = xc[ic][leg];
						y[ic] = yc[ic][leg];
						f[ic] = ff[yc[ic][leg]][ix + xc[ic][leg]] - fcont[icont];
					}
					p = 0;
					if (f[0] * f[1] < 0.0) {
						xp[p] = (double) ((f[0] * x[1] - f[1] * x[0]) / (f[0] - f[1]));
						yp[p++] = (double) ((f[0] * y[1] - f[1] * y[0]) / (f[0] - f[1]));
					}
					if (f[1] * f[2] < 0.0) {
						xp[p] = (double) ((f[1] * x[2] - f[2] * x[1]) / (f[1] - f[2]));
						yp[p++] = (double) ((f[1] * y[2] - f[2] * y[1]) / (f[1] - f[2]));
					}
					if (f[2] * f[0] < 0.0) {
						xp[p] = (double) ((f[2] * x[0] - f[0] * x[2]) /  (f[2] - f[0]));
						yp[p++] = (double) ((f[2] * y[0] - f[0] * y[2]) / (f[2] - f[0]));
					}
					if (p == 2) {
						for (p = 0; p < 2; p++) {
							xp[p] = (0.5 + ix + xp[p]) * width / (N1 - 1);
							yp[p] = (0.5 + iy + yp[p]) * height / (N2 - 1);
						}
						psline(xp[0], yp[0], xp[1], yp[1]);
					}
				}
			}
		}
	}	
}


