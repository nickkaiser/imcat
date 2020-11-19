/*
 * painttriangle.c
 */

#include <stdio.h>
#include <math.h>
#include "average.h"


#define	MAX(a,b)	((a) > (b) ? (a) : (b))
#define	MIN(a,b)	((a) < (b) ? (a) : (b))

int	pointintriangle(double x, double y, double **xt);

void	painttriangle(float ***fim, int fdim, double x1, double x2, int Nx, double y1, double y2, int Ny, double **xt, double **ft, int favgmode)
{
	double		x, y, xmin, xmax, ymin, ymax, m[2][2], minv[2][2], detm;
	static double	*f0 = NULL, **gradf = NULL, *df01 = NULL, *df12 = NULL;
	int		i, j, ix, ix1, ix2, iy, iy1, iy2;

	/* get min and max x */
	xmin = MIN(MIN(xt[0][0], xt[1][0]), xt[2][0]);
	xmax = MAX(MAX(xt[0][0], xt[1][0]), xt[2][0]);
	ymin = MIN(MIN(xt[0][1], xt[1][1]), xt[2][1]);
	ymax = MAX(MAX(xt[0][1], xt[1][1]), xt[2][1]);

	/* allocate f0 (and gradf if necessary) */
	if (!f0) {
		f0 = (double *) calloc(fdim, sizeof(double));
		df01 = (double *) calloc(fdim, sizeof(double));
		df12 = (double *) calloc(fdim, sizeof(double));
		if (favgmode == F_INTERP) {
			gradf = (double **) calloc(fdim, sizeof(double *));
			for (i = 0; i < fdim; i++) {
				gradf[i] = (double *) calloc(2, sizeof(double));
			}
		}
	}

	/* compute f0 = mean, median or intercept (and gradient for interpolation) */
	if (favgmode == F_MEAN || favgmode == F_MEDIAN) {
		for (i = 0; i < fdim; i++) {
			f0[i] = (favgmode == F_MEAN ? mean(ft, i) : median(ft, i));
		}
	} else {
		m[0][0] = xt[0][0] - xt[1][0];
		m[0][1] = xt[0][1] - xt[1][1];
		m[1][0] = xt[1][0] - xt[2][0];
		m[1][1] = xt[1][1] - xt[2][1];
		detm = m[0][0] * m[1][1] - m[0][1] * m[1][0];
		minv[0][0] =  m[1][1] / detm;
		minv[0][1] = -m[0][1] / detm;
		minv[1][0] = -m[1][0] / detm;
		minv[1][1] =  m[0][0] / detm;
		for (i = 0; i < fdim; i++) {
			df01[i] = ft[0][i] - ft[1][i];
			df12[i] = ft[1][i] - ft[2][i];
			for (j = 0; j < 2; j++) {
				gradf[i][j] = minv[j][0] * df01[i] + minv[j][1] * df12[i];
			}
			f0[i] = ft[0][i] - gradf[i][0] * xt[0][0] - gradf[i][1] * xt[0][1];
		}
	}

	/* paint the pixels */
	ix1 = MAX(0,  (int) floor(Nx * (xmin - x1) / (x2 - x1)));
	ix2 = MIN(Nx, (int) ceil(Nx * (xmax - x1) / (x2 - x1)));
	iy1 = MAX(0,  (int) floor(Ny * (ymin - y1) / (y2 - y1)));
	iy2 = MIN(Ny, (int) ceil(Ny * (ymax - y1) / (y2 - y1)));
	for (ix = ix1; ix < ix2; ix++) {
		x = x1 + ix * (x2 - x1) / Nx;
		for (iy = iy1; iy < iy2; iy++) {
			y = y1 + iy * (y2 - y1) / Ny;
			if (pointintriangle(x, y, xt)) {
				for (i = 0; i < fdim; i++) {
					fim[i][iy][ix] = f0[i];
					if (favgmode == F_INTERP) {
						fim[i][iy][ix] += gradf[i][0] * x + gradf[i][1] * y;
					}
				}
			}
		}
	}
}

int	pointintriangle(double x, double y, double **xt)
{
	int i, j, c = 0;

	for (i = 0, j = 2; i < 3; j = i++) {
		if ((((xt[i][1]<=y) && (y<xt[j][1])) ||
			((xt[j][1]<=y) && (y<xt[i][1]))) &&
			(x < (xt[j][0] - xt[i][0]) * (y - xt[i][1]) / (xt[j][1] - xt[i][1]) + xt[i][0])) {
				c = !c;
		}
	}
	return c;
}

