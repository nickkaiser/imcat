/*
 * getsky.c
 */

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include "../../utils/error.h"
#include "../../imlib/fits.h"
#include "getsky.h"

#define	MAX_PIX		100000
#define	MIN_PIX		16

#define X_POS		2
#define	X_NEG		1
#define Y_POS		3
#define Y_NEG		0

#define MAGIC	FLOAT_MAGIC

static	int qlist[4];
static	double	xbar[4], ybar[4];
static	float	fmode[4];

int	getsky(double *fb0, double *dfb, double r1, double r2, double *x, 
		float **f, int N1, int N2)
{
	static	double 	dx[4][MAX_PIX], dy[4][MAX_PIX];
	static	float	ff[4][MAX_PIX];
	static	int	npix[4];
	int	ix, iy, dix, diy, dimax, pix;
	int	qp, qm, q;
	float   median, lquart, uquart, sigma;
	int	goodquads, iq;
	double 	r;

	ix = (int) floor(x[0]);
	iy = (int) floor(x[1]);

	/* shrink outer annulus radius if too large */
	if ((r2 * r2) > MAX_PIX)
		r2 = (int) floor(sqrt(MAX_PIX));
	dimax = 1 + r2;

	/* generate pix lists */
	for (q = 0; q < 4; q++)
		npix[q] = 0;
	for (diy = -dimax; diy <= dimax; diy++) {
		if (((iy + diy) < 0) || ((iy + diy) >= N2))
			continue;
		for (dix = -dimax; dix <= dimax; dix++) {
			if (((ix + dix) < 0) || ((ix + dix) >= N1))
				continue;
			r = dix * dix + diy * diy;
			if (r < r1 * r2 || r > r2 * r2)
				continue;
			qp = dix + diy;
			qm = dix - diy;
			if ((qp == 0) || (qm == 0))
				continue;
			q = 2 * (qp > 0) + (qm > 0);
			if (f[iy + diy][ix + dix] != MAGIC) {
				dx[q][npix[q]] = (double) dix;
				dy[q][npix[q]] = (double) diy;
				ff[q][npix[q]] = f[iy + diy][ix + dix];
				npix[q]++;
			}
			if (npix[q] >= MAX_PIX)
				error_exit("getsky: too many pixels in quadrant\n");
		}
	}

	/* calculate centroids, modes */
	iq = 0;
	goodquads = 4;
	for (q = 0; q < 4; q++) {
		if (npix[q] < MIN_PIX) {
			goodquads--;
			continue;
		} else {
			qlist[iq++] = q;
		}
		xbar[q] = ybar[q] = 0.0;
		for (pix = 0; pix < npix[q]; pix++) {
			xbar[q] += dx[q][pix];
			ybar[q] += dy[q][pix];
		}
		xbar[q] /= npix[q];
		ybar[q] /= npix[q];
		liststats(ff[q], npix[q], &(fmode[q]), &median, &lquart, &uquart, &sigma);
	}

	dfb[0] = dfb[1] = *fb0 = 0.0;
	switch (goodquads) {
		case 4:
			pruneextreme();
		case 3:
			getplane(fb0, dfb);
			break; 
		case 2:
			*fb0 = (double) (fmode[qlist[0]] + fmode[qlist[1]]) / 2.0; 
			break; 
		case 1: 
			*fb0 = (double) fmode[qlist[0]]; 
			break; 
		case 0:
			*fb0 = 0.0; 
			break;
		default:
			error_exit("getsky: bad case\n");
	}
	return (goodquads);
}


void	pruneextreme(void)
{
	int i;
	float	median;

	qsort(qlist, 4, sizeof(int), fmodecmp);
	median = 0.5 * (fmode[qlist[1]] + fmode[qlist[2]]);
	if ((median - fmode[qlist[0]]) > (fmode[qlist[3]] - median)) {
		qlist[0] = qlist[3];
	}
}



int	getplane(double *fb0, double *dfb)
{
	double	det, x[3], y[3], f[3], cof[3][3];
	int	i;

	for (i = 0; i < 3; i++) {
		x[i] = xbar[qlist[i]];
		y[i] = ybar[qlist[i]];
		f[i] = (double) (fmode[qlist[i]]);
	}

	cof[0][0] = x[1] * y[2] - x[2] * y[1];
	cof[0][1] = y[1] - y[2];
	cof[0][2] = x[2] - x[1];
	cof[1][0] = x[2] * y[0] - x[0] * y[2];
	cof[1][1] = y[2] - y[0];
	cof[1][2] = x[0] - x[2];
	cof[2][0] = x[0] * y[1] - x[1] * y[0];
	cof[2][1] = y[0] - y[1];
	cof[2][2] = x[1] - x[0];
	det = cof[0][0] + cof[1][0] + cof[2][0];
	if (det == 0.0) {
		return (0);
	}
	*fb0 = (cof[0][0] * f[0] + cof[1][0] * f[1] + cof[2][0] * f[2]) / det;
	dfb[0] = (cof[0][1] * f[0] + cof[1][1] * f[1] + cof[2][1] * f[2]) / det;
	dfb[1] = (cof[0][2] * f[0] + cof[1][2] * f[1] + cof[2][2] * f[2]) / det;
	return (1);
}



int	fmodecmp(const void *ptr1, const void *ptr2)
{
	int	*q1, *q2;

	q1 = (int *) ptr1;
	q2 = (int *) ptr2;
	if (fmode[*q1] > fmode[*q2]) {
		return (1);
	} else {
		if (fmode[*q1] < fmode[*q2]) {
			return (-1);
		} else {
			return (0);
		}
	}
}
