/*
 * getshape.c
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "getshape.h"
#include "../../utils/error.h"
#include "../../imlib/fits.h"
#define MAGIC FLOAT_MAGIC

/*
7/31/97 - fixed bug in KSB
8/5/97 - non-integer pixel stuff
*/

#define R_MAX_FACTOR 4

int	getshape(double *x, double *fb0, double *dfb, double rc, double *flux, double *d, double *e, 
		double **psm, double **psh, float **f, int N1, int N2)
{
	float		q[2][2], denom;
	int		i0, j0, i, j, di, dj, rmax, l, m, negflux;
	double		r, dx, dy;
	double		W, Wp, Wpp, fc, DD, DD1, DD2;
	double		Xsm[2][2], Xsh[2][2], em[2], eh[2];
	double		temp;

	rmax = (int) ceil(rc * R_MAX_FACTOR);
	i0 = (int) floor(x[1]);
	j0 = (int) floor(x[0]);
	for (l = 0; l < 2; l++) {
		eh[l] = em[l] = 0.0;
		e[l] = 0.0;
		d[l] = 0.0;
		for (m = 0; m < 2; m++) {
			Xsh[l][m] = Xsm[l][m] = q[l][m] = 0.0;
			psm[l][m] = psh[l][m] = 0.0;
		}
	}
	if (rc <= 0.0) {
		return(0);
	}
	for (i =i0 - rmax; i <= i0 + rmax; i++) { 
		for (j =j0 - rmax; j <= j0 + rmax; j++) {
			if (i >= 0 && i < N2 && j >= 0 && j < N1) {
				di = i - i0;
				dj = j - j0;
				dx = j + 0.5 - x[0];
				dy = i + 0.5 - x[1];
/*
				dx = dj;
				dy = di;
*/
				r = sqrt((double) (dx * dx + dy * dy));
				if (r <= rmax) {
					if (f[i][j] != MAGIC) {
						W = exp(-0.5 * r * r / (rc * rc));
						Wp = -0.5 * W / (rc * rc);
						Wpp = 0.25 * W / (rc * rc * rc * rc);
						if (fb0) {
							fc = f[i][j] - fb0[0] - dx * dfb[0] - dy * dfb[1];
						} else {
							fc = f[i][j];
						}
						d[0] += W * fc * dx;
						d[1] += W * fc * dy;
						q[0][0] += fc * W * dx * dx;
						q[1][1] += fc * W * dy * dy;
						q[0][1] += fc * W * dx * dy;
						q[1][0] += fc * W * dx * dy;
						DD = di * di + dj * dj;
						DD1 = dj * dj - di * di;
						DD2 = 2 * di * dj;
						Xsm[0][0] += (2 * W + 4 * Wp * DD + 2 * Wpp * DD1 * DD1) * fc; 
						Xsm[1][1] += (2 * W + 4 * Wp * DD + 2 * Wpp * DD2 * DD2) * fc; 
						Xsm[0][1] += 2 * Wpp * DD1 * DD2 * fc; 
						Xsm[1][0] += 2 * Wpp * DD1 * DD2 * fc;
						em[0] += (4 * Wp + 2 * Wpp * DD) * DD1 * fc;
						em[1] += (4 * Wp + 2 * Wpp * DD) * DD2 * fc;
						Xsh[0][0] += (2 * W * DD + 2 * Wp * DD1 * DD1) * fc;
						Xsh[1][1] += (2 * W * DD + 2 * Wp * DD2 * DD2) * fc;
						Xsh[0][1] += 2 * Wp * DD1 * DD2 * fc;
						Xsh[1][0] += 2 * Wp * DD1 * DD2 * fc;
						eh[0] += 2 * Wp * DD * DD1 * fc;
						eh[1] += 2 * Wp * DD * DD2 * fc;
					}
				}
			}
		}
	}
	/* normalise d */
	if (*flux > 0) {
		for (l = 0; l < 2; l++) {
			d[l] /= *flux;
		}
		negflux = 0;
	} else {
		negflux = 1;
	}
	/* calculate ellipticities */
	denom = q[0][0] + q[1][1];			/* can change to sqrt(det) if desired */
	if (denom > 0) {
		e[0] = (q[0][0] - q[1][1]) / denom;
		e[1] = (q[0][1] + q[1][0]) / denom;
		em[0] /= denom;
		em[1] /= denom;
		eh[0] /= denom;
		eh[1] /= denom;
		eh[0] += 2 * e[0];
		eh[1] += 2 * e[1];
		for (l = 0; l < 2; l++) {
			for (m = 0; m < 2; m++) {
				/* factor 2 fix goes in here */
				psm[l][m] = 0.5 * (Xsm[l][m] / denom - e[l] * em[m]);
				psh[l][m] = Xsh[l][m] / denom - e[l] * eh[m];
			}
		}
		return (negflux ? 0 : 1);
	} else {
		return (0);
	}
}





