/*
 * getshape2.c
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "getshape2.h"
#include "../../utils/error.h"
#include "../../imlib/fits.h"

#define MAGIC FLOAT_MAGIC

int	getshape2(double *x, double *fb0, double *dfb, double *F, double *p, double **P, double *R,
	float **f, int N1, int N2, float **w, float ***W, float ***RR, float ****K, int M1, int M2)
{
	int	m1, m2, ix0, iy0, ix, iy, kx, ky, l, m, dx, dy;
	float	fc;

	m1 = M1 / 2;
	m2 = M2 / 2;
	if ((2 * m1 != M1) || (2 * m2 != M2)) {
		error_exit("getshape2: psf file must have even dimensions\n");
	}

	/* we work to nearest pixel precision */
	ix0 = (int) floor(x[0]);
	iy0 = (int) floor(x[1]);

	/* initialise to zero */
	*F = 0.0;
	for (m = 0; m < 2; m++) {
		R[m] = 0;
	}
	for (l = 0; l < 3; l++) {
		p[l] = 0.0;
		for (m = 0; m < 2; m++) {
			P[l][m] = 0.0;
		}
	}

	/* accumulate */
	for (ky = 0; ky < M2; ky++) {
		iy = iy0 + ky - m2;
		for (kx = 0; kx < M1; kx++) {
			ix = ix0 + kx - m1;
			if (ix < 0 || ix >= N1 || iy < 0 || iy >= N2)
				continue;
			if (f[iy][ix] != MAGIC) {
				if (fb0) {
					dx = ix - ix0;
					dy = iy - iy0;
					fc = f[iy][ix] - fb0[0] - dx * dfb[0] - dy * dfb[1];
				} else {
					fc = f[iy][ix];
				}
				*F += w[ky][kx] * fc;
				for (m = 0; m < 2; m++) {
					R[m] += RR[m][ky][kx] * fc;
				}
				for (l = 0; l < 3; l++) {
					p[l] += W[l][ky][kx] * fc;
					for (m = 0; m < 2; m++) {
						P[l][m] += K[l][m][ky][kx] * fc;
					}
				}
			}
		}
	}

	/* normalise */
	if (*F != 0.0) {
		for (m = 0; m < 2; m++) {
			R[m] /= *F;
		}
		for (l = 0; l < 3; l++) {
			p[l] /= *F;
			for (m = 0; m < 2; m++) {
				P[l][m] /= *F;
			}
		}
	}
	return (1);
}





