/*
 * getshape3.c
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "getshape3.h"
#include "../../utils/error.h"
#include "../../imlib/fits.h"

#define MAGIC FLOAT_MAGIC

/** this is changed wrt getshape2 in that we use the smoothed image to measure polarisation **/

int	getshape3(double *x, double *fb0, double *dfb, double *F, double *q0, double *q, double *P0, double **P, double **Z, double *R,
	float **f, float **fs, int N1, int N2, float **w, float ***W, float ***RR, float ****K, float ****ZZ, int M1, int M2)
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
	*q0 = 0.0;
	for (l = 0; l < 2; l++) {
		R[l] = 0.0;
		q[l] = 0.0;
		P0[l] = 0.0;
		for (m = 0; m < 2; m++) {
			P[l][m] = 0.0;
			Z[l][m] = 0.0;
		}
	}

	/* accumulate */
	for (ky = 0; ky < M2; ky++) {
		iy = iy0 + ky - m2;
		for (kx = 0; kx < M1; kx++) {
			ix = ix0 + kx - m1;
			if (ix < 0 || ix >= N1 || iy < 0 || iy >= N2)
				continue;
			if (fs[iy][ix] != MAGIC) {
                                if (fb0) {
                                        dx = ix - ix0;
                                        dy = iy - iy0;
                                        fc = fs[iy][ix] - fb0[0] - dx * dfb[0] - dy * dfb[1];
                                } else {
                                        fc = fs[iy][ix];
                                }
                                *F += w[ky][kx] * fc;
				*q0 += W[0][ky][kx] * fc;
                                for (l = 0; l < 2; l++) {
                                        q[l] += W[1+l][ky][kx] * fc;
					for (m = 0; m < 2; m++) {
						Z[l][m] += ZZ[1+l][m][ky][kx] * fc;
					}
                                }
                        }
			if (f[iy][ix] != MAGIC) {
				if (fb0) {
					dx = ix - ix0;
					dy = iy - iy0;
					fc = f[iy][ix] - fb0[0] - dx * dfb[0] - dy * dfb[1];
				} else {
					fc = f[iy][ix];
				}
				for (l = 0; l < 2; l++) {
					R[l] += RR[l][ky][kx] * fc;
					P0[l] += K[0][l][ky][kx] * fc;
					for (m = 0; m < 2; m++) {
						P[l][m] += K[1+l][m][ky][kx] * fc;
					}
				}
			}
		}
	}

	/* normalise */
	if (*F != 0.0) {
		*q0 /= *F;
		for (l = 0; l < 2; l++) {
			R[l] /= *F;
			q[l] /= *F;
			P0[l] /= *F;
			for (m = 0; m < 2; m++) {
				P[l][m] /= *F;
				Z[l][m] /= *F;
			}
		}
	}
	return (1);
}





