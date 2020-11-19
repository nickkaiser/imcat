/*
 * Getshape.c
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
new version with q_alpha and R_alpha_beta = repsonse of q's
*/

#define R_MAX_FACTOR 	4
#define	TINY_R 		1.e-10
#define Q_DIM		3

int	Getshape(double *x, double *fb0, double *dfb, double rc, double *flux, 
		double *q, double **R, double **P, float **f, int N1, int N2)
{
	float		F0, Q[2][2], R4[2][2][2][2], P4[2][2][2][2], Z[2][2][2][2], denom;
	int		i0, j0, ii, jj, i, j, di, dj, rmax, l, m, a, b;
	double		r, dx[2];
	double		W, fc;
	double		M[Q_DIM][2][2] = {{{1,0},{0,1}}, {{1,0},{0,-1}}, {{0,1},{1,0}}};

	rmax = (int) ceil(rc * R_MAX_FACTOR);
	i0 = (int) floor(x[1]);
	j0 = (int) floor(x[0]);
	
	/* initialise moments */
	F0 = 0.0;
	for (l = 0; l < 2; l++) {
		for (m = 0; m < 2; m++) {
			Q[l][m] = 0.0;
			for (i = 0; i < 2; i++) {
				for (j = 0; j < 2; j++) {
					Z[l][m][i][j] = 0.0;
				}
			}
		}
	}

	/* check that the filter radius is positive */
	if (rc <= 0.0) {
		return(0);
	}

	/* accumulate moments F0, Q[2][2], Z[2][2][2][2] */
	for (ii = i0 - rmax; ii <= i0 + rmax; ii++) { 
		for (jj =j0 - rmax; jj <= j0 + rmax; jj++) {
			if (ii >= 0 && ii < N2 && jj >= 0 && jj < N1) {
				di = ii - i0;
				dj = jj - j0;
				dx[0] = jj + 0.5 - x[0];
				dx[1] = ii + 0.5 - x[1];
				r = sqrt((double) (dx[0] * dx[0] + dx[1] * dx[1]));
				if (r <= rmax) {
					if (f[ii][jj] != MAGIC) {
						W = exp(-0.5 * r * r / (rc * rc));
						if (fb0) {
							fc = f[ii][jj] - fb0[0] - dx[0] * dfb[0] - dx[1] * dfb[1];
						} else {
							fc = f[ii][jj];
						}
						fc /= *flux;
						F0 += W * fc;
						dx[0] /= rc;
						dx[1] /= rc;
						for (l = 0; l < 2; l++) {
							for (m = 0; m < 2; m++) {
								Q[l][m] += dx[l] * dx[m] * W * fc;
								for (i = 0; i < 2; i++) {
									for (j = 0; j < 2; j++) {
										Z[l][m][i][j] += dx[l] * dx[m] * dx[i] * dx[j] * W * fc;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	/* compute R4, P4 */
	for (l = 0; l < 2; l++) {
		for (m = 0; m < 2; m++) {
			for (i = 0; i < 2; i++) {
				for (j = 0; j < 2; j++) {
					R4[l][m][i][j] = 0.5 * Z[l][m][i][j];
					P4[l][m][i][j] = -Z[l][m][i][j];
					if (j == l && i == m) R4[l][m][i][j] += 0.5 * F0;
					if (i == l && j == m) R4[l][m][i][j] += 0.5 * F0;
					if (j == l) R4[l][m][i][j] -= 0.5 * Q[m][i];
					if (j == m) R4[l][m][i][j] -= 0.5 * Q[l][i];
					if (i == l) R4[l][m][i][j] -= 0.5 * Q[m][j];
					if (i == m) R4[l][m][i][j] -= 0.5 * Q[l][j];
					if (i == j) R4[l][m][i][j] -= 0.5 * Q[l][m];
					if (j == l) P4[l][m][i][j] += Q[m][i];
					if (j == m) P4[l][m][i][j] += Q[l][i];
					if (i == j) P4[l][m][i][j] += Q[l][m];
				}
			}
		}
	}

	/* now compute contractions with M-matrixes */
	for (a = 0; a < Q_DIM; a++) {
		q[a] = 0.0;
		for (i = 0; i < 2; i++) {
			for (j = 0; j < 2; j++) {
				q[a] += M[a][i][j] * Q[i][j] * rc * rc;
			}
		}
		for (b = 0; b < Q_DIM; b++) {
			R[a][b] = P[a][b] = 0.0;
			for (i = 0; i < 2; i++) {
				for (j = 0; j < 2; j++) {
					for (l = 0; l < 2; l++) {
						for (m = 0; m < 2; m++) {
							R[a][b] += 0.5 * M[a][l][m] * M[b][i][j] * R4[l][m][i][j];
							P[a][b] += 0.5 * M[a][l][m] * M[b][i][j] * P4[l][m][i][j];
						}
					}
				}
			}
			P[a][b] *= rc * rc;
		}
	}
}





