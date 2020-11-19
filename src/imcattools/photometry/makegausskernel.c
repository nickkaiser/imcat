/*
 * makegausskernel.c - calculate the kernel for the polarizability
 */

#define usage "\n\n\
NAME\n\
	makegausskernel --- calculate the kernel for gaussian psf\n\
\n\
SYNOPSIS\n\
	makegausskernel alpha beta sigma rf [N1 N2]\n\
\n\
DESCRIPTION\n\
	'makegausskernel' computes kernel K_alpha beta for the polarisation\n\
	for gaussian psf.\n\
	rf is the scale length for the gaussian weight function.\n\
	sigma is the scale length for the psf.\n\
	By default N1 - N2 = 256.\n\
\n"

#include <stdio.h>
#include <math.h>
#include "../../fftlib/myfft.h"
#include "../../imlib/fits.h"
#include "../../utils/arrays.h"



main (int argc, char *argv[])
{
	double 	rf, rf2, r2, r[2], w, dwlm, ddwlm, sigma, sigma2;
	int	N1, N2;
	fitsheader	*fits;
	float	**K;
	int	alpha, beta, ix, iy, i, j, l, m;
	int	M[2][2][2] = {{{1,0}, {0,-1}}, {{0,1}, {1,0}}};
	int	delta[2][2] = {{1, 0}, {0, 1}};

	N1 = N2 = 256;

	if (argc < 5) {
		fprintf(stderr, usage);
		exit(-1);
	}
	sscanf(argv[1], "%d", &alpha);
	sscanf(argv[2], "%d", &beta);
	if (!goodindex(alpha) || ! goodindex(alpha)) {
		fprintf(stderr, "makegausskernel: indices must be 0 or 1\n");
		exit(-1);
	}		
	sscanf(argv[3], "%lf", &sigma);
	sigma2 = sigma * sigma;
	sscanf(argv[4], "%lf", &rf);
	rf2 = rf * rf;
	if (argc > 5) {
		sscanf(argv[5], "%d", &N1);
		sscanf(argv[6], "%d", &N2);
	}

	allocFloatArray(&K, N1, N2);
	for (iy = 0; iy < N2; iy++) {
		r[1] = iy - N2 / 2;
		for (ix = 0; ix < N1; ix++) {
			r[0] = ix - N1 / 2;
			r2 = r[0] * r[0] + r[1] * r[1];
			w = exp(-0.5 * r2 / rf2);
			for (i = 0; i < 2; i++) {
				for (j = 0; j < 2; j++) {
					for (l = 0; l < 2; l++) {
						for (m = 0; m < 2; m++) {
							dwlm = w * (delta[j][l] * r[m] + delta[j][m] * r[l] - 
								r[j] * r[l] * r[m] / rf2);
							ddwlm = w * (
								delta[j][l] * delta[i][m] +
								delta[j][m] * delta[i][l] -
								(
									delta[i][j] * r[l] * r[m] +
									delta[i][l] * r[j] * r[m] +
									delta[i][m] * r[j] * r[l]
								) / rf2
							) - r[i] * dwlm / rf2;
							K[iy][ix] += M[alpha][i][j] * M[beta][l][m] * 
								(r[i] * dwlm - sigma2 * ddwlm);
						}
					}
				}
			}	
		}
	}
	fits = new2Dfitsheader(N1, N2, FLOAT_PIXTYPE);
	add_comment(argc, argv, fits);
	write2Dfloatimage(K, fits);
	exit(0);
}




int	goodindex(int x)
{
	return (x == 0 || x == 1);
}
