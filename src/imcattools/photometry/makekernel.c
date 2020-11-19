/*
 * makekernel.c - calculate the kernel for the polarizability
 */

#define usage "\n\n\
NAME\n\
	makekernel --- calculate the kernel for the polarizability\n\
\n\
SYNOPSIS\n\
	makekernel [-w rf | -W alpha rf | -K alpha beta rf | -Z alpha beta rf | -R alpha rf]\n\
\n\
DESCRIPTION\n\
	'makekernel' reads a psf image from stdin and computes various\n\
	kernels for the polarisation\n\
\n\
	Options are:\n\
		-u			# print this message\n\
		-w rf			# compute w\n\
		-W alpha rf		# compute W_alpha\n\
		-K alpha beta rf	# compute K_alpha beta\n\
		-Z alpha beta rf	# compute Z_alpha beta\n\
		-R alpha rf		# compute R_alpha   \n\
	rf is the scale length for the gaussian weight function.\n\
\n\
	With -w, -W, -Z options, only the size of the psf image is used.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n"

#include <stdio.h>
#include <math.h>
#include "../../fftlib/myfft.h"
#include "../../imlib/fits.h"
#include "../../utils/arrays.h"


#define wMODE 0
#define KMODE 1
#define WMODE 2
#define RMODE 3
#define ZMODE 4


main (int argc, char *argv[])
{
	int	arg = 1, mode;
	double 	rf, rf2, r2, r[2], w, dwrr, dwrrr;
	int	N1, N2;
	fitsheader	*fits;
	float	**g, **A[2], **B[2], **C, **K, **W, **R, gsum;
	int	alpha, beta, ix, iy, i, j, l, m;
	int	M[3][2][2] = {{{1, 0}, {0, 1}}, {{1,0}, {0,-1}}, {{0,1}, {1,0}}};
	int	delta[2][2] = {{1, 0}, {0, 1}};
	fft_type	gk, Ck, Bk[2], Wk;

	if (argc < 2) {
		fprintf(stderr, usage);
		exit(-1);
	}

	while (arg < argc) {
		if (argv[arg][0] != '-') {
			fprintf(stderr, usage);
			exit(-1);
		}
		switch (argv[arg++][1]) {
			case 'w':
				sscanf(argv[arg++], "%lf", &rf);
				rf2 = rf * rf;
				mode = wMODE;
				break;
			case 'W':
				sscanf(argv[arg++], "%d", &alpha);
				if (!goodindex(alpha)) {
					fprintf(stderr, "makekernel: index must be 0, 1 or 2\n");
					exit(-1);
				}
				sscanf(argv[arg++], "%lf", &rf);
				rf2 = rf * rf;
				mode = WMODE;
				break;
			case 'K':
				sscanf(argv[arg++], "%d", &alpha);
				sscanf(argv[arg++], "%d", &beta);
				if (!goodindex(alpha) || !goodindex(beta)) {
					fprintf(stderr, "makekernel: index must be 0, 1 or 2\n");
					exit(-1);
				}
				sscanf(argv[arg++], "%lf", &rf);
				rf2 = rf * rf;
				mode = KMODE;
				break;
			case 'Z':
				sscanf(argv[arg++], "%d", &alpha);
				sscanf(argv[arg++], "%d", &beta);
				if (!goodindex(alpha) || !goodindex(beta)) {
					fprintf(stderr, "makekernel: index must be 0, 1 or 2\n");
					exit(-1);
				}
				sscanf(argv[arg++], "%lf", &rf);
				rf2 = rf * rf;
				mode = ZMODE;
				break;
			case 'R':
				sscanf(argv[arg++], "%d", &alpha);
				if (!goodindex(alpha)) {
					fprintf(stderr, "makekernel: index must be 0, 1 or 2\n");
					exit(-1);
				}
				sscanf(argv[arg++], "%lf", &rf);
				rf2 = rf * rf;
				mode = RMODE;
				break;
			default:
				fprintf(stderr, usage);
				exit(-1);
				break;
		}
	}


	/* read the psf */
	read2Dfloatimage(&g, &N1, &N2, &fits, stdin);
	/* normalise the psf to unit mean */
	gsum = 0.0;
	for (iy = 0; iy < N2; iy++) {
		for (ix = 0; ix < N1; ix++) {
			gsum += g[iy][ix];
		}
	}
	gsum /= (N1 * N2);
	for (iy = 0; iy < N2; iy++) {
		for (ix = 0; ix < N1; ix++) {
			g[iy][ix] /= gsum;
		}
	}

	if (mode == wMODE) {
		allocFloatArray(&W, N1, N2);
		for (iy = 0; iy < N2; iy++) {
			r[1] = iy - N2 / 2;
			for (ix = 0; ix < N1; ix++) {
				r[0] = ix - N1 / 2;
				r2 = r[0] * r[0] + r[1] * r[1];
				w = exp(-0.5 * r2 / rf2);
				W[iy][ix] = w;
			}
		}
		add_comment(argc, argv, fits);
		write2Dfloatimage(W, fits);
		exit(0);
	}

	if (mode == WMODE) {
		allocFloatArray(&W, N1, N2);
		for (iy = 0; iy < N2; iy++) {
			r[1] = iy - N2 / 2;
			for (ix = 0; ix < N1; ix++) {
				r[0] = ix - N1 / 2;
				r2 = r[0] * r[0] + r[1] * r[1];
				w = exp(-0.5 * r2 / rf2);
				for (l = 0; l < 2; l++) {
					for (m = 0; m < 2; m++) {
						W[iy][ix] += 0.5 * M[alpha][l][m] * w * r[l] * r[m];
					}
				}
				
			}
		}
		add_comment(argc, argv, fits);
		write2Dfloatimage(W, fits);
		exit(0);
	}

	if (mode == ZMODE) {
		allocFloatArray(&K, N1, N2);
		for (iy = 0; iy < N2; iy++) {
			r[1] = iy - N2 / 2;
			for (ix = 0; ix < N1; ix++) {
				r[0] = ix - N1 / 2;
				r2 = r[0] * r[0] + r[1] * r[1];
				w = exp(-0.5 * r2 / rf2);
				for (l = 0; l < 2; l++) {
					for (m = 0; m < 2; m++) {
						for (i = 0; i < 2; i++) {
							for (j = 0; j < 2; j++) {
								K[iy][ix] += M[alpha][l][m] * M[beta][i][j] * w * r[i] * r[j] * r[l] * r[m];
							}
						}
					}
				}
				
			}
		}
		add_comment(argc, argv, fits);
		write2Dfloatimage(K, fits);
		exit(0);
	}

	if (mode == RMODE) {
		allocFloatArray(&(A[0]), N1, N2);
		allocFloatArray(&(A[1]), N1, N2);
		allocFloatArray(&(B[0]), N1, N2);
		allocFloatArray(&(B[1]), N1, N2);
		allocFloatArray(&C, N1, N2);
		allocFloatArray(&R, N1, N2);

		alloc_fft(&gk, N1, N2);
		alloc_fft(&Ck, N1, N2);
		alloc_fft(&(Bk[0]), N1, N2);
		alloc_fft(&(Bk[1]), N1, N2);

		for (iy = 0; iy < N2; iy++) {
			r[1] = iy - N2 / 2;
			for (ix = 0; ix < N1; ix++) {
				r[0] = ix - N1 / 2;
				r2 = r[0] * r[0] + r[1] * r[1];
				w = exp(-0.5 * r2 / rf2);
				for (j = 0; j < 2; j++) {
					for (i = 0; i < 2; i++) {
						A[j][iy][ix] += M[alpha][i][j] * r[i];
					}
					B[j][iy][ix] = -w * r[j] / rf2;
					C[iy][ix] += A[j][iy][ix] * B[j][iy][ix];
				}
			}
		}

		forward_fft(g, N1, N2, gk);
		forward_fft(C, N1, N2, Ck);
		ccf(gk, Ck, N1, N2, C, N1 / 2, N2 / 2);
		for (j = 0; j < 2; j++) {
			forward_fft(g, N1, N2, gk);
			forward_fft(B[j], N1, N2, Bk[j]);
			ccf(gk, Bk[j], N1, N2, B[j], N1 / 2, N2 / 2);
		}
		for (iy = 0; iy < N2; iy++) {
			for (ix = 0; ix < N1; ix++) {
				for (j = 0; j < 2; j++) {
					R[iy][ix] += 2 * A[j][iy][ix] * B[j][iy][ix];
				}
				R[iy][ix] -= C[iy][ix];
			}
		}

		add_comment(argc, argv, fits);
		write2Dfloatimage(R, fits);
		exit(0);
	}

	allocFloatArray(&(A[0]), N1, N2);
	allocFloatArray(&(A[1]), N1, N2);
	allocFloatArray(&(B[0]), N1, N2);
	allocFloatArray(&(B[1]), N1, N2);
	allocFloatArray(&C, N1, N2);
	allocFloatArray(&K, N1, N2);

	alloc_fft(&gk, N1, N2);
	alloc_fft(&Ck, N1, N2);
	alloc_fft(&(Bk[0]), N1, N2);
	alloc_fft(&(Bk[1]), N1, N2);

	for (iy = 0; iy < N2; iy++) {
		r[1] = iy - N2 / 2;
		for (ix = 0; ix < N1; ix++) {
			r[0] = ix - N1 / 2;
			r2 = r[0] * r[0] + r[1] * r[1];
			w = exp(-0.5 * r2 / rf2);
			for (j = 0; j < 2; j++) {
				for (i = 0; i < 2; i++) {
					A[j][iy][ix] += M[beta][i][j] * r[i];
				}
				for (l = 0; l < 2; l++) {
					for (m = 0; m < 2; m++) {
						dwrr = w * (delta[j][l] * r[m] + delta[j][m] * r[l] -
							r[j] * r[l] * r[m] / rf2);
						B[j][iy][ix] += M[alpha][l][m] * dwrr;
					}
				}
				C[iy][ix] += A[j][iy][ix] * B[j][iy][ix];
			}
		}
	}

	forward_fft(g, N1, N2, gk);
	forward_fft(C, N1, N2, Ck);
	ccf(gk, Ck, N1, N2, C, N1 / 2, N2 / 2);
	for (j = 0; j < 2; j++) {
		forward_fft(g, N1, N2, gk);
		forward_fft(B[j], N1, N2, Bk[j]);
		ccf(gk, Bk[j], N1, N2, B[j], N1 / 2, N2 / 2);
	}
	for (iy = 0; iy < N2; iy++) {
		for (ix = 0; ix < N1; ix++) {
			for (j = 0; j < 2; j++) {
				K[iy][ix] += 2 * A[j][iy][ix] * B[j][iy][ix];
			}
			K[iy][ix] -= C[iy][ix];
			K[iy][ix] *= 0.5;
		}
	}

	add_comment(argc, argv, fits);
	write2Dfloatimage(K, fits);
	exit(0);

}




int	goodindex(int x)
{
	return (x == 0 || x == 1 || x == 2);
}
