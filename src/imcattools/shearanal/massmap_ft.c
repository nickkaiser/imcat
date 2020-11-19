/*
 * massmap_ft.c
 */

#define usage "\n\
NAME\n\
	massmap_ft --- ks93 reconstruction using FFT\n\
\n\
SYNOPSIS\n\
	massmap_ft [options...]\n\
		-m mode	 	# mode (0)\n\
		-r R	 	# pad by factor R (2)\n\
\n\
DESCRIPTION\n\
	'massmap_ft' reads a shear image in the format produced by 'makeshearimage'\n\
	from stdin and reconstructs density field by fourier-space\n\
	version of KS93.\n\
	The input data are zero-padded onto an internal image which is\n\
	by default twice the size of the original.\n\
	You can use a slightly different estimator with -m option:\n\
		mode = 0 for KS93 chi estimator\n\
		mode = 1 for zeta estimator.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n"

#include <stdio.h>
#include <math.h>

#include "../../imlib/fits.h"
#include "../../utils/arrays.h"
#include "../../fftlib/myfft.h"

#define	EPS	1.e-10

float	chi1(float ki, float kj);
float	chi2(float ki, float kj);
float	zeta1(float ki, float kj);
float	zeta2(float ki, float kj);

#define CHI_MODE		0
#define ZETA_MODE		1


main (int argc, char *argv[])
{
	int		arg = 1, M, N, mode, i, j;
	float		**g_in, **g1, **g2, **gamma1, **gamma2, **kappa, **kappa1, **kappa2;
	fft_type 	Gamma1, Gamma2;
	int		R;
	fitsheader	*fits;

	/* defaults */
	mode = CHI_MODE;
	R = 2;

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			fprintf(stderr, "%s", usage);
			exit(-1);
		}
		switch (argv[arg++][1]) {
			case 'm':
				sscanf(argv[arg++], "%d", &mode);
				if (mode != CHI_MODE && mode != ZETA_MODE) {
					fprintf(stderr, "%s", usage);
					exit(-1);
				}
				break;
			case 'r':
				sscanf(argv[arg++], "%d", &R);
				break;
			default:
				fprintf(stderr, "%s", usage);
				exit(-1);
				break;
		}
	}

	fits = readfitsheader(stdin);
	N = fits->n[0];
	M = 2 * fits->n[1];
	allocFloatArray(&g_in, N, M);
	for (i = 0; i < M; i++) {
		readfitsline(g_in[i], fits);
	}
	g1 = g_in;
	g2 = g_in + N;

	allocFloatArray(&gamma1, R * N, R * N);
	allocFloatArray(&gamma2, R * N, R * N);
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			gamma1[i][j] = g1[i][j];
			gamma2[i][j] = g2[i][j];
		}
	}

	N *= R;

	alloc_fft(&Gamma1, N, N);
	alloc_fft(&Gamma2, N, N);
	allocFloatArray(&kappa, N, N);
	allocFloatArray(&kappa1, N, N);
	allocFloatArray(&kappa2, N, N);

	forward_fft(gamma1, N, N, Gamma1);
	forward_fft(gamma2, N, N, Gamma2);

	switch (mode) {
		case CHI_MODE:
			filter(Gamma1, N, N, chi1);
			filter(Gamma2, N, N, chi2);
			break;
		case ZETA_MODE:
			filter(Gamma1, N, N, zeta1);
			filter(Gamma2, N, N, zeta2);
			break;
		default:
			fprintf(stderr, "fourier: bad mode\n");
			break;
	}

	inverse_fft(Gamma1, N, N, kappa1);
	inverse_fft(Gamma2, N, N, kappa2);

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			kappa[i][j] = kappa1[i][j] + kappa2[i][j];
		}
	}

	N /= R;
	
	add_comment(argc, argv, fits);
	fits->ndim = 2;
	fits->n[0] = fits->n[1] = N;
	write2Dfloatimage(kappa, fits);
	exit(0);
}


float	chi1(float ki, float kj)
{
	float	kk;

	kk = ki * ki + kj * kj + EPS;
	return ((ki * ki - kj * kj) / kk);
}


float	chi2(float ki, float kj)
{
	float	kk;

	kk = ki * ki + kj * kj + EPS;
	return (2 * ki * kj / kk);
}


float	zeta1(float ki, float kj)
{
	float	fabski, fabskj;

	fabski = fabs(ki);
	fabskj = fabs(kj);

	if (fabski == 0.0 && fabskj == 0.0)
		return (0.0);
	else
		return (fabski > fabskj ? 1.0 : -1.0);
}


float	zeta2(float ki, float kj)
{
	float	fabski, fabskj;

	fabski = fabs(ki);
	fabskj = fabs(kj);

	if (fabski == 0.0 && fabskj == 0.0)
		return (0.0);
	else
		return (fabski > fabskj ? kj / ki : ki / kj);
}




#undef EPS



