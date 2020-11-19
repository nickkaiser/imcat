/*
 * makeshear.c
 */

#define usage "\n\
NAME\n\
	maketestshear --- make a test shear field\n\
\n\
SYNOPSIS\n\
	maketestshear [options...]\n\
		-n N	 	# image size (256)\n\
		-o i0 j0	# lens position (128,128)\n\
		-c rc		# core radius (50)\n\
		-k		# generate kappa\n\
		-g		# generate gradkappa\n\
\n\
DESCRIPTION\n\
	'maketestshear' generates a float format fits image of shear\n\
	for kappa = exp(-r^2 / 2 r_c^2) model\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n"


#include <stdio.h>
#include <math.h>

#include "../../imlib/fits.h"
#include "../../utils/arrays.h"

#define EPS	1.e-10

#define	GAMMA_MODE	0
#define KAPPA_MODE	1
#define GRAD_KAPPA_MODE	2

main (int argc, char *argv[])
{
	int		N, i0, j0, arg = 1, i, j, mode;
	double		x, y, rc, Y, kappabar, kappaav, gammaT, r2;
	float		**gamma1, **gamma2, **kappa, **kappa1, **kappa2, kappasum = 0.0;
	fitsheader	*fits;

	/* defaults */
	N = 256;
	i0 = j0 = 128;
	rc = 50.0;
	mode = GAMMA_MODE;

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			fprintf(stderr, "%s", usage);
			exit(-1);
		}
		switch (argv[arg++][1]) {
			case 'n':
				sscanf(argv[arg++], "%d", &N);
				break;
			case 'o':
				sscanf(argv[arg++], "%d", &i0);
				sscanf(argv[arg++], "%d", &j0);
				break;
			case 'c':
				sscanf(argv[arg++], "%lf", &rc);
				break;
			case 'k':
				mode = KAPPA_MODE;
				break;
			case 'g':
				mode = GRAD_KAPPA_MODE;
				break;
			default:
				fprintf(stderr, "%s", usage);
				exit(-1);
				break;
		}
	}
	
	allocFloatArray(&gamma1, N, 2 * N);
	gamma2 = gamma1 + N;
	kappa1 = kappa = gamma1;
	kappa2 = gamma2;

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			x = (j - j0);
			y = (i - i0);
			Y = 0.5 * (x * x + y * y) / (rc * rc);
			kappaav = exp(-Y);
			switch (mode) {
				case KAPPA_MODE:
					kappa[i][j] = kappaav;
					kappasum += kappaav;
					break;
				case GAMMA_MODE:
					kappabar = (Y > 0.0 ? (1 - kappaav) / Y : 1.0);
					gammaT = kappabar - kappaav;
					r2 = x * x + y * y + EPS;
					gamma1[i][j] = -gammaT * (x * x - y * y) / r2;
					gamma2[i][j] = -gammaT * 2 * x * y / r2;
					break;
				case GRAD_KAPPA_MODE:
					kappa1[i][j] = -x * kappaav / (rc * rc);
					kappa2[i][j] = -y * kappaav / (rc * rc);
					break;
				default:
					fprintf(stderr, "makeshear: bad mode\n");
					exit(-1);
					break;
			}
		}
	}

	fits = new2Dfitsheader(N, N, FLOAT_PIXTYPE);
	add_comment(argc, argv, fits);
	if (mode == KAPPA_MODE) {
		kappasum /= (double) (N * N);
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				kappa[i][j] -= kappasum;
			}
		}
		fits->n[1] = N;
		write2Dfloatimage(kappa, fits);
		exit(0);
	} else {
		fits->ndim = 3;
		fits->n[2] = 2;
		writefitsheader(fits);
		if (mode == GAMMA_MODE) {
			for (i = 0; i < 2 * N; i++) {
				writefitsline(gamma1[i], fits);
			}
		} else {
			for (i = 0; i < 2 * N; i++) {
				writefitsline(kappa1[i], fits);
			}
		}
		writefitstail(fits);
	}
	exit(0);
}


#undef EPS


