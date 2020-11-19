/*
 * apmass.c
 *
 * mkII aperture mass statistic 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kernels.h"
#include "../../utils/error.h"


#define usage "\n\
NAME\n\
	apmass - calculate mkII style aperture mass profile\n\
\n\
SYNOPSIS\n\
	apmass [options]\n\
		-c xc yc	# centre of coords (1024,1024)\n\
		-n NX NY	# box dimensions (2048, 2048)\n\
		-e eps		# softening parameter (0.01)\n\
		-l lossfac	# signal loss factor (1.0)\n\
\n\
DESCRIPTION\n\
	Aperture mass statistic: reads cat file; outputs info analogous to etprofile\n\
	but for rectangular aperture of linear size alpha times box size using mkII kernel.\n\
	Remember the centre of coords is the point at which lines from corners of box\n\
	through corners of aperture meet --- NOT the centroid of the aperture.\n\
	Calculates the mean kappa in aperture rel to mean kappa in surrounding strip\n\
	which depends only on shear estimates outside the aperture.\n\
	The mean surface density nbar is also calculated using only galaxies outside aperture.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n"

#define	Pi	M_PI
#define N	10

main (int argc, char *argv[])
{
	int	arg = 1;
	int	ix, iy, NX, NY, ngal[N], ialpha, nalpha;
	double	lossfac, K[2], xg, yg, x0, y0, gamma1, gamma2, eps, area, dalpha;
	double	alpha, alpha_ap[N], kappabar[N], kappaerr[N], ggsum[N], KKsum[N], gamma_rms[N];
	FILE	*inputpipe;
	char	line[1024];
 

	/* defaults */
	NX = NY = 2048;
	iy = ix = 1024;
	lossfac = 1.0;
	eps = 0.01;
	
	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			fprintf(stderr, usage);
			exit(1);
		}
		switch (argv[arg++][1]) {
			case 'n':
				sscanf(argv[arg++], "%d", &NX);
				sscanf(argv[arg++], "%d", &NY);
				break;
			case 'c':
				sscanf(argv[arg++], "%d", &ix);
				sscanf(argv[arg++], "%d", &iy);
				break;
			case 'l':
				sscanf(argv[arg++], "%lf", &lossfac);
				break;
			case 'e':
				sscanf(argv[arg++], "%lf", &eps);
				break;
			case 'u':
			default:
				fprintf(stderr, usage);
				exit(1);
				break;
		}
	}

	/* set up alpha_ap array */
	nalpha = 10;
	dalpha = 0.1;
	for (ialpha = 0; ialpha < nalpha; ialpha++)
		alpha_ap[ialpha] = dalpha * ialpha;
		
	x0 = (ix - NX / 2);
	y0 = (iy - NY / 2);
	area = (NX * NY);
	for (ialpha = 1; ialpha < nalpha; ialpha++) {
		kappabar[ialpha] = ggsum[ialpha] = KKsum[ialpha] = 0.0;
		ngal[ialpha] = 0;
	}

	inputpipe = popen("lc -o x e", "r");
	if (!inputpipe)
		error_exit("apmass: failed to open lc-pipe for input\n");

	while (fgets(line, 1024, inputpipe)) {
		sscanf(line, "%lf %lf %lf %lf", &xg, &yg, &gamma1, &gamma2);
		xg -= NX / 2;
		yg -= NY / 2;
		gamma1 /= lossfac;
		gamma2 /= lossfac;
		rectkernel(K, &alpha, x0, y0, xg, yg, 0.5 * NX, 0.5 * NY, eps);
		for (ialpha = 1; ialpha < nalpha; ialpha++) {
			if (alpha > alpha_ap[ialpha]) {
				ngal[ialpha]++;
				kappabar[ialpha] += (gamma1 * K[0] + gamma2 * K[1]);
				ggsum[ialpha] += gamma1 * gamma1 + gamma2 * gamma2;
				KKsum[ialpha] += K[0] * K[0] + K[1] * K[1];
			}
		}
	}

	/* now we have to normalise: remember rectkernel really returns (area * kernel) */
	/* so we just divide kappabar by ngal */
	for (ialpha = 1; ialpha < nalpha; ialpha++) {
		if (ngal[ialpha]) {
			kappabar[ialpha] /= ngal[ialpha];
			gamma_rms[ialpha] = sqrt(0.5 * ggsum[ialpha] / ngal[ialpha]);
			kappaerr[ialpha] = gamma_rms[ialpha] * sqrt(KKsum[ialpha]) / ngal[ialpha];
		} else {
			gamma_rms[ialpha] = kappaerr[ialpha] = 0.0;
		}
	}

        /* and output */
	fprintf(stdout, "#    alpha   kappabar      error       ngal\n");
	for (ialpha = 1; ialpha < nalpha; ialpha++) {
		fprintf(stdout, "%10.2lf %10.5lf %10.5lf %10d\n", 
			alpha_ap[ialpha], kappabar[ialpha], kappaerr[ialpha], ngal[ialpha]);
	}
	exit(0);
}



