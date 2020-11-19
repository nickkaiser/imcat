#define usage "\n\n\n\
NAME\n\
	make2Dotf --- make 2-D fast guiding otf with offset guide star\n\
SYNOPSIS\n\
	make2Dotf [-N N] [-R r_outer] [-r r0] [-D D] [-o opfile] [-z zmax] [-i imname]\n\
\n\
DESCRIPTION\n\
	'make2Dotf' computes the OTF gk(z) for perfect fast guiding.\n\
\n\
	Options are\n\
		-N N		# image size in pixels (512)\n\
		-r r0		# Fried length in m (0.2)\n\
		-D D		# telescope diameter (1.0)\n\
		-z zmax		# upper limit for integerized z (N/2)\n\
		-d dz		# step in integer z (2)\n\
		-p nphi		# number of rays in azimuthal angle (5)\n\
		-g rg		# distance to guide star (1.0)\n\
\n\
	Output is a fits image of the OTF.\n\
\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "imlib/fits.h"
#include "utils/error.h"
#include "utils/args.h"

int		main(int argc, char *argv[])	
{
	char		*flag, lcstring[512];
	int		N, dz, zmax, iphi, nphi, z, nz, ix, iy, ik1, ik2, iphi1, iphi2;
	double		**gkpolar;
	float		**gk;
	double		kx, ky, k, rg, xg, yg, phi, dphi, r0, D, deltaphi, deltak, gk1, gk2;
	fitsheader	*fits;
	FILE		*lcpipe;

	/* defaults */
	N 	= 512;
	r0	= 0.2;
	D	= 1.0;
	zmax	= 0;
	dz	= 2;
	nphi	= 5;
	rg 	= 1.0;

	/* parse args */
	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'N':
				N = getargi();
				break;
			case 'r':
				r0 = getargd();
				break;
			case 'D':
				D = getargd();
				break;
			case 'z':
				zmax = getargi();
				break;
			case 'g':
				rg = getargd();
				break;
			case 'd':
				dz = getargi();
				break;
			case 'p':
				nphi = getargi();
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}

	if (!zmax) {
		zmax = N / 2;
	}

	/* compute nz, and allocate arrays */
	nz = 0;
	for (z = 0; z < zmax; z += dz) {
		nz++;
	}
	gkpolar = (double **) calloc(nphi, sizeof(double *));
	for (iphi = 0; iphi < nphi; iphi++) {
		gkpolar[iphi] = (double *) calloc(nz, sizeof(double));
	}
	allocFloatArray(&gk, N, N);

	/* call makeotf repeatedly to generate gkpolar */
	dphi = M_PI / (nphi - 1);
	for (iphi = 0; iphi < nphi; iphi++) {
		phi = iphi * dphi;
		xg = rg * cos(phi);
		yg = rg * sin(phi);
		sprintf(lcstring, "makeotf -N %d -r %lg -D %lg -z %d -d %d -g %lg %lg | lc -b -o gktilt", 
			N, r0, D, zmax, dz, xg, yg);
		lcpipe = popen(lcstring, "r");
		fprintf(stderr, "executing: %s\n", lcstring);
		fread(gkpolar[iphi], sizeof(double), nz, lcpipe);
		close(lcpipe);
	}

	/* interpolate gkpolar onto gk */
	for (iy = 0; iy < N; iy++) {
		ky = (iy - 0.5 * N) / dz;
		for (ix = 0; ix < N; ix++) {
			kx = (ix - 0.5 * N) / dz;
			k = sqrt(kx * kx + ky * ky);
			ik1 = (int) floor(k);
			deltak = k - ik1;
			ik2 = ik1 + 1;
			phi = atan2(fabs(ky), kx);
			iphi1 = (int) floor(phi / dphi);
			if (iphi1 == nphi - 1) {
				iphi1--;
			}
			deltaphi = phi / dphi - iphi1;
			iphi2 = iphi1 + 1;
			if (ik2 < nz) {
				gk1 = (1 - deltak) * gkpolar[iphi1][ik1] + deltak * gkpolar[iphi1][ik2];
				gk2 = (1 - deltak) * gkpolar[iphi2][ik1] + deltak * gkpolar[iphi2][ik2];
				gk[iy][ix] = (1 - deltaphi) * gk1 + deltaphi * gk2;
			}
		}
	}
	write2Dfloatimage(gk, new2Dfitsheader(N, N, FLOAT_PIXTYPE));
	exit(0);
}


