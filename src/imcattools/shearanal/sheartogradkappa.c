/*
 * sheartogradkappa.c
 */

#define usage "\n\
NAME\n\
	sheartogradkappa --- calculate gradkappa image from shear image\n\
\n\
SYNOPSIS\n\
	sheartogradkappa\n\
\n\
DESCRIPTION\n\
	'sheartogradkappa' reads shear from a fits file on stdin and calculates\n\
	grad - kappa which is sent to stdout\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n"

#include <stdio.h>
#include <math.h>
#include "../../imlib/fits.h"
#include "../../utils/arrays.h"

main (int argc, char *argv[])
{
	int		arg = 1, N1, N2, i, j, pixtype;
	float		**gamma1, **gamma2, **gradkappa1, **gradkappa2;
	float		gamma11, gamma12, gamma21, gamma22;
	fitsheader	*fitsin, *fitsout;

	/* defaults */

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			fprintf(stderr, "%s", usage);
			exit(-1);
		}
		switch (argv[arg++][1]) {
			default:
				fprintf(stderr, "%s", usage);
				exit(-1);
				break;
		}
	}

	fitsin = readfitsheader(stdin);
	N1 = fitsin->n[0];
	N2 = fitsin->n[1];
	allocFloatArray(&gamma1, N1, 2 * N2);
	for (i = 0; i < 2 * N2; i++) {
		readfitsline(gamma1[i], fitsin);
	}
	gamma2 = gamma1 + N2;
/*
	N1--;
	N2--;
*/
	fitsout = copyfitsheader(fitsin);
	allocFloatArray(&gradkappa1, N1, 2 * N2);
	gradkappa2 = gradkappa1 + N2;

	/* calculate gradient */
	for (i = 0; i < N2 - 1; i++) {
		for (j = 0; j < N1 - 1; j++) {
			gamma11 = gamma1[i][j+1] - gamma1[i][j];
			gamma21 = gamma2[i][j+1] - gamma2[i][j];
			gamma12 = gamma1[i+1][j] - gamma1[i][j];
			gamma22 = gamma2[i+1][j] - gamma2[i][j];
			gradkappa1[i][j] = gamma11 + gamma22;
			gradkappa2[i][j] = gamma21 - gamma12;
		}
	}
	add_comment(argc, argv, fitsout);
	writefitsheader(fitsout);
	for (i = 0; i < 2 * N2; i++) {
		writefitsline(gradkappa1[i], fitsout);
	}
	writefitstail(fitsout);
	exit(0);
}





