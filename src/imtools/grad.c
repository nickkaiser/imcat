/*
 * grad.c
 */

#define usage "\n\
NAME\n\
	grad --- take the gradient of a fits image\n\
\n\
SYNOPSIS\n\
	grad [-p]\n\
\n\
DESCRIPTION\n\
	'grad' reads a N1 by N2 2-D fits file f[][] from stdin and calculates discrete\n\
	second derivatives.  The output is sent to stdout in the form of a \n\
	N1 by N2 by 2 image where the zeroth plane is\n\
		df[0][y][x] = f[y][x+1] - f[y][x]\n\
	and the second plane is\n\
		df[1][y][x] = f[y+1][x] - f[y][x]\n\
	By default the rightmost column of the first half and the last\n\
	row are zero. Use -p option to use periodic bc's.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n"

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "../imlib/fits.h"
#include "../utils/arrays.h"

#define MAGIC FLOAT_MAGIC

main (int argc, char *argv[])
{
	int		arg = 1, N1, N2, i, j, pixtype, periodic;
	float		**fx, **fy, **f;
	fitsheader	*fits;

	/* defaults */
	periodic = 0;

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			fprintf(stderr, "%s", usage);
			exit(-1);
		}
		switch (argv[arg++][1]) {
			case 'p':
				periodic = 1;
				break;
			default:
				fprintf(stderr, "%s", usage);
				exit(-1);
				break;
		}
	}

	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);

	allocFloatArray(&fx, N1, 2 * N2);
	fy = fx + N2;

	/* take the difference */
	for (i = 0; i < N2; i++) {
		for (j = 0; j < N1; j++) {
			if (f[i][j] == MAGIC) {
				fx[i][j] = fy[i][j] = MAGIC;
			} else {
				if (j < N1 - 1) {
					if (f[i][j+1] == MAGIC) {
						fx[i][j] = MAGIC;
					} else {
						fx[i][j] = f[i][j+1] - f[i][j];
					}
				} else {
					if (periodic) {
						if (f[i][0] == MAGIC) {
							fx[i][j] = MAGIC;
						} else {
							fx[i][j] = f[i][0] - f[i][j];
						}
					}
				}
				if (i < N2 - 1) {
					if (f[i+1][j] == MAGIC) {
						fy[i][j] = MAGIC;
					} else {
						fy[i][j] = f[i+1][j] - f[i][j];
					}
				} else {
					if (periodic) {
						if (f[0][j] == MAGIC) {
							fy[i][j] = MAGIC;
						} else {
							fy[i][j] = f[0][j] - f[i][j];
						}
					}
				}
			}
		}
	}
	add_comment(argc, argv, fits);
	fits->ndim = 3;
	fits->n[2] = 2;
	writefitsheader(fits); 
	for (i = 0; i < 2 * N2; i++) {
		writefitsline(fx[i], fits);
	}
	writefitstail(fits); 
	exit(0);
}





