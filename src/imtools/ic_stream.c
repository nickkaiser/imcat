#define usage "\n\n\n\
NAME\n\
	ic_stream --- apply simple arithmetic to a stream of images\n\
\n\
SYNOPSIS\n\
	ic_stream [plus | minus | times | divide] operandfits\n\
\n\
DESCRIPTION\n\
	'ic_stream plus myfits' reads a 2D image myfits and then\n\
	reads planes from a higher dimensional image, adds the value\n\
	from myfits and outputs the result.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


#include "../imlib/fits.h"
#include "../utils/error.h"

#define	ADD_MODE	0
#define SUB_MODE	1
#define	MUL_MODE	2
#define	DIV_MODE	3

int		main(int argc, char *argv[])	
{
	int		mode, N1, N2, nplanes, i, x, y;
	FILE		*operandf;
	fitsheader	*operandfits, *fits;
	float		**f, **foperand;

	/* parse args */
	if (argc != 3) {
		error_exit(usage);
	}
	switch (argv[1][0]) {
		case 'p':
			mode = ADD_MODE;
			break;
		case 'm':
			mode = SUB_MODE;
			break;
		case 't':
			mode = MUL_MODE;
			break;
		case 'd':
			mode = DIV_MODE;
			break;
		default:
			error_exit(usage);
	}

	/* open the operand image */
	operandf = fopen(argv[2], "r");
	if (!operandf) {
		error_exit("ic_stream : failed to open operand image\n");
	}
				
	/* read the operand header */
	operandfits = readfitsheader(operandf);
	if (operandfits->ndim != 2) {
		error_exit("ic_stream: operand must be a 2 dimensional fits image\n");
	}
	N1 = operandfits->n[0];
	N2 = operandfits->n[1];

	/* read the source fits header */
	fits = readfitsheader(stdin);
	if (fits->ndim < 2) {
		error_exit("ic_stream: operand must be a 2 or larger dimensional fits image\n");
	}
	if ((fits->n[0] != N1) || (fits->n[1] != N2)) {
		error_exit("ic_stream: image plane size mismatch\n");
	}
	nplanes = 1;
	for (i = 2; i < fits->ndim; i++) {
		nplanes *= fits->n[i];
	} 

	/* allocate and read the operand image */
	allocFloatArray(&foperand, N1, N2);
	readfitsplane((void *) foperand, operandfits);
	fclose(operandf);

	/* allocate the input image */
	allocFloatArray(&f, N1, N2);

	/* add history and output header */
	add_comment(argc, argv, fits);
	writefitsheader(fits);

	for (i = 0; i < nplanes; i++) {
		readfitsplane((void *) f, fits);
		switch (mode) {
			case ADD_MODE:
				for (y = 0; y < N2; y++) {
					for (x = 0; x < N1; x++) {
						f[y][x] += foperand[y][x];
					}
				}
				break;			
			case SUB_MODE:
				for (y = 0; y < N2; y++) {
					for (x = 0; x < N1; x++) {
						f[y][x] -= foperand[y][x];
					}
				}
				break;			
			case MUL_MODE:
				for (y = 0; y < N2; y++) {
					for (x = 0; x < N1; x++) {
						f[y][x] *= foperand[y][x];
					}
				}
				break;			
			case DIV_MODE:
				for (y = 0; y < N2; y++) {
					for (x = 0; x < N1; x++) {
						f[y][x] /= foperand[y][x];
					}
				}
				break;
			default:
				error_exit("ic_stream : internal error : bad mode value\n");		
		}
		writefitsplane((void *) f, fits);
	}
	return (0);
}


