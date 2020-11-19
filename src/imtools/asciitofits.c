#define usage "\n\n\n\
NAME\n\
	asciitofits --- converts an ascii image to fits format\n\
\n\
SYNOPSIS\n\
	asciitofits [option...]\n\
\n\
DESCRIPTION\n\
	'asciitofits' reads an ascii format image from stdin\n\
	and writes a fits image to stdout\n\
		-f	# create float format fits file\n\
		-i	# create 32 bit int format fits file\n\
		-p	# data in 2-hex digits a la postscript   \n\
	ascii format is # N2 N1 then N2 x N1 pixvals \n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../imlib/fits.h"
#include "../utils/error.h"
#include "../utils/arrays.h"

int		main(int argc, char *argv[])	
{
	int		N1, N2, i, j, arg = 1, readpostscript, ok, itemp, pixtype;
	fitsheader	*fits;
	float		**f, temp;
	
	/* defaults */
	readpostscript = 0;
	pixtype = FLOAT_PIXTYPE;

	while (arg < argc) {
		if (argv[arg][0] != '-')
			error_exit(usage);
		switch (argv[arg++][1]) {
			case 'f':
				pixtype = FLOAT_PIXTYPE;
				break;
			case 'i':
				pixtype = INT_PIXTYPE;
				break;
			case 'p':
				readpostscript = 1;
				break;
			default:
				error_exit(usage);
		}
	}
	
	if (2 != fscanf(stdin, "# %d %d\n", &N2, &N1))
		error_exit("asciitofits: ascii file format problem\n");
/*	fprintf(stderr, "%d x %d image\n", N1, N2);*/
	allocFloatArray(&f, N1, N2);
	for (i = 0; i < N2; i++) {
		for (j = 0; j < N1; j++) {
			if (readpostscript) {
				ok = fscanf(stdin, "%2x", &itemp);
				temp = (float) itemp;
			} else {
				ok = fscanf(stdin, "%f", &temp);
			}
			if (!ok)
				error_exit("asciitofits: ascii file read problem\n");
			f[i][j] = temp;
		}
	}
	fits = new2Dfitsheader(N1, N2, pixtype);
	add_comment(argc, argv, fits);
	write2Dfloatimage(f, fits);
	exit(0);
}


