#define usage "\n\n\n\
NAME\n\
	makesubimage --- extracts a subimage from a fits file\n\
\n\
SYNOPSIS\n\
	makesubimage x y dx dy [options...]\n\
		-o	# set unmapped pixels to zero\n\
		-n N3max	# output only the first N3max planes if processing 3D image\n\
\n\
DESCRIPTION\n\
	'makesubimage' extracts a dx x dy subimage with\n\
	origin at (x,y) from the input image.  By default,\n\
	if the subimage extends beyond the source image\n\
	then extra pixels are set to MAGIC.  Use -o\n\
	option to force these to zero.\n\
	Reads from stdin, writes to stdout.\n\
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

#define MAGIC FLOAT_MAGIC

int		main(int argc, char *argv[])	
{
	int		arg = 1, N1, N2, N3 = 1, M1, M2, i0, j0, i, j, jj, N3max = 0;
	fitsheader	*fitsin , *fitsout;
	float		*fin, *fout, defaultvalue;

	/* defaults */
	defaultvalue = MAGIC;

	if (argc < 5)
		error_exit(usage);

	if (1 != sscanf(argv[arg++], "%d", &j0) ||
		1 != sscanf(argv[arg++], "%d", &i0) || 
		1 != sscanf(argv[arg++], "%d", &M1) || 
		1 != sscanf(argv[arg++], "%d", &M2)) error_exit(usage);

	while (arg < argc) {
		if (argv[arg][0] != '-')
			error_exit(usage);
		switch (argv[arg++][1]) {
			case 'o':
				defaultvalue = 0.0;
				break;
			case 'n':
				if (1 != sscanf(argv[arg++], "%d", &N3max)) error_exit(usage);
				break;
			default:
				error_exit(usage);
		}
	}


	fitsin = readfitsheader(stdin);
	N1 = fitsin->n[0];
	N2 = fitsin->n[1];
	if (fitsin->ndim == 3) {
		if (N3max && (fitsin->n[2] > N3max)) {
			fitsin->n[2] = N3max;
		}
		N3 = fitsin->n[2];
	}
	fin = (float *) calloc(N1, sizeof(float));
	fout = (float *) calloc(M1, sizeof(float));

	fitsout = copyfitsheader(fitsin);
	add_comment(argc, argv, fitsout);
	fitsout->n[0] = M1;
	fitsout->n[1] = M2;
	writefitsheader(fitsout);

	while (N3 > 0) {
		N3--;
		if (i0 > 0) {
			skiplines(fitsin, i0);
		}
		for (i = i0; i < i0 + M2; i++) {
			if (i >= N2 || i < 0)
				for (j = 0; j < N1; j++)
					fin[j] = defaultvalue;
			else
				readfitsline(fin, fitsin);
			for (jj = 0; jj < M1; jj++)
				fout[jj] = (jj + j0 >= N1 ? defaultvalue : (jj + j0 < 0 ? defaultvalue : fin[jj + j0]));
			writefitsline(fout, fitsout);
		}
		if ((N2 - i) > 0) {
			skiplines(fitsin, N2 - i);
		}
	}
	writefitstail(fitsout);
	exit(0);
}


