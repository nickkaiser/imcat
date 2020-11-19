#define usage "\n\n\n\
NAME\n\
	unscrunch - expand a FITS image by a factor 2\n\
\n\
SYNOPSIS\n\
	unscrunch\n\
\n\
DESCRIPTION\n\
	expands a fits image by a factor 2.\n\
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

int		main(int argc, char *argv[])	
{
	int		i, j, N1in, N2in, N1out, N2out;
	fitsheader	*fitsin, *fitsout;
	float		*fin, *fout;
	
	if (argc > 1)
		error_exit(usage);	
	fitsin = readfitsheader(stdin);
	N1in = fitsin->n[0];
	N2in = fitsin->n[1];
	N1out = N1in * 2;
	N2out = N2in * 2;
	fitsout = copyfitsheader(fitsin);
	fitsout->n[0] = N1out;
	fitsout->n[1] = N2out;
	add_comment(argc, argv, fitsout);
	writefitsheader(fitsout);
	fin = (float *) calloc(N1in, sizeof(float));
	fout = (float *) calloc(N1out, sizeof(float));
	for (i = 0; i < N2in; i++) {
		readfitsline(fin, fitsin);
		for (j = 0; j < N1in; j++) {
			fout[2 * j] = fout[2 * j + 1] = fin[j];
		}
		writefitsline(fout, fitsout);
		writefitsline(fout, fitsout);
	}
	writefitstail(fitsout);
	exit(0);
}


