#define usage "\n\n\n\
NAME\n\
	getplane --- extracts an image plane from multi-dimensional fits image\n\
\n\
SYNOPSIS\n\
	getplane p\n\
\n\
DESCRIPTION\n\
	\"getplane\" reads from stdin a N-dimensional fits image\n\
		f_in[Nn][Nn-1]....[N2][N1]\n\
	and writes to stdout a N-1 dimensional image\n\
		f_out[Nn-1]....[N2][N1]\n\
	with values\n\
		f_out[xn-1]...[x2][x1] = f_in[p][xn-1]...[x2][x1]\n\
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

int		main(int argc, char *argv[])	
{
	int		p, dim, nlines, nskip;
	fitsheader	*fits;
	void		*f;

	/* parse args */
	if (argc != 2) {
		error_exit(usage);
	}
	if (argv[1][0] == '-') {
		error_exit(usage);
	}
	if (1 != sscanf(argv[1], "%d", &p)) {
		error_exit(usage);
	}

	/* read and write headers */
	fits = readfitsheader(stdin);
	if (p < 0 || p >= fits->n[fits->ndim-1]) {
		error_exit("getplane: p out of allowed range\n");
	}
	add_comment(argc, argv, fits);
	fits->ndim--;
	fits->intpixtype = fits->extpixtype;
	writefitsheader(fits);


	/* figure number of lines to read and how many to skip*/
	nlines = 1;
	for (dim = 1; dim < fits->ndim; dim++) {
		nlines *= fits->n[dim];
	}
	nskip = p * nlines;

	/* allocate line buffer */
	f = calloc(fits->n[0] * pixsize(fits->extpixtype), sizeof(char));
	skiplines(fits, nskip);
	while (nlines--) {
		readfitsline(f, fits);
		writefitsline(f, fits);
	} 
	writefitstail(fits);
	free(f);

	exit(0);
}
