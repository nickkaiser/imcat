#define usage "\n\n\n\
NAME\n\
	colorbox - generate example color image\n\
\n\
SYNOPSIS\n\
	colourbox N\n\
\n\
DESCRIPTION\n\
	\"colourbox\" generates a N x N 2-colour image  which spans the spectrum\n\
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

static	short 	*f;
static	int	ncolours;

int		main(int argc, char *argv[])	
{
	int	Naxes, N, colour;
	char	argstring[COM_LENGTH], *color[1];
	fitsheader	*fits;
	fitscomment	*com;
	int	i, j;
	int	NAXIS[3];
	
	if (argc != 2 || argv[1][0] == '-')
		error_exit(usage);
	
	NAXIS[2] = ncolours = 2;
	sscanf(argv[1], "%d", &N);
	NAXIS[0] = NAXIS[1] = N;
        Naxes = 3;

	fits = (fitsheader *) calloc(1, sizeof(fitsheader));
	fits->extpixtype = SHORT_PIXTYPE;
	fits->intpixtype = SHORT_PIXTYPE;
	fits->ndim = 3;
	fits->n[0] = NAXIS[0];
	fits->n[1] = NAXIS[1];
	fits->n[2] = NAXIS[2];
	fits->bscaling = 0;
	fits->linebuffer = 0;
	fits->opstream = stdout;

	add_comment(argc, argv, fits);
	com = newtextcomment("HISTORY", "color: red green", NULL);
	appendcomment(com, fits);

	writefitsheader(fits);
	f = (short *) calloc(N, sizeof(short));
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++)
			f[j] = (short) j;
		writefitsline(f, fits);
	}
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++)
			f[j] = (short) i;
		writefitsline(f, fits);
	}
	writefitstail(fits);
	exit (0);
}

