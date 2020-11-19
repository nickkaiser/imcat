#define usage "\n\n\n\
NAME\n\
	fitstoascii --- convert a fits image to ascii format\n\
\n\
SYNOPSIS\n\
	fitstoascii [option...]\n\
		-f fmt	# use this format string (%%13.8g)\n\
		-l	# 1 col list output format\n\
		-t 	# x,y,pixval triplets\n\
		-g	# generate gnuplot format list\n\
\n\
DESCRIPTION\n\
	converts a fits image to ascii format\n\
	format must be valid for float (i.e. f,e or g)\n\
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

#define	DEFAULT_MODE	0
#define	ONE_COL_MODE	1
#define	TRIPLET_MODE	2
#define	GNUPLOT_MODE	3

int		main(int argc, char *argv[])	
{
	int		N1, N2, i, j, arg = 1, w, p;
	int		mode;
	fitsheader	*fits;
	char		format[32];
	float		**f;
	
	/* defaults */
	sprintf(format,"%%13.8g");
	mode = DEFAULT_MODE;

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-')
			error_exit(usage);
		switch (argv[arg++][1]) {
			case 'f':
				sscanf(argv[arg++], "%s", format);
				break;
			case 'l':
				mode = ONE_COL_MODE;
				break;
			case 't':
				mode = TRIPLET_MODE;
				break;
			case 'g':
				mode = GNUPLOT_MODE;
				break;
			default:
				error_exit(usage);
				break;
		}
	}
	
	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);
	fprintf(stdout, "# %3d %3d\n", N2, N1);
	for (i = 0; i < N2; i++) {
		for (j = 0; j < N1; j++) {
			if (mode == TRIPLET_MODE)
				fprintf(stdout, "%5d %5d ", j, i);
			fprintf(stdout, format, f[i][j]);
			if (mode != DEFAULT_MODE)
				fprintf(stdout, "\n");
		}
			if (mode == DEFAULT_MODE || mode == GNUPLOT_MODE)
				fprintf(stdout, "\n");
	}
	exit(0);
}




