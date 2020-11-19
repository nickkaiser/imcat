#define usage "\n\n\n\
NAME\n\
	cycleimage - translate a FITS image\n\
\n\
SYNOPSIS\n\
	cycleimage dx dy\n\
\n\
DESCRIPTION\n\
	\"cycleimage\" cycles an image moving each point\n\
	a distance dx, dy (in pixels).\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


#include "../fftlib/myfft.h"
#include "../imlib/fits.h"
#include "../utils/error.h"

int		main(int argc, char *argv[])	
{
	int		N1, N2, dx, dy, pixtype;
	fitsheader	*fits;
	float		**f;

	if (argc != 3)
		error_exit(usage);
	if ((1 != sscanf(argv[1], "%d", &dx)) || (1 != sscanf(argv[2], "%d", &dy)))
		error_exit(usage);
	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);
	cycleimage(f, N1, N2, dx, dy);
	add_comment(argc, argv, fits);
	write2Dfloatimage(f, fits);
	exit(0);
}



