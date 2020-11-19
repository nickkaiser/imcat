#define	usage "\n\n\n\
NAME\n\
	fitstocat --- make lc format list of fits image values\n\
\n\
SYNOPSIS\n\
	fitstocat [-u] [-b]\n\
\n\
DESCRIPTION\n\
	'fitstocat' reads a FITS image from stdin and sends a list of\n\
	pixel values to stdout as a lc format catalogue containing\n\
	pixel indices x[] and values f.\n\
\n\
	Options:\n\
		-u		# print this message\n\
		-b		# generate a binary format catalogue\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n\n\n"		


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../imlib/fits.h"
#include "../utils/error.h"
#include "../utils/arrays.h"
#include "../utils/args.h"

static int		binarycat, xdim, *idim;
static fitsheader	*fits;
static double		*x, f;
static float		*F;

void	writeimage(int level);

main(int argc, char *argv[])	
{
	char		*flag, sysstring[512], argstring[512];

	/* defaults */
	binarycat = 0;
	
	/* parse args */
	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'b':
				binarycat = 1;
				break;
			default:
				error_exit(usage);
		}
	}

	/* read the fits file header */
	fits = readfitsheader(stdin);

	/* allocate the pixel index */
	xdim = fits->ndim;
	x = (double *) calloc(xdim, sizeof(double));
	idim = (int *) calloc(xdim, sizeof(int));

	/* allocate the image line */
	F = (float *) calloc(fits->n[0], sizeof(float));

	/* output the cat header */
	argsToString(argc, argv, argstring);
	sprintf(sysstring, "lc -C -N '1 %d x' -n f -x -a '%s' < /dev/null", xdim, argstring);
	if (binarycat) {
		strcat(sysstring, " -b");
	}
	system(sysstring);

	/* write the image recursively */
	writeimage(fits->ndim);

	exit(0);
}


void	writeimage(int level)
{
	int	i, j;
	double	Fd;

	if (level > 1) {
		level--;
		for (idim[level] = 0; idim[level] < fits->n[level]; idim[level]++) {
			x[level] = (double) idim[level];
			writeimage(level);
		} 
	} else {
		readfitsline(F, fits);
		for (i = 0; i < fits->n[0]; i++) {
			x[0] = (double) i;
			if (binarycat) {
				fwrite(x, sizeof(double), xdim, stdout);
				Fd = (double) F[i];
				fwrite(&Fd, sizeof(double), 1, stdout);
			} else {
				fprintf(stdout, " ");
				for (j = 0; j < xdim; j++) {
					fprintf(stdout, " %14.8lg", x[j]);
				}
				fprintf(stdout, " %14.8g\n", F[i]);
			}
		}
	}
}

