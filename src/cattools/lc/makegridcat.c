/*
 * makerandcat
 */


#include <stdio.h>
#include "../../imlib/fits.h"
#include "../../utils/args.h"
#include "../../catlib/cat.h"

#define usage "\n\n\
NAME\n\
        makegridcat - generate catalogue containing a grid of points\n\
\n\
SYNOPSIS\n\
        makegridcat [ -u | -b ] nx [ny [nz ....]]\n\
\n\
DESCRIPTION\n\
	Makegridcat generates a lc-format catalogue containing a unit spaced\n\
	grid of points for testing purposes.\n\
\n\
	If invoked with a single argument nx, makegridcat generates a catalogue\n\
	containing nx objects consisting of a numerical scalar x with values\n\
	0,1...(nx-1).\n\
	If invoked with two arguments nx, ny, it generates a catalogue\n\
	containing nx * ny objects consisting of a 2-vector x[]\n\
	where x[0] runs from 0 through nx-1 and x[1] runs from 0 through ny-1.\n\
	If invoked with three arguments nx, ny, nz, it generates a catalogue\n\
	containing nx * ny * nz objects consisting of a 3-vector x[]\n\
	where x[0] runs from 0 through nx-1 and x[1] runs from 0 through ny-1\n\
	and x[2] runs from 0 through nz-1, and so on.\n\
\n\
	Default is text format output.  Use -b option to generate a binary\n\
	format catalogue.\n\
\n\
AUTHOR\n\
        Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n\n"

static int	ndim, *dim, binarycat;
static double	*x;

void	writegridcat(int level);

main(int argc, char *argv[]){
	int	i, argtype;
	char	lcstring[128], argstring[128], *flag;

	/* defaults */
	binarycat = 0;
	ndim = argc - 1;

	/* parse args */
	argsinit(argc, argv, usage);
	argtype = nextargtype();
	if (argtype == NO_ARG) {
		error_exit(usage);
	}	
	if (argtype == FLAG_ARG) {
		flag = getflag();
		if (flag[0] == 'b') {
			binarycat = 1;
			ndim--;
		} else {
			error_exit(usage);
		}
	} 

	/* set up the dim[] array */
	dim = (int *) calloc(ndim, sizeof(int));
	for (i = 0; i < ndim; i++) {
		dim[i] = getargi();
	}

	/* set up the x[] array */
	x = (double *) calloc(ndim, sizeof(double));

	/* write the header */
	argsToString(argc, argv, argstring);
        sprintf(lcstring, "lc -C -x -a \"history: %s\" -N '1 %d x' < /dev/null", argstring, ndim);
	if (binarycat) {
		strcat(lcstring, " -b");
	}
	system(lcstring);

	/* now we recursively generate the output */
	writegridcat(ndim);	

	exit(0);
}

void	writegridcat(int level)
{
	int	i, lev;

	level--;
	for (i = 0; i < dim[level]; i++) {
		x[level] = (double) i;
		if (level) {
			writegridcat(level);
		} else {
			if (binarycat) {
				fwrite(x, ndim, sizeof(double), stdout);
			} else {
				fprintf(stdout, " ");
				for (lev = 0; lev < ndim; lev++) {
					fprintf(stdout, LC_NUM_FMT, x[lev]);
				}
				fprintf(stdout, "\n");
			}
		}
	}
}	


