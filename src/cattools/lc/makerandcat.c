/*
 * makerandcat
 */


#include <stdio.h>
#include <stdlib.h>
#include "../../utils/error.h"
#include "../../utils/args.h"
#include "../../catlib/cat.h"

#define usage "\n\n\
NAME\n\
        makerandcat - generate catalogue containing random values\n\
\n\
SYNOPSIS\n\
        makerandcat n [-seed seed] [-dim dim] [-b]\n\
\n\
DESCRIPTION\n\
	Makerandcat generates a lc-format catalogue containing random values\n\
	for testing purposes.  By default the catalogue contains n objects consisting\n\
	of a 2-vector x[2], with values x[0], x[1] uniformly distributed\n\
	on the interval 0,1 and is in text format.\n\
\n\
	Options:\n\
		-seed 	seed	# supply an integer seed for random number generator (1)\n\
		-dim	dim	# specify alternative dimension for x-vector (2)\n\
		-b		# generate a binary format catalogue.\n\
\n\
AUTHOR\n\
        Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n\n"

double	drand48();

main(int argc, char *argv[]){
	int	seed, n, dim, obj, i, binarycat;
	double	*x;
	char	*flag, lcstring[128], argstring[512];

	/* defaults */
	dim = 2;
	seed = 1;
	binarycat = 0;

	/* parge args */
	argsinit(argc, argv, usage);
	if (nextargtype() == FLAG_ARG) {
		error_exit(usage);
	}
	n = getargi();
	while (flag = getflag()) {
		if (!strcmp(flag, "seed")) {
			seed = getargi();
		} else if (!strcmp(flag, "dim")) {
			dim = getargi();
		} else if (!strcmp(flag, "b")) {
			binarycat = 1;
		} else {
			error_exit(usage);
		}
	}

	/* seed the random number generator */
	srand48((long int) seed);

	/* allocate the data */
	x = (double *) calloc(dim, sizeof(double));

	argsToString(argc, argv, argstring);
	sprintf(lcstring, "lc -C -x -a \"history: %s\" -N '1 %d x' < /dev/null", argstring, dim);
	if (binarycat) {
		strcat(lcstring, " -b");
	}
	system(lcstring);
		
	for (obj = 0; obj < n; obj++) {
		for (i = 0; i < dim; i++) {
			x[i] = drand48();
		}
		if (binarycat) {
			fwrite(x, sizeof(double), dim, stdout);
		} else {
			fprintf(stdout, " ");
			for (i = 0; i < dim; i++) {
				fprintf(stdout, LC_NUM_FMT, x[i]);
			}
			fprintf(stdout, "\n");
		}
	}

	exit(0);
}
