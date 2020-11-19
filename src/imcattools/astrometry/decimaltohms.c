/*
 * dectohms.c
 */

#define usage "\n\
NAME\n\
	decimaltomhs --- convert an angle from decimal degrees to h:m:s format\n\
\n\
SYNOPSIS\n\
	decimaltohms angle\n\
\n\
DESCRIPTION\n\
	decimaltohms takes a decimal format angle as argument and\n\
	generates a string corresponding to that argument in h:m:s\n\
	format\n\
\n\
AUTHOR\n\
	Nick Kaiser	kaiser@hawaii.edu\n\
\n\n"


#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#include "utils/error.h"
#include "radecio.h"

main (int argc, char *argv[])
{
	double	angle;
	char	hmsstring[64];

	if (argc != 2) {
		error_exit(usage);
	}
	if (argv[1][0] == '-' && argv[1][1] == 'u') {
		error_exit(usage);
	}
	if (1 != sscanf(argv[1], "%lf", &angle)) {
		error_exit(usage);
	}
	decimaltoxms(angle / 15.0, hmsstring);
	fprintf(stdout, "%s\n", hmsstring);
	exit(0);
}
