/*
 * decimaltodms.c
 */

#define usage "\n\
NAME\n\
	decimaltodms --- convert an angle from decimal degrees to d:m:s format\n\
\n\
SYNOPSIS\n\
	decimaltodms angle\n\
\n\
DESCRIPTION\n\
	decimaltodms takes a decimal format angle (in degrees) as argument and\n\
	generates a string corresponding to that argument in d:m:s\n\
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
	char	dmsstring[64];

	if (argc != 2) {
		error_exit(usage);
	}
	if (argv[1][0] == '-' && argv[1][1] == 'u') {
		error_exit(usage);
	}
	if (1 != sscanf(argv[1], "%lf", &angle)) {
		error_exit(usage);
	}
	decimaltoxms(angle , dmsstring);
	fprintf(stdout, "%s\n", dmsstring);
	exit(0);
}
