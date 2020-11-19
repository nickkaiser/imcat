/*
 * hmstodec.c
 */

#define usage "\n\
NAME\n\
	dmstodecimal --- convert an angle from d:m:s to decimal degrees format\n\
\n\
SYNOPSIS\n\
	dmstodecimal angle\n\
\n\
DESCRIPTION\n\
	dmstodecimal takes a string as argument; decipers it\n\
	using \"%%d:%%d:%%lf\" format specification to d, m, s;\n\
	checks that m, s lie in range 0-60 and writes the \n\
	angle\n\
		theta = (sign(d)) * (fabs(d) + m / 60.0 + s / 3600.0)\n\
	to standard output.\n\
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

	if (argc != 2) {
		error_exit(usage);
	}
	if (argv[1][0] == '-' && argv[1][1] == 'u') {
		error_exit(usage);
	}
	if (!xmstodecimal(argv[1], &angle)) {
		error_exit(usage);
	}
	fprintf(stdout, "%13.8lf\n", angle);
	exit(0);
}
