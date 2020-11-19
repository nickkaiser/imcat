/*
 * hmstodec.c
 */

#define usage "\n\
NAME\n\
	hmstodecimal --- convert an angle from h:m:s to decimal degrees format\n\
\n\
SYNOPSIS\n\
	hmstodec [angle]\n\
\n\
DESCRIPTION\n\
	hmstodecimal takes a string as argument; decipers it\n\
	using \"%%d:%%d:%%lf\" format specification to h, m, s;\n\
	checks that m, s lie in range 0-60 and writes the \n\
	angle\n\
		theta = (sign(h)) * (15 * fabs(h) + m / 60.0 + s / 3600.0)\n\
	to standard output.\n\
\n\
	If no arguments are supplied, it reads a from stdin.  Input\n\
	stream should contain one angle per line.\n\
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
	char	line[1024];

	if (argc != 2) {
		if (argc == 1) {
			while (fgets(line, 1024, stdin)) {
				if (!xmstodecimal(line, &angle)) { 
                			error_exit("failed to convert input");	
				}
				fprintf(stdout, "%13.8lf\n", 15 * angle);
			}
			exit(0);
		} else {
			error_exit(usage);
		}
	}
	if (argv[1][0] == '-' && argv[1][1] == 'u') {
		error_exit(usage);
	}
	if (!xmstodecimal(argv[1], &angle)) {
		error_exit(usage);
	}
	fprintf(stdout, "%13.8lf\n", 15 * angle);
	exit(0);
}
