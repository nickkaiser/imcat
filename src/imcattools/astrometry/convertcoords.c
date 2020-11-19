/*
 * convertcoords.c
 */

#define usage "\n\
NAME\n\
	convertcoords --- transform between different sky coordinate systems\n\
\n\
SYNOPSIS\n\
	convertcoords convtype angle1 angle2\n\
\n\
DESCRIPTION\n\
	convertcoords invokes Doug Mink's utilities for conversion\n\
	between different sky coordinate systems.\n\
\n\
	Angles can be expressed in either decimal or [hd]:m:s format.\n\
\n\
	The string 'convtype' specifies the type of conversion and can be\n\
	one of\n\
		fk524 		# Convert J2000(FK5) to B1950(FK4) coordinates\n\
		fk425		# Convert B1950(FK4) to J2000(FK5) coordinates\n\
		fk42gal		# Convert B1950(FK4) to galactic coordinates\n\
		fk52gal		# Convert J2000(FK5) to galactic coordinates\n\
		gal2fk4		# Convert galactic coordinates to B1950(FK4)\n\
		gal2fk5		# Convert galactic coordinates to J2000<FK5)   \n\
\n\
	Output angles are expressed in decimal degrees.\n\
AUTHOR\n\
	Nick Kaiser	kaiser@hawaii.edu\n\
\n\n"


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


#include "utils/error.h"
#include "radecio.h"

main (int argc, char *argv[])
{
	double		ang1, ang2;
	int		ang1type, ang2type;
	char		*convtype;


	if (argc != 4) {
		error_exit(usage);
	}
	if (!strcmp("-u", argv[1])) {
		error_exit(usage);
	}
	convtype = argv[1];
	ang1 = getangle(argv[2], &ang1type);	
	if (ang1type == HMS_ANGLE_TYPE) {
		ang1 *= 15.0;
	}
	ang2 = getangle(argv[3], &ang2type);	


	if (!strcmp(convtype, "fk524")) {
		fk524 (&ang1, &ang2);
	}
	if (!strcmp(convtype, "fk425")) {
		fk425 (&ang1, &ang2);
	}
	if (!strcmp(convtype, "fk42gal")) {
		fk42gal (&ang1, &ang2);
	}
	if (!strcmp(convtype, "fk52gal")) {
		fk52gal (&ang1, &ang2);
	}
	if (!strcmp(convtype, "gal2fk4")) {
		gal2fk4 (&ang1, &ang2);
	}
	if (!strcmp(convtype, "gal2fk5")) {
		gal2fk5 (&ang1, &ang2);
	}
	fprintf(stdout, "%13.8lf %13.8lf\n", ang1, ang2);
}
