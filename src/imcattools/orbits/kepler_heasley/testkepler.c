/*
 * testkepler - test kepler.c transformations

 * Nick Kaiser 02/09/03
 */

#include <stdio.h>
#include <math.h>
#include "kepler.h"

#define usage "\nNAME\n\
	testkepler\n\
\n\
SYNOPSIS\n\
	testkepler \n\
\n\
DESCRIPTION\n\
	testkepler reads r[3], v[3] values from stdin in lc-format\n\
	and converts these to Kepler elements and back again.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"


main (int argc, char *argv[])
{
	double	*buff, *r, *v;
	FILE	*ipf, *opf;
	keplerorbit	*theorbit;

	if (argc > 1) {
		fprintf(stderr, usage);
		exit(-1);
	}

	buff = (double *) calloc(6, sizeof(double));
	r = buff;
	v = buff + 3;
	theorbit = (keplerorbit *) calloc(sizeof(keplerorbit));

	/* open pipe to input catalog */
	ipf = popen("lc -b -o r v", "r");
	if (!ipf) {
		fprintf(stderr, "testkepler: failed to open pipe for input\n");
		exit(-1);
	}

	/* generate output heaader */
	system("lc -C -b -N '1 3 r_in' -N '1 3 v_in' -N '1 3 r_out' -N '1 3 v_out' -n a -n e -n i -n omega -n Omega -n M < /dev/null");

	while(6 == fread(buff, sizeof(double), 6, ipf)) {
		fwrite(r, sizeof(double), 3, stdout);
		fwrite(v, sizeof(double), 3, stdout);
		cartesiantokepler(r, v, theorbit);
		keplertocartesian(theorbit, r, v);
		fwrite(r, sizeof(double), 3, stdout);
		fwrite(v, sizeof(double), 3, stdout);
		fillbufferfromkeplerorbit(theorbit, buff);
		fwrite(buff, sizeof(double), 6, stdout);
	}

	pclose(ipf);
	exit(0);
}

