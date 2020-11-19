/*
 * testpgplot.c -- quick test of cpgplot stuff
 */


#include "cpgplot.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define error_exit(a) fprintf(stderr, (a)); exit(-1);

main (int argc, char *argv[])
{
	int	nx = 40, ny = 40;

	if (cpgbeg(0, "/xwindow", 1, 1) != 1) {
		error_exit("cpgbeg failed\n");
	}
/*	cpgask(0);*/
	cpgpage();
	cpgsvp(0.05, 0.95, 0.05, 0.95);
	cpgswin(1.0, (float) nx, 1.0, (float) ny);
	cpgbox("bcts", 0.0, 0, "bcts", 0.0, 0);
	cpgend();
	exit(0);
}
