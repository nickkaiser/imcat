/*
 * planetephem - compute position of a planet at a given time
 *
 * Nick Kaiser 06/24/02
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jpleph.h"

#define usage "\nNAME\n\
	planetephem - compute position of a planet at a given time\n\
\n\
SYNOPSIS\n\
	planetephem timetype time target center\n\
\n\
DESCRIPTION\n\
	planetephem is a wrapper for the JPL subroutine PLEPH.\n\
\n\
	timetype should be either 'JD' (Julian date) or 'MJD' (mean\n\
	Julian date) or 'Y2K' (time in days relative to start of Y2K).\n\
\n\
	The numbering convention for 'target' and 'center' is:\n\
		1 = mercury\n\
		2 = venus\n\
		3 = earth\n\
		4 = mars\n\
		5 = jupiter\n\
		6 = saturn\n\
		7 = uranus\n\
		8 = neptune\n\
		9 = pluto\n\
		10 = moon\n\
		11 = sun\n\
		12 = solar-system barycenter\n\
		13 = earth-moon barycenter\n\
\n\
	If nutations are wanted, set target = 14, center = 0\n\
	For librations, set target = 15, center = 0.\n\
\n\
	Output are the cartesian position (in AU) and velocities\n\
	in (AU per (yr / 2 PI)) in the ecliptic coordinate frame.\n\
\n\
	Planetephem expects to find the path to the binary\n\
	ephemeris file in the environment variable JPL_EPHEMERIDES\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"

/* a day in units of t_dyn */
#define K 0.01720209895

/* MJD of start of Y2K */
#define MJD_Y2K	51544.0
/* Julian Date of MJD=0 */
#define JD0	2400000.5

/* obliqity of the Earth */
#define OBLIQUITYINDEG 23.43929111

main (int argc, char *argv[])
{
	char	*ephemfilename, names[400][6], *timetype;
	double	vals[400], *r, *v, tJD, c, s, y, z;
	int	itarget, icenter, calc_velocity, i, it;
	void	*ephem;

	if (argc != 5) {
		fprintf(stderr, usage);
		exit(-1);
	}
	timetype = argv[1];
	if ((1 != sscanf(argv[2], "%lf", &tJD)) || (1 != sscanf(argv[2], "%d", &itarget)) 
	    || (1 != sscanf(argv[2], "%d", &icenter))) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (!strcmp(timetype, "MJD")) {
		tJD += JD0;
	}
	if (!strcmp(timetype, "Y2K")) {
		tJD += JD0 + MJD_Y2K;
	}

	if (!(ephemfilename = getenv("JPL_EPHEMERIDES"))) {
		fprintf(stderr, "%s : please setenv JPL_EPHEMERIDES to the binary ephemeris file\n", argv[0]);
		exit(-1);
	}

	if (!(ephem = jpl_init_ephemeris( ephemfilename, names, vals))) {
		fprintf(stderr, "%s : jpl_init_ephemeris() failed with ephem file '%s'\n", argv[0], ephemfilename);
		exit(-1);
	}

	/* allocate result */
	r = (double *) calloc(6, sizeof(double));
	v = r + 3;

	itarget = 3;
	icenter = 11;

	jpl_pleph(ephem, tJD, itarget, icenter, r, calc_velocity = 1);
	/* rotate into ecliptic */
	c = cos(OBLIQUITYINDEG * M_PI / 180.0);
	s = sin(OBLIQUITYINDEG * M_PI / 180.0);
	y = r[1];
	z = r[2];
	r[1] = c * y + s * z;
	r[2] = c * z - s * y;
	y = v[1];
	z = v[2];
	v[1] = c * y + s * z;
	v[2] = c * z - s * y;
	for (i = 0; i < 3; i++) {
		fprintf(stdout, "%14.8lg ", r[i]);
	}
	for (i = 0; i < 3; i++) {
		fprintf(stdout, "%14.8lg ", v[i] / K);
	}
	fprintf(stdout, "\n");

	exit(0);
}

