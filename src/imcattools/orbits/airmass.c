#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../orbitutils/vectors.h"
#include "../orbitutils/airmassmin.h"

#define usage "\nNAME\n\
	airmass - compute minimum air-mass for observation of a target\n\
\n\
SYNOPSIS\n\
	airmass lon lat min_zd_solar min_zd_lunar dphase nsteps\n\
\n\
DESCRIPTION\n\
	airmass computes, as a function of time,\n\
	the minimum angle from zenith and airmass for a target at\n\
	helio-ecliptic coords lon lat when the solar angle from zenith is\n\
	at least min_zd_solar degrees and the lunar angle from the\n\
	zenith is at least min_zd_lunar degrees.\n\
\n\
	The phase increment dphase is given in degrees relative to winter.\n\
\n\
	We output the time (relative to midnight in hours) and the\n\
	and minimum zenith distance and air-mass of the target.\n\
\n\
	If nsteps is negative, then we output the results for a single\n\
	time = -nsteps * dphase.\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"


main (int argc, char *argv[])
{
	double	phaseofyear, dphase, phase0;
	double	time, zd_solar_min, zd_lunar_min, zd_target_min, airmass_min;
	double	lon, lat, date;
	int	nphi, nsteps, iphase;
	FILE	*opf;

	/* hard-wired parameters */
	nphi = 240;		/* 5 minute steps */

	/* parse args */
	if (argc != 7) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[1], "%lf", &lon) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[2], "%lf", &lat) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[3], "%lf", &zd_solar_min) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[4], "%lf", &zd_lunar_min) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[5], "%lf", &dphase) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[6], "%d", &nsteps) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}

	/* convert angles to radians */
	lon *= M_PI / 180.0;
	lat *= M_PI / 180.0;
	zd_solar_min *= M_PI / 180.0;
	zd_lunar_min *= M_PI / 180.0;
	dphase *= M_PI / 180.0;

	/* in case we're doing a single shot */
	if (nsteps < 0) {
		phase0 = - nsteps * dphase;
		nsteps = 1;
	} else {
		phase0 = 0.0;
	}


	/* open the output file */
	opf = popen("lc -C -n phaseofyear -n time -n zd_min -n air_mass", "w");
	if (!opf) {
		fprintf(stderr, "airmass : failed to open output pipe\n");
		exit(-1);
	}

	for (iphase = 0; iphase < nsteps; iphase++) {
		phaseofyear = phase0 + iphase * dphase;
		airmass_min = airmassmin(lat, lon, phaseofyear, zd_solar_min, zd_lunar_min, nphi, &time, &zd_target_min);
		fprintf(opf, "%14.8lg %14.8lg %14.8lg %14.8lg\n", 
			phaseofyear * 180.0 / M_PI, time * 12.0 / M_PI, zd_target_min * 180.0 / M_PI, airmass_min);
	}
	
	pclose(opf);

	exit(0);
}

