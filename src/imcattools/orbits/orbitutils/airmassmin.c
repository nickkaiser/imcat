#include <stdio.h>
#include <math.h>
#include "vectors.h"
#include "planetdata.h"
#include "obsdata.h"

/* 
	airmassmin() returns the minimum air-mass for a target at rotated ecliptic latitude lat lon
	at phaseofyear given that the minimum solar and lunar angles from the zenith are
	zd_solar_min and zd_lunar_min.  It also returns timeofday in hours (midnight = 0.0)
	at which the minimum air-mass is achieved, and the target angle from the zenith zd_target_min.
	It simply computes the relevant quantity at nstepsperday intervals and keeps track of
	the min airmass.  All double args are in radians.
*/


double	airmassmin(double lat, double lon, double phaseofyear, double zd_solar_min, double zd_lunar_min, 
	int nstepsperday, double *phaseofday, double *zd_target_min)
{
	double	phi, dphi, latobs, axistilt, lunarmonth;
	double	zd_solar, zd_lunar, zd_target, airmass_min;
	static double	*rearth, *robs, *nvert, *ntarget, *nsun, *nmoon;
	static int first = 1;

	/* hard-wired parameters */
	/* lunar month in year */
	lunarmonth = 1.0 / 13.0;
	/* latitude of observatory */
	latobs = LATOBSINDEG * M_PI / 180.0;
	/* tilt of earth's rotation axis */
	axistilt= AXISTILTINDEG * M_PI / 180.0;
	
	/* rotation per time-step */
	dphi = 2 * M_PI / nstepsperday;

	if (first) {
		/* allocate space for vectors */
		rearth	= (double *) calloc(3, sizeof(double));
		robs	= (double *) calloc(3, sizeof(double));
		nvert	= (double *) calloc(3, sizeof(double));
		ntarget	= (double *) calloc(3, sizeof(double));
		nsun	= (double *) calloc(3, sizeof(double));
		nmoon	= (double *) calloc(3, sizeof(double));
		first = 0;
	}

	/* compute the earth's heliocentric coords */
	assign(rearth, cos(phaseofyear), sin(phaseofyear), 0.0);
	/* and the direction to the moon */
	assign(nmoon, cos(phaseofyear / lunarmonth), sin(phaseofyear / lunarmonth), 0.0);
	/* compute the target direction */
	assign(ntarget, cos(lat) * cos(phaseofyear + lon), cos(lat) * sin(phaseofyear + lon), sin(lat));
	/* and the direction to the sun */
	copy(nsun, rearth);
	scale(nsun, -1.0 / length(nsun));
	/* initialise zd_target_min to worst case (on horizon) */
	*zd_target_min = 0.5 * M_PI;
	/* sensible value for min air_mass */
	airmass_min = 99.99;
	/* and for time */
	*phaseofday = 0.0;
	/* loop over times of day */
	for (phi = - M_PI; phi < M_PI; phi += dphi) {
		/* geocentric coords of the observatory */
		assign(robs, cos(latobs) * cos(phaseofyear + phi), cos(latobs) * sin(phaseofyear + phi), sin(latobs));	
		/* tilt the Earth's axis  */
		roty(robs, axistilt);
		/* compute the normal to the Earth's surface */
		copy(nvert, robs);
		scale(nvert, 1.0 / length(nvert));
		/* compute the solar zenith angle */
		zd_solar = acos(dot(nvert, nsun));
		/* compute the lunar zenith angle */
		zd_lunar = acos(dot(nvert, nmoon));
		if ((zd_solar > zd_solar_min) && (zd_lunar > zd_lunar_min)) {
			zd_target = acos(dot(nvert, ntarget));
			if (zd_target < *zd_target_min) {
				*zd_target_min = zd_target;
				/* time in hours */
				*phaseofday = phi;
				/* minimum air mass */
				airmass_min = 1.0 / dot(nvert, ntarget);
			}
		}		
	}
	return(airmass_min);
}

