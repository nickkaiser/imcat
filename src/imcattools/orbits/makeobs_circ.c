#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"
#include "Ffunc.h"
#include "orbitutils/gaussdev.h"

#define usage "\nNAME\n\
	makeobs_circ - make observations for a circular orbit\n\
\n\
SYNOPSIS\n\
	makeobs_circ parfile\n\
\n\
DESCRIPTION\n\
	makeobs_circ generates a catalog containing triplets of observations\n\
	consisting of:\n\
		t     the time t\n\
		re[]  the position of the earth\n\
		ve[]  velocity of the earth\n\
		ra[]  position of the asteroid\n\
		va[]  velocity of the asteroid va\n\
		rho[] the position of the observor relative to re[]\n\
		n[]   the direction vector\n\
\n\
	It generates these according to parameters in 'parfile', another lc cat.\n\
	These are:\n\
		pe    the phase of the earths orbit (degrees, winter = 0.0)\n\
		t0    the local time of the central observation (hours, midnight = 0.0)\n\
		dt    the interval between the observations (hours)\n\
		ra    the radius of the asteroid's orbit (AU)\n\
		ia    the inclination (deg)\n\
		la    the longitude of the ascending node (deg)\n\
		pa    the phase of the asteroid in its orbit (deg)\n\
		sigma the uncertainty astrometry (arcsec)\n\
		nreal the number of realizations\n\
\n\
	Both the earth and asteroid are in circular orbits.\n\
\n\
	Both observatory latitude and tilt of earth's axis are\n\
	hard-wired to 20 deg.\n\
\n\
	We use an idealised model where the year has precisely 360 days\n\
\n\
	Output times are given in units of earths dynamical time (approx 58 days).\n\
\n\
	Distances are in AU.\n\
\n\
SEE ALSO\n\
	makeobs_inertial laplace3\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"


main (int argc, char *argv[])
{
	/* input parameters */
	double	pe, t0, dt, Ra, ia, la, pa, sigma;
	/* output values */
	double	t, *re, *ve, *ra, *va, *rho, *n, tmpd;
	/* internals */
	double	tdyn, Va, Re, etilt, obslat;
	/* indices */
	int	i, it, nt, nreal;
	/* the parameter filename */
	char	*parfilename;
	/* lc-command */
	char	lccom[1024];
	/* pipes for input and output */
	FILE	*ipf, *opf;

	/* parse args */
	if (argc != 2) {
		fprintf(stderr, usage);
		exit(-1);
	}
	parfilename = argv[1];
	if (!strcmp(parfilename, "-u")) {
		fprintf(stderr, usage);
		exit(-1);
	}		

	/* number of opbservations */
	nt = 3;

	/* allocate space for re[3], ve[3], ra[3], va[3], rho[3], n[3] */
	re = (double *) calloc(3, sizeof(double));
	ve = (double *) calloc(3, sizeof(double));
	ra = (double *) calloc(3, sizeof(double));
	va = (double *) calloc(3, sizeof(double));
	rho = (double *) calloc(3, sizeof(double));
	n = (double *) calloc(3, sizeof(double));

	/* read the input parameters */
	sprintf(lccom, "lc -b -o pe t0 dt ra ia la pa sigma nreal < %s", parfilename);
	ipf = popen(lccom, "r");
	if (!ipf) {
		fprintf(stderr, "makeobs_circ : failed to open lc-pipe for input\n");
		exit(-1);
	}
	fread(&pe, sizeof(double), 1, ipf);
	fread(&t0, sizeof(double), 1, ipf);
	fread(&dt, sizeof(double), 1, ipf);
	fread(&Ra, sizeof(double), 1, ipf);
	fread(&ia, sizeof(double), 1, ipf);
	fread(&la, sizeof(double), 1, ipf);
	fread(&pa, sizeof(double), 1, ipf);
	fread(&sigma, sizeof(double), 1, ipf);
	fread(&tmpd, sizeof(double), 1, ipf);
	nreal = (int) tmpd;
	
	pclose(ipf);

	fprintf(stderr, "# makeobs_circ: parameters: pe=%lg t0=%lg dt=%lg ra=%lg ia=%lg la=%lg pa=%lg sigma=%lg nreal=%d\n",
		pe, t0, dt, Ra, ia, la, pa, sigma, nreal);

	/* convert angles to radians */
	pe *= M_PI / 180.0;
	ia *= M_PI / 180.0;
	la *= M_PI / 180.0;
	pa *= M_PI / 180.0;

	/* convert sigma to radians */
	sigma *= M_PI / (180.0 * 60.0 * 60.0);

	/* dynamical time in hours */
	tdyn = 24.0 * 180.0 / M_PI;

	/* asteroid velocity */
	Va = pow(Ra, -0.5);

	/* radius of the earth */
	if (0) {
		Re = 4.e-5;
	} else {
		Re = 0.0;
	}

	/* tilt of rotation axis */
	etilt = M_PI * 20.0 / 180.0;

	/* observatory latitude */
	obslat = M_PI * 20.0 / 180.0;

	/* open lc-pipe for output */
	opf = popen("lc -C -b -n t -n sigma -N '1 3 re' -N '1 3 ve' -N '1 3 ra' -N '1 3 va' -N '1 3 rho' -N '1 3 n'", "w");
	if (!opf) {
		fprintf(stderr, "makeobs_circ : failed to open lc-pipe for output\n");
		exit(-1);
	}

	while (nreal--) {
	    for (it = -1; it <= 1; it++) {
		t = (t0 + it * dt) / tdyn;
		fprintf(opf, "%14.8lg %14.8lg\n", t, sigma);
		/* generate position, velocity vectors */
		assign(re, 1.0, 0.0, 0.0);
		assign(ve, 0.0, 1.0, 0.0);
		assign(ra, Ra, 0.0, 0.0);
		assign(va, 0.0, Va, 0.0);
		/* rotate the earth */
		rotz(re, pe + t);
		rotz(ve, pe + t);
		/* rotate the asteroid along its orbit*/
		rotz(ra, pa + Va * t / Ra);
		rotz(va, pa + Va * t / Ra);
		/* apply the inclination */
		rotx(ra, ia);
		rotx(va, ia);
		/* rotate the ascending node */
		rotz(ra, la);
		rotz(va, la);
		/* compute un-rotated observatory position rho */
		assign(rho, cos(obslat), 0.0, sin(obslat));
		scale(rho, Re);
		/* rotate to proper local time */
		rotz(rho, pe + t + 360.0 * t);
		/* tilt the rotation axis */
		roty(rho, etilt);
		/* compute the angles */
		for (i = 0; i < 3; i++) {
			n[i] = ra[i] - re[i] - rho[i];
		}
		scale(n, 1.0 / length(n));
		/* add error */
		for (i = 0; i < 3; i++) {
			n[i] += gaussdev() * sigma;
		}
		scale(n, 1.0 / length(n));
		/* output the vectors */
		fprintvec(re, opf);
		fprintvec(ve, opf);
		fprintvec(ra, opf);
		fprintvec(va, opf);
		fprintvec(rho, opf);
		fprintvec(n, opf);
	    }
	}

	pclose(opf);

	exit(0);
}

