#define usage "\n\n\n\
NAME\n\
	otf2psf --- compute PSF from makeotf output\n\
SYNOPSIS\n\
	otf2psf dtheta thetamax [-l lambda]\n\
\n\
DESCRIPTION\n\
	'otf2psf' reads a file containing z gk(z) from stdin and\n\
	computes\n\
\n\
		g(theta) = 2 pi int dz z J0(2 pi theta z / lambda) gk(z) / lambda^2\n\
\n\
	for theta = 0 - thetamax with increment dtheta.\n\
	The quantity output is g(theta) x (1 radian / 1 arcsec)^2.\n\
	Command line args are given in units of arcseconds.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "utils/error.h"
#include "utils/args.h"
#include "utils/ipbuff.h"
#include "imlib/fits.h"

int		main(int argc, char *argv[])	
{
	char		*flag, argstring[256];
	char		lcstr[256];
	int		nz, iz, ix;
	double		**ipbuff, theta, dtheta, thetamax, lambda, arcsec, *z, g, dz, *gk;
	FILE		*ipf;

	/* defaults */
	lambda = 8.e-7;

	/* parse args */
	argsinit(argc, argv, usage);
	dtheta 		= getargd();
	thetamax 	= getargd();
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'l' :
				lambda = getargd();
				break;
			default:
				error_exit(usage);
		}
	}

	arcsec = M_PI / (180.0 * 3600.0);

	/* open the output stream */
	argsToString(argc, argv, argstring);
	sprintf(lcstr, "lc -C -x -a 'history: %s' -n theta -n g < /dev/null", argstring);
	system(lcstr);

	/* read the input buffer */
	ipf = popen("lc -b -o z g", "r");
	if (!ipf) {
		error_exit("otf2psf: failed to open input pipe\n");
	}
	ipbuff = readdoublebuff(2, ipf, &nz);
	z = (double *) calloc(nz, sizeof(double));
	gk = (double *) calloc(nz, sizeof(double));
	for (iz = 0; iz < nz; iz++) {
		z[iz] = ipbuff[iz][0];
		gk[iz] = ipbuff[iz][1];	
	}
	dz = z[1];

	theta = 0.0;
	while (theta < thetamax) {
		g = 0.0;
		for (iz = 0; iz < nz; iz++) {
			g += dz * z[iz] * j0(2 * M_PI * theta * arcsec * z[iz] / lambda) * gk[iz];
		}
		fprintf(stdout, "  %14.8lg %14.8lg\n", theta, 2 * M_PI * g * arcsec * arcsec / (lambda * lambda));
		theta += dtheta;
	}

	exit(0);
}

