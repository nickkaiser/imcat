#define	usage "\n\n\n\
NAME\n\
	etprofile --- calculates tangential alignment profile\n\
\n\
SYNOPSIS\n\
	etprofile [option...] < catfile > asciifile\n\
		-o io jo	# origin about which we do profile (2048, 2048)\n\
		-d dlnr		# log bin size 0.25\n\
		-r rmin rmax	# min and max radii (200, 2000)\n\
		-l lossfactor	# multiply e by 1/ lossfactor\n\
		-e ename	# name for 2-vector ellipticity (e)\n\
		-x xname	# name for 2-vector spatial coordinate (x)\n\
DESCRIPTION\n\
	\"etprofile\" calculates tangential alignment profile\n\
	from a catalogue.\n\
	Also calculates kappabar = 2 int d ln r eT\n\
	with geometrical boost factor 1 / ( 1 - r^2 / rmax^2).\n\
	Error bars calculated using orthogonal shear component.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"		

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../utils/error.h"
#include "../../imlib/fits.h"

main(int argc, char *argv[])	
{
	int		arg = 1; 
	int		nbins, bin, nobj = 0;
	double		dx, dy, xo, yo, dlnr, rmax, rmin, r, c2, s2, x[2], e[2];
	double		*eTsum, *eT, *nsum, *rbin, *sigmabar, *sigmabarerror;
	double		eX, ee, eesum = 0.0, *eTerror, nu;
	double		lossfactor, boostfactor, rsigma;
	char		defename[32] = "e", defxname[32] = "x", *ename, *xname, lcstring[128];
	char		line[256], argstring[1024];
	FILE		*lcpipe;

	/* defaults */
	xo = yo = 2048;
	dlnr = 0.25;
	rmin = 200.0;
	rmax = 2000;
	lossfactor = 1.0;
	ename = defename;
	xname = defxname;
	
	while (arg < argc) {
		if (argv[arg][0] != '-')
			error_exit(usage);
		switch (argv[arg++][1]) {
			case 'o':
				if (1 != sscanf(argv[arg++], "%lf", &xo))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%lf", &yo))
					error_exit(usage);
				break;
			case 'r':
				if (1 != sscanf(argv[arg++], "%lf", &rmin))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%lf", &rmax))
					error_exit(usage);
				break;
			case 'd':
				if (1 != sscanf(argv[arg++], "%lf", &dlnr))
					error_exit(usage);
				break;
			case 'l':
				if (1 != sscanf(argv[arg++], "%lf", &lossfactor))
					error_exit(usage);
				break;
			case 'e':
				ename = argv[arg++];
				break;
			case 'x':
				xname = argv[arg++];
				break;
			default:
				error_exit(usage);
				break;
		}
	}

	nbins = (int) ceil(log(rmax / rmin) / dlnr);
	dlnr = log(rmax / rmin) / nbins;
	eT = (double *) calloc(nbins + 1, sizeof(double));
	eTerror = (double *) calloc(nbins + 1, sizeof(double));
	eTsum = (double *) calloc(nbins, sizeof(double));
	nsum = (double *) calloc(nbins, sizeof(double));
	rbin = (double *) calloc(nbins, sizeof(double));
	sigmabar = (double *) calloc(nbins + 1, sizeof(double));
	sigmabarerror = (double *) calloc(nbins + 1, sizeof(double));

	for (bin = 0; bin < nbins; bin++) {
		rbin[bin] = rmin * exp((bin + 0.5) * dlnr);
	}
	
	sprintf(lcstring, "lc -o %s %s", xname, ename);
	if (!(lcpipe = popen(lcstring, "r"))) {
		error_exit("etprofile: failed to open lc-pipe for input\n");
	}
	while (fgets(line, 255, lcpipe)) {
		sscanf(line, "%lf %lf %lf %lf", &(x[0]), &(x[1]), &(e[0]), &(e[1]));
		dx = (x[0] - xo);
		dy = (x[1] - yo);
		r = dx * dx + dy * dy;
		if (r > 0.0) {
			r = sqrt(r);
			bin = floor(log(r / rmin) / dlnr);
			/* bin = (bin < 0 ? 0 : bin); */
			if (bin >= 0 && bin < nbins) {
				nobj++;
				nsum[bin] += 1.0;
				e[0] /= lossfactor;
				e[1] /= lossfactor;
				c2 = (dx * dx - dy * dy) / (r * r);
				s2 = 2 * dx * dy / (r * r);
				eX = s2 * e[0] + c2 * e[1];
				eesum += eX * eX;
				eTsum[bin] -= (c2 * e[0] + s2 * e[1]);
			}
		}
	}
	pclose(lcpipe);

	for (bin = nbins - 1; bin >= 0; bin--) {
		if (nsum[bin] > 0) {
			eT[bin] = eTsum[bin] / nsum[bin];
			eTerror[bin] = sqrt(eesum / (nobj * nsum[bin]));
		}
		sigmabar[bin] = sigmabar[bin+1] + 2 * dlnr * eT[bin];
		sigmabarerror[bin] += sigmabarerror[bin+1] + 4 * dlnr * dlnr *
			eTerror[bin] * eTerror[bin];
	}

	argsToString(argc, argv, argstring);
	sprintf(lcstring, "lc -C -n bin -n r -n ngals -n et -n eterror -n rkappa -n kappa -n kappaerror -n nu -x -a '%s'", argstring);
	if (!(lcpipe = popen(lcstring, "w"))) {
		error_exit("etprofile: failed to open lc-pipe for output\n");
	}
	for (bin = 0; bin < nbins; bin++) {
		rsigma = rbin[bin] * exp(-0.5 * dlnr);
		boostfactor = 1.0 / (1 - rsigma * rsigma / (rmax * rmax));
		if (sigmabarerror[bin] > 0) {
			sigmabarerror[bin] = sqrt(sigmabarerror[bin]);
			nu = sigmabar[bin] / sigmabarerror[bin];
		} else {
			nu = 0.0;
		}
		sigmabar[bin] *= boostfactor;
		sigmabarerror[bin] *= boostfactor;
		fprintf(lcpipe, "%5d %10.3lf %5d %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf\n", 
			bin, rbin[bin], (int) nsum[bin], eT[bin], eTerror[bin], 
			rsigma, sigmabar[bin], sigmabarerror[bin], nu);
	}
	pclose(lcpipe);
	exit(0);
}



