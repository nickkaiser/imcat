#define usage "\n\n\n\
NAME\n\
	makemockimage --- generate mock deep CCD image\n\
\n\
SYNOPSIS\n\
	makemockimage [option...]\n\
		-n 	N1 N2		# size of image (1024, 1024)\n\
		-p	pixsize		# size of a pixel in arcsec (0.2)\n\
		-x 	xname		# get position from 2 vector input variable 'xname'\n\
		-phi 	phiname		# get position angle from input variable 'phiname'\n\
		-mu 	muname		# get position angle from input variable 'muname'\n\
		-seed	seed		# seed for random numbers (1)\n\
\n\
DESCRIPTION\n\
	\"makemockimage\" reads a catalogue of containing size\n\
	'theta' and central surface brightness 'csb' (perhaps generated\n\
	by 'makecosmocat') information and generates a fits\n\
	image containing exponential disks of random orientation.\n\
\n\
	By default, makemockimage will generate uniform random position\n\
	x, position angle phi and orientation mu (the latter being the\n\
	cosine of the angle between the disk polar axis and the line\n\
	of sight) but you can use -x, -phi and -mu options to\n\
	read these from the input catalogue instead.\n\
\n\
	For dust free disk galaxies one would expect the surface\n\
	brighness to vary with orientation as 1 / mu.  However, images\n\
	made in this way look quite unrealistic, with way too many bright\n\
	edge on things, so in the default mode we take the 'optically\n\
	thick' approximation and don't scale the csb.\n\
\n\
	If used with makecosmocat with default parameters, the flux of\n\
	objects in the resulting images with number density like that of\n\
	R=24 galaxies (about 5e4 per square degree) is approximately 0.4.\n\
	Equivalently, running 'apphot -z 23.0' on a catalog generated\n\
	from such an image should give sensible galaxy counts.\n\
\n\
SEE ALSO\n\
	makecosmocat\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "imlib/fits.h"
#include "utils/error.h"
#include "utils/arrays.h"

#define  THETA_MAX 8

#define	RAND	((double)rand())/RAND_MAX

double	drand48();
void	addgalaxy(double x0, double y0, double theta, double csb, double phi, double mu, float **f, int N1, int N2);

int		main(int argc, char *argv[])	
{
	int	arg = 1, N1, N2, readxval, readphival, readmuval, ipbuffsize;
	fitsheader	*fits;
	float	**f;
	double	*theta, *csb, ipbuff[6], *x, *y, *phi, *mu, X, Y, PHI, MU, pixsize;
	FILE	*lcpipe;
	char	*xname, *phiname, *muname, lcstring[64];
	int	seed;
	
	/* defaults */
	N1 = 1024;
	N2 = 1024;
	pixsize = 0.2;
	readxval = 0;
	readphival = 0;
	readmuval = 0;
	x = &X;
	y = &Y;
	phi = &PHI;
	mu = &MU;
	seed = 1;
	
	while (arg < argc) {
		if (!strcmp(argv[arg], "-n")) {
			if (++arg >= argc - 1) error_exit(usage);
			sscanf(argv[arg++], "%d", &N1);
			sscanf(argv[arg++], "%d", &N2);
		} else if (!strcmp(argv[arg], "-p")) {
			if (++arg >= argc) error_exit(usage);
			sscanf(argv[arg++], "%lf", &pixsize);
		} else if (!strcmp(argv[arg], "-phi")) {
			if (++arg >= argc) error_exit(usage);
			phiname = argv[arg++];
			readphival = 1;
		} else if (!strcmp(argv[arg], "-mu")) {
			if (++arg >= argc) error_exit(usage);
			muname = argv[arg++];
			readmuval = 1;
		} else if (!strcmp(argv[arg], "-x")) {
			if (++arg >= argc) error_exit(usage);
			xname = argv[arg++];
			readxval = 1;
		} else if (!strcmp(argv[arg], "-seed")) {
			if (++arg >= argc) error_exit(usage);
			sscanf(argv[arg++], "%d", &seed);
		} else {
			error_exit(usage);
		}
	}

	/* seed the random number generator */
	seed++;
	if (seed <= 0) {
		error_exit("makemockimage: seed must be non -ve\n");
	}
	
	srand(seed);

	allocFloatArray(&f, N1, N2);	
	fits = new2Dfitsheader(N1, N2, FLOAT_PIXTYPE);
	add_comment(argc, argv, fits);
	sprintf(lcstring, "lc -b -o theta csb ");
	theta = ipbuff;
	csb = ipbuff + 1;
	ipbuffsize = 2;
	if (readxval) {
		strcat(lcstring, xname);
		strcat(lcstring, " ");
		x = ipbuff + ipbuffsize++;
		y = ipbuff + ipbuffsize++;
	}
	if (readphival) {
		strcat(lcstring, phiname);
		strcat(lcstring, " ");
		phi = ipbuff + ipbuffsize++;
	}
	if (readmuval) {
		strcat(lcstring, muname);
		strcat(lcstring, " ");
		mu = ipbuff + ipbuffsize++;
	}
	lcpipe = popen(lcstring, "r");
	while (fread(ipbuff, ipbuffsize, sizeof(double), lcpipe)) {
		if (!readxval) {
			*x = N1 * RAND;
			*y = N2 * RAND;
		}
		if (!readphival) {
			*phi = 2 * M_PI * RAND;
		}
		if (!readmuval) {
			*mu = RAND;
		}
		addgalaxy(*x, *y, *theta / pixsize, *csb, *phi, *mu, f, N1, N2);
	}
	pclose(lcpipe);
	write2Dfloatimage(f, fits);
	exit(0);
}

void	addgalaxy(double x0, double y0, double theta, double csb, double phi, double mu, float **f, int N1, int N2)
{
	int	ix, iy, ix0, iy0;
	double	dx, dy, x, y, f0;
	
	ix0 = (int) floor(x0);
	iy0 = (int) floor(y0);
	f0 = csb;
	for (iy = iy0 - THETA_MAX * theta; iy <= iy0 + THETA_MAX * theta; iy++) {
		if (iy < 0 || iy >= N2) {
			continue;
		}
		for (ix = ix0 - THETA_MAX * theta; ix <= ix0 + THETA_MAX * theta; ix++) {
			if (ix < 0 || ix >= N1) {
				continue;
			}
			dx = ix + 0.5 - x0;
			dy = iy + 0.5 - y0;
			x = dx * cos(phi) - dy * sin(phi);
			y = dx * sin(phi) + dy * cos(phi);
			f[iy][ix] += f0 * exp(-sqrt(x * x + y * y / (mu * mu)) / theta);
		}
	}
}






