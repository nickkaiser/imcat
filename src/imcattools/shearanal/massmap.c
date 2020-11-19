#define	usage "\n\n\n\
NAME\n\
	massmap --- direct summation massmap\n\
SYNOPSIS\n\
	massmap [option...] \n\
		-g ng			# output grid size (64)\n\
		-s rs			# gaussian smoothing radius in output grid size units (2.0)\n\
		-n N			# input image size	(2048)\n\
		-o xo yo		# origin of cat in pixels (0,0)\n\
		-R R			# radius for determining n_bar (N/4)\n\
		-e pol			# output smoothed (n * e[pol] / nbar) map\n\
		-d 			# output D.s\n\
		-c 			# output D x s\n\
		-l lossfactor		# divide final mass-map by lossfactor (1.0)\n\
		-E ename		# name for 2-vector ellipticity (e)\n\
		-x xname		# name for 2-vector spatial coordinate (x)\n\
\n\
DESCRIPTION\n\
	\"massmap\" reads e[2], x[2] from a catalogue and calculates foreground\n\
	surface mass density field from background galaxy ellipticities a la KS.\n\
	Calculates mean galaxy number density in disk radius R around field centre\n\
	Outputs Sigma / Sigma_crit unless -e option set in which case it\n\
	outputs a gaussian smoothed map of (n * e[pol] / nbar) for pol = 0 or 1\n\
	Use lossfactor (< 1.0) to correct for seeing\n\
	Uses W(theta) = (0.1 y^4 / (1 + 0.1 y^4)) / theta^2 (where y = theta / sigma)\n\
	as truly excellent approximation to bessel function window for gaussian T(k)\n\
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

#define	TINY	1.e-8

#define DO_MASS_MAP	0
#define DO_E_MAP	1
#define DO_DOT_MAP	2
#define DO_CROSS_MAP	3

#define PI M_PI


main(int argc, char *argv[])	
{
	int			arg = 1;
	FILE		*catf;
	char		errorstring[1024];
	int		ng, N, xg, yg;
	double		xo, yo, x, y, e[2];
	int		useweight;
	double		rs;								/* smoothing radius */
	float		**sigma;							/* surface density */
	double		scale;								/* grid points = pixels x scale */
	double		dd, dx, dy;							/* separation in grid units */
	double		w;								/* the weight */
	double		z, ee, eT;							/* e, e^2 and tangential ellipticity */
	double		cos2, sin2;
	int		ngood;	
	double		R, nbar, lossfactor, g;
	int		pol, mode;
	char		defename[32] = "e", defxname[32] = "x", *ename, *xname, lcstring[128];
	char		line[256];
	FILE		*lcpipe;
	fitsheader	*fits;

	/* defaults */
	ng = 64;
	N = 2048;
	xo = yo = 0;
	rs = 2.0;
	R = 0.25 * N;
	mode = DO_MASS_MAP;
	lossfactor = 1.0;
	xname = defxname;
	ename = defename;
	
	/* get the optional arguments */
	while (arg < argc) {
		if (argv[arg][0] == '-') {			/* an option argument */
			switch (argv[arg++][1]) {
				case 'g':
					if (1 != sscanf(argv[arg++], "%d", &ng))
						error_exit(usage);
					break;					
				case 's':
					if (1 != sscanf(argv[arg++], "%lf", &rs))
						error_exit(usage);
					break;					
				case 'n':
					if (1 != sscanf(argv[arg++], "%d", &N))
						error_exit(usage);
					break;					
				case 'o':
					if (1 != sscanf(argv[arg++], "%lf", &xo))
						error_exit(usage);
					if (1 != sscanf(argv[arg++], "%lf", &yo))
						error_exit(usage);
					break;					
				case 'R':
					if (1 != sscanf(argv[arg++], "%lf", &R))
						error_exit(usage);
					break;					
				case 'l':
					if (1 != sscanf(argv[arg++], "%lf", &lossfactor))
						error_exit(usage);
					break;					
				case 'd':
					mode = DO_DOT_MAP;
					break;					
				case 'c':
					mode = DO_CROSS_MAP;
					break;					
				case 'e':
					mode = DO_E_MAP;
					if (1 != sscanf(argv[arg++], "%d", &pol))
						error_exit(usage);
					if (pol != 0 && pol != 1)
						error_exit(usage);
					break;					
				case 'E':
					ename = argv[arg++];
					break;
				case 'x':
					xname = argv[arg++];
					break;
				default:
					error_exit(usage);
					break;
			}
		} else {
			error_exit(usage);
		}
	}

	/* set option dependent variables */
	scale = (double) ng / (double) N;

	/* set up the output image */
	fits = new2Dfitsheader(ng, ng, FLOAT_PIXTYPE);
	allocFloatArray(&sigma, ng, ng);
	
	/* read catalogue and accumulate weighted sum of eT for each grid point */
	/* also accumulate gals in central disk */
	ngood = 0;
	nbar = 0.0;
	sprintf(lcstring, "lc -o %s %s", xname, ename);
	if (!(lcpipe = popen(lcstring, "r"))) {
		error_exit("etprofile: failed to open lc-pipe for input\n");
	}
	while (fgets(line, 255, lcpipe)) {
		ngood ++;
		if (4 != sscanf(line, "%lf %lf %lf %lf", &x, &y, &(e[0]), &(e[1])))
			error_exit("massmap: input format error\n");
		if (((x - N / 2) * (x - N / 2) + (y - N / 2) * (y - N / 2)) < (R * R))
			nbar += 1.0; 
		for (yg = 0; yg < ng; yg++) {
			for (xg = 0; xg < ng; xg++) {
				dx = scale * (x + xo) - xg;
				dy = scale * (y + yo) - yg;
				dd = dx * dx + dy * dy + TINY;
				if ((dd < 16 * rs * rs) && mode != DO_MASS_MAP)
					g = exp(-0.5 * dd / (rs * rs)) / (2 * PI * rs * rs);
				if (mode != DO_MASS_MAP && (dd >= 16 * rs * rs))
					continue;
				switch (mode) {
					case DO_MASS_MAP:
						cos2 = (dx * dx - dy * dy) / dd;
						sin2 = 2 * dx * dy / dd;
						eT = - (cos2 * e[0] + sin2 * e[1]);
						z = dd / (rs * rs);
						w = 0.1 * z / ((1 + 0.1 * z * z) * rs * rs);
						sigma[yg][xg] += eT * w;
						break;
					case DO_E_MAP:
						sigma[yg][xg] += e[pol] * g;
						break;
					case DO_DOT_MAP:
						sigma[yg][xg] += g * ((dx * dx - dy * dy) * e[0]
							+ 2 * dx * dy * e[1]) / (rs * rs);
						break;
					case DO_CROSS_MAP:
						sigma[yg][xg] += g * ((dx * dx - dy * dy) * e[1]
							- 2 * dx * dy * e[0]) / (rs * rs);
						break;
					default:
						error_exit("massmap: bad mode\n");
				}
			}
		}
	}
	nbar /= PI * (R * R);
	
	/* make scaled version for output*/
	for (yg = 0; yg < ng; yg++) {
		for (xg = 0; xg < ng; xg++) {
			switch (mode) {
				case DO_MASS_MAP:
					sigma[yg][xg] *= scale * scale / (nbar * PI);
					break;
				default:
					sigma[yg][xg] *= scale * scale / nbar;
			}
			sigma[yg][xg] /= lossfactor;
		}
	}

	/* and output */
	add_comment(argc, argv, fits);
	write2Dfloatimage(sigma, fits);
	return(0);
}


#undef PI

