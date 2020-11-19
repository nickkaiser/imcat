#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "imlib/fits.h"
#include "utils/ipbuff.h"
#include "utils/arrays.h"

#define usage "\nNAME\n\
	orbanim - generate an animation of orbit trajectories\n\
\n\
SYNOPSIS\n\
	orbanim np nframes Nx Ny [-r rname] [-d decay_factor] [-b bgfits]\n\
\n\
DESCRIPTION\n\
	orbabim reads from stdin a catalog containing a set of nframes\n\
	sets of np particles with coordinates x[2].  It generates a 3D fits\n\
	image f[nframes][N][N] showing the positions of the particles\n\
	in the frames.\n\
\n\
	The -r option is used to modulate the brightness of points.\n\
	If r is less than 1, a single pixel is painted with the value r\n\
	and if r is greater than one, a disk of radius r is painted.\n\
\n\
	With -d option, displayed intensity decays exponentially, with\n\
	successive frames reduced by a factor decay_factor.\n\
\n\
	With -b option we use the image 'bgfits' as a background image.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"


main (int argc, char *argv[])
{
	double	*ipbuff, r;
	int	np, nframes, Nx, Ny, ip, i, x, y, ix, iy, Dx, Dy, dim[3];
	int	doradius, dodecay, arg, planet_radius, dobg;
	double	decay_factor;
	float	**f, **fb;
	FILE	*ipf;
	char	lccom[1024], opt, *rname, *bgfitsname;
	fitsheader	*fits;

	/* defaults */
	doradius = 0;
	dodecay = 0;
	dobg = 0;

	/* parse args */
	if ((argc < 5)) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[1], "%d", &np) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[2], "%d", &nframes) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[3], "%d", &Nx) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[4], "%d", &Ny) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}
	arg = 5;
	while (arg < argc) {
		if (strncmp(argv[arg], "-", 1)) {
			fprintf(stderr, usage);
			exit(-1);
		}
		opt = argv[arg][1];
		switch (opt) {
			case 'd':
				arg++;
				if (arg >= argc) {
					fprintf(stderr, usage);
					exit(-1);
				}
				if (1 != sscanf(argv[arg], "%lf", &decay_factor)) {
					fprintf(stderr, usage);
					exit(-1);
				} 
				dodecay = 1;
				break;
			case 'r':
				arg++;
				if (arg >= argc) {
					fprintf(stderr, usage);
					exit(-1);
				}
				rname = argv[arg];
				doradius = 1;
				break;
			case 'b':
				arg++;
				if (arg >= argc) {
					fprintf(stderr, usage);
					exit(-1);
				}
				bgfitsname = argv[arg];
				dobg = 1;
				break;
			default:
				fprintf(stderr, "unknown option: %s\n", argv[arg]);
				exit(-1);
		}
		arg++;
	}

	if (dobg) {
		/* read the background fits image */
		ipf = fopen(bgfitsname, "r");
		if (!ipf) {
			fprintf(stderr, "orbanim : failed to open background fits image for input\n");
			exit(-1);
		}
		fits = readfitsheader(ipf);
		if ((fits->ndim != 2) || (fits->n[0] != Nx) || (fits->n[1] != Ny)) {
			fprintf(stderr, "orbanim : background image has wrong size\n");
		}
		allocFloatArray(&fb, Nx, Ny);
		readfitsplane((void **) fb, fits);
		fclose(ipf);
	}

	/* open lc pipe for input */
	if (doradius) {
		sprintf(lccom, "lc -b -o x %s", rname);
	} else {
		sprintf(lccom, "lc -b -o x");
	}
	ipf = popen(lccom, "r");
	if (!ipf) {
		fprintf(stderr, "orbanim : failed to open lc-pipe for input\n");
		exit(-1);
	}

	ipbuff = (double *) calloc((doradius ? 3 : 2), sizeof(double));

	allocFloatArray(&f, Nx, Ny);

	dim[0] = Nx;
	dim[1] = Ny;
	dim[2] = nframes;
	fits = newfitsheader(3, dim, FLOAT_PIXTYPE);
	writefitsheader(fits);	

	for (i = 0; i < nframes; i++) {
		for (y = 0; y < Ny; y++) {
			for (x = 0; x < Nx; x++) {
				if (dodecay) {
					f[y][x] *= decay_factor;
				} else {
					f[y][x] = 0.0;
				}
				if (dobg) {
					if (fb[y][x] > 0.0) {
						f[y][x] = fb[y][x];
					}
				}
			}
		}
		for (ip = 0; ip < np; ip++) {
			fread(ipbuff, sizeof(double), (doradius ? 3 : 2), ipf);
			x = (int) floor(ipbuff[0] + 0.5);
			y = (int) floor(ipbuff[1] + 0.5);
			if (doradius) {
				r = ipbuff[2];
				planet_radius = (int) floor(r);
			} 
			if (doradius && planet_radius) {
				for (Dx = -planet_radius; Dx <= planet_radius; Dx++) {
					for (Dy = -planet_radius; Dy <= planet_radius; Dy++) {
						if (sqrt(Dx * Dx + Dy * Dy) < r) {
							ix = x + Dx;
							iy = y + Dy;
							if (ix >= 0 && ix < Nx && iy >= 0 && iy < Ny) {
								f[iy][ix] = 1.0;
							}
						}
					}
				}
			} else {
				if (x >= 0 && x < Nx && y >= 0 && y < Ny) {
					f[y][x] = (doradius ? r : 1.0);
				}
			}
		}
		writefitsplane((void **) f, fits);
	}
	writefitstail(fits);
	pclose(ipf);
	exit(0);
}

