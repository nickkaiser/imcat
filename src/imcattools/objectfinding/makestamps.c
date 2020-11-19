#define	usage "\n\n\n\
NAME\n\
	makestamps --- extract 'postage stamp images'\n\
SYNOPSIS\n\
	makestamps [option...]\n\
		-b b		# box side (32)\n\
		-n 		# normalise images\n\
		-M mag0		# normalisation magnitude\n\
		-f fitsfile	# specify source fits\n\
		-x xname	# name for position 2-vector (x)\n\
		-m magname	# name for magnitude (mag)\n\
		-N normvalname	# normalise to unit 'normvalname'\n\
		-c		# generate catalog with x0[2], x[2], f, i values\n\
\n\
DESCRIPTION\n\
	\"makestamps\" creates a set of tiny images of\n\
	patches of sky around objects in the catalogue.\n\
	It reads a catalogue containing positions and optionally\n\
	magnitudes of nobj objects from stdin and writes to\n\
	stdout a 3-D fits image of dimensions b x b x nobj.\n\
\n\
	It uses the images named in the catalogue header\n\
	by default.  Use -f option to override this.\n\
\n\
	With -n option the surface brightness will be\n\
	scaled by a factor 10^(0.4 * (mag - mag0))\n\
	where mag0 is the magnitude of the first object\n\
	unless you supply a value by hand with -M option.\n\
\n\
	Use -N option to normalise by dividing by 'normval':\n\
	e.g. do\n\
		makestamps -N flux ....\n\
	to normalise to unit flux.\n\
\n\
	Use -c option to generate lc-catalog format output with\n\
	pixel values f, object number i, object coords x0[], and pixel center coords x[]\n\
	with respect to the object coords.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n\n\n"		

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../imlib/fits.h"
#include "../../utils/error.h"
#include "../../utils/arrays.h"
#include "../../utils/ipbuff.h"

#define MAGIC FLOAT_MAGIC

void	getstamp(float **f, int N1, int N2, float **fstamp, int box, int ix, int iy, double ffactor);

main(int argc, char *argv[])	
{
	float 		**f, **fstamp;
	int		arg = 1, N1, N2, needfits, box, normalise, ix0, iy0, ix, iy, iobj, nobj, makecatop; 
	FILE		*fitsipf, *fitsopf, *lcpipe;
	char		lcstring[128] = "", argstring[1024];
	char		*xname, defxname[32] = "x", *magname, defmagname[32] = "mag", line[1024];
	char		sysstring[128], fitsfilename[128];
	double		**ipbuff, *x, mag, mag0, ffactor, opbuff[6];
	int		hasmag0, fluxnormalise;
	fitsheader	*fitsin, *fitsout;
	
	/* defaults */
	box = 32;
	needfits = 1;
	normalise = 0;
	xname = defxname;
	magname = defmagname;
	hasmag0 = 0;
	fluxnormalise = 0;
	makecatop = 0;

	/* parse args */
	while (arg < argc) {
		if (*argv[arg] != '-')
			error_exit(usage);
		switch (*(argv[arg++]+1)) {
			case 'f':
				if (1 != sscanf(argv[arg++], "%s", fitsfilename))
					error_exit(usage);
				needfits = 0;
				break;
			case 'b':
				if (1 != sscanf(argv[arg++], "%d", &box))
					error_exit(usage);
				break;
			case 'n':
				normalise = 1;
				break;
			case 'M':
				if (1 != sscanf(argv[arg++], "%lf", &mag0))
					error_exit(usage);
				hasmag0 = 1;
				normalise = 1;
				break;
			case 'x':
				xname = argv[arg++];
				break;
			case 'm':
				magname = argv[arg++];
				break;
			case 'N':
				fluxnormalise = 1;
				normalise = 1;
				magname = argv[arg++];
				break;
			case 'c':
				makecatop = 1;
				break;
			default:
				error_exit(usage);
				break;
		}
	}
	
	/* read the catalogue */
	if (needfits) {
		strcat(lcstring, "lc -b -P fits_name -o ");
	} else {
		strcat(lcstring, "lc -b -o ");
	}
	strcat(lcstring, xname);
	strcat(lcstring, " ");
	if (normalise) {
		strcat(lcstring, magname);
	}
	if (!(lcpipe = popen(lcstring, "r")))
		error_exit("makestamps: failed to open lc-pipe for input\n");
	if (needfits) {
		fgets(line, 255, lcpipe);
		sscanf(line, "%s", fitsfilename);	
	}
	if (normalise) {
		ipbuff = readdoublebuff(3, lcpipe, &nobj);
	} else {
		ipbuff = readdoublebuff(2, lcpipe, &nobj);
	}

	if (!(fitsipf = fopen(fitsfilename, "r"))) {
		error_exit("makestamps: unable to open fits file for input\n");
	}

	/* read the input image */
	read2Dfloatimage(&f, &N1, &N2, &fitsin, fitsipf);

	if (makecatop) {
		argsToString(argc, argv, argstring);
		sprintf(lcstring, "lc -C -b -n i -n f -N '1 2 x0' -N '1 2 x' -a '%s' -x < /dev/null", argstring);
		system(lcstring);
	} else {
		/* generate the output fitsheader */
		fitsout = copyfitsheader(fitsin);
		if (normalise) {
	        	fitsout->extpixtype = FLOAT_PIXTYPE;
		}
		fitsout->ndim = 3;
		fitsout->n[0] = box;
		fitsout->n[1] = box;
		fitsout->n[2] = nobj;	
		add_comment(argc, argv, fitsout);
		writefitsheader(fitsout);
	}

	/* allocate the postage stamp */
	allocFloatArray(&fstamp, box, box);

	for (iobj = 0; iobj < nobj; iobj++) {
		x = ipbuff[iobj];
		ix0 = (int) floor(x[0]);
		iy0 = (int) floor(x[1]);
		ffactor = 1.0;
		if (normalise) {
			mag = *(ipbuff[iobj] + 2);
			if (fluxnormalise) {
				ffactor = 1.0 / mag;
			} else {
				if ((!iobj) && !hasmag0)  {
					mag0 = mag;
				}
				ffactor = pow(10.0, 0.4 * (mag - mag0));
			}
		}
		getstamp(f, N1, N2, fstamp, box, ix0, iy0, ffactor);
		if (makecatop) {
			opbuff[0] = (double) iobj;
			for (iy = 0; iy < box; iy++) {
				for (ix = 0; ix < box; ix++) {
					if (fstamp[iy][ix] != MAGIC) {
						opbuff[1] = (double) fstamp[iy][ix];
						opbuff[2] = x[0];
						opbuff[3] = x[1];
						opbuff[4] = 0.5 + ix0 + ix - box / 2;
						opbuff[5] = 0.5 + iy0 + iy - box / 2;
						fwrite(opbuff, sizeof(double), 6, stdout);
					}
				}
			}
		} else {
			for (iy = 0; iy < box; iy++) {
				writefitsline(fstamp[iy], fitsout);
			}
		}
	}
	if (!makecatop) {
		writefitstail(fitsout);
	}
	exit(0);
}



void	getstamp(float **f, int N1, int N2, float **fstamp, int box, int x0, int y0, double ffactor)
{
	int	x, y, dx, dy;

	x0 -= box/2;
	y0 -= box/2;

	for (dy = 0; dy < box; dy++) {
		y = y0 + dy;
		for (dx = 0; dx < box; dx++) {
			x = x0 + dx;
			if (y < 0 || y >= N2 || x < 0 || x >= N1) {
				fstamp[dy][dx] = MAGIC;
			} else {
				fstamp[dy][dx] = (f[y][x] == MAGIC ? MAGIC : ffactor * f[y][x]);
			}
		}
	}
}














