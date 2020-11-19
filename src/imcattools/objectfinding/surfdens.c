#define	usage "\n\n\n\
NAME\n\
	surfdens --- generate a fits image from a catalogue\n\
\n\
SYNOPSIS\n\
	surfdens [option...]\n\
		-x xname		# name for position vector ('x')\n\
		-r x1 x2 y1 y2		# range of coordinates\n\
		-n N1 N2		# image size in pixels\n\
		-w weight		# weight points by this value\n\
		-s 			# generate 16 bit image\n\
		-d 			# print input lc filter string an quit\n\
		-m			# set points to MAGIC value\n\
		-f fitsfile		# supply source image\n\
\n\
DESCRIPTION\n\
	\"surfdens\" reads a catalogue from stdin and calculates\n\
	a surface density map by binning counts onto an image.\n\
	By default it looks for 2-vector entry 'x', gets the image size\n\
	from the catalogue header, sets the range of\n\
	coordinates to be 0-N1, 0-N2, and uses unit weight per object.\n\
	Specify name for weight with -w option.\n\
	Fits image is output to stdout.\n\
\n\
	With the -m option, occupied pixels are set to MAGIC value and\n\
	-w option, if present, is ignored.\n\
\n\
	With the -f option we initialise the image to 'fitsfile' and\n\
	the -n and -s options, if present, are ignored (the output image\n\
	inheriting the pixel type of the source image.\n\
\n\
NOTES\n\
	If you want to use an rpn expression in the input lc filter\n\
	you will need to double quote it.  E.g.:\n\
		surfdens -w \"'w = %%x[0]'\"\n\
	to get a weight proportional to x-position.\n\
	To use a pair of scalar catalogue object values do e.g.:\n\
		surfdens -x 'fs nu'.\n\
\n\
	See also 'makedensity' which can handle arbitrary dimension\n\
	FITS file and position variables.\n\
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
#include "../../utils/arrays.h"

#define MAGIC	FLOAT_MAGIC


main(int argc, char *argv[])	
{
	int		arg = 1;
	int		useweight, needfitssize, setscale, debug, extpixtype;
	int		N1, N2, ix, iy;
	float		**f;
	fitsheader	*fits;
	double		x, y, w, x1, x2, y1, y2, xscale, yscale;
	char		defaultposname[2] = "x", *posname, *weightname, lcstring[128], line[128];
	FILE		*lcpipe, *fitsfile;
	double		ipbuffer[3];
	int		nel, setmagic, readfitsfile;
	char		*fitsfilename;

	/* defaults */
	useweight = 0;
	extpixtype = FLOAT_PIXTYPE;
	needfitssize = 1;
	posname = defaultposname;
	setscale = 0;
	debug = 0;
	x1 = y1 = 0.0;
	w = 1.0;
	readfitsfile = 0;
	setmagic = 0;

	/* get the optional arguments */
	while (arg < argc) {
		if (argv[arg][0] == '-') {			/* an option argument */
			switch (argv[arg++][1]) {
				case 'x':
					posname = argv[arg++];
					break;					
				case 'r':
					if (1 != sscanf(argv[arg++], "%lf", &x1))
						error_exit(usage);
					if (1 != sscanf(argv[arg++], "%lf", &x2))
						error_exit(usage);
					if (1 != sscanf(argv[arg++], "%lf", &y1))
						error_exit(usage);
					if (1 != sscanf(argv[arg++], "%lf", &y2))
						error_exit(usage);
					if ((y2 == y1) || (x2 == x1))
						error_exit("surfdens: zero range specified\n");
					setscale = 1;
					break;					
				case 'n':
					if (1 != sscanf(argv[arg++], "%d", &N1))
						error_exit(usage);
					if (1 != sscanf(argv[arg++], "%d", &N2))
						error_exit(usage);
					needfitssize = 0;
					break;					
				case 'w':
					useweight = 1;
					weightname = argv[arg++];
					break;					
				case 's':
					extpixtype = SHORT_PIXTYPE;
					break;					
				case 'd':
					debug = 1;
					break;					
				case 'm':
					setmagic = 1;
					break;
				case 'f':
					readfitsfile = 1;
					needfitssize = 0;
					fitsfilename = argv[arg++];
					break;				
				default:
					error_exit(usage);
					break;
			}
		} else {
			error_exit(usage);
		}
	}
	if (setmagic) {
		useweight = 0;
	}

	/* create the input filter string */
	strcpy(lcstring, "lc -b -o ");
	if (needfitssize)
		strcat(lcstring, "-P fits_size ");
	strcat(lcstring, posname);
	if (useweight) {
		strcat(lcstring, " ");
		strcat(lcstring, weightname);
	}

	if (debug) {
		fprintf(stdout, "%s\n", lcstring);
		exit(0);
	}

	if (!(lcpipe = popen(lcstring, "r")))
		error_exit("surfdens: failed to open lc for input filter\n");

	if (needfitssize)
		fgets(line, 128, lcpipe);
	sscanf(line, "%d %d", &N1, &N2);

	/* set the scale */
	if (setscale) {
		xscale = N1 / (x2 - x1);
		yscale = N2 / (y2 - y1);
	} else {
		xscale = yscale = 1.0;
	}

	/* create the output image */
	if (readfitsfile) {
		fitsfile = fopen(fitsfilename, "r");
		if (!fitsfile) {
			error_exit("surfdens: failed to open source image\n");
		}
		read2Dfloatimage(&f, &N1, &N2, &fits, fitsfile);
	} else {	
		fits = new2Dfitsheader(N1, N2, extpixtype);
		allocFloatArray(&f, N1, N2);
	}
	
	nel = (useweight ? 3 : 2);	
	while(fread(ipbuffer, sizeof(double), nel, lcpipe)) {
		if (useweight) {
			w = ipbuffer[2];
		}
		x = ipbuffer[0];
		y = ipbuffer[1];		
		ix = floor(xscale * (x - x1));
		iy = floor(yscale * (y - y1));
		if (ix < 0 || ix >= N1 || iy < 0 || iy >= N2)
			continue;
		if (setmagic) {
			f[iy][ix] = MAGIC;
		} else {
			f[iy][ix] += w;
		}
	}
	add_comment(argc, argv, fits);
	write2Dfloatimage(f, fits);
	exit(0);
}


