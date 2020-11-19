#define	usage "\n\n\n\
NAME\n\
	findpeaks --- simple object finder\n\
SYNOPSIS\n\
	findpeaks fitsimage [option...] >catalogue \n\
		-r	sigma	# gaussian filter radius (default = 2)\n\
		-n 	nu	# threshold (default = 3)\n\
		-e		# find all extrema\n\
		-m		# find all minima\n\
		-s sigma mode	# sky statistics\n\
		-d 		# add 2nd derivative information\n\
		-o		# output FITS file followed by catalog\n\
\n\
DESCRIPTION\n\
	\"findpeaks\" fft gaussian filters image and finds peaks above a threshold nu.\n\
	Creates a catalogue with limited information\n\
	With -e or -m option, nu parameter ignored.\n\
	The position 'x' is measured relative to the bottom left\n\
	corner of the bottom left pixel (so e.g. a single 'hot' pixel\n\
	at (ix,iy) = (23,67), would generate an object with x = (23.5, 67.5)\n\
\n\
	Use -r option to control the smoothing radius.  With negative\n\
	sigma we don't smooth at all.\n\
\n\
	'findpeaks outputs\n\
		x[2]		# peak position\n\
		fs		# smoothed image value at the peak\n\
		nu		# significance value\n\
		maximum		# true if extremum is a maximum\n\
\n\
	and with '-d' option it also outputs\n\
		ddfs[2][2]	# discretized 2nd derivative at peak\n\
		detddfs		# determinant\n\
\n\
	'findpeaks' uses iostream library.\n\
\n\
	With -o option we send the source image to stdout (but without\n\
	the tail) followed by the catalog.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"		

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../imlib/fits.h"
#include "../../utils/error.h"
#include "../../imlib/filters.h"
#include "../../utils/arrays.h"
#include "../../utils/stats_stuff.h"
#include "../../utils/iostream.h"
#include "../../catlib/cat.h"

#define	DEFAULT_NU	3
#define	DEFAULT_RF	2

#define FINDPEAKS	0
#define FINDEXTREMA	1
#define FINDMINIMA	2

#define MAGIC FLOAT_MAGIC

#define PI		M_PI

#define TINY	1.e-20

main(int argc, char *argv[])	
{
	float 		**f, **fs, f_thresh, magicsubstitute, fn, fsn;
	float		nu, rf, Fp, Fpp;
	int		arg = 2, N1, N2, margin = 3; 
	fstatsrec	srec;
	fitsheader	*fits;
	char		*fitsfilename, lcstring[1024], tempstring[128];
	int		i, j;
	FILE		*fitsfile, *lcpipe;
	cathead		*thecathead;
	object		*theobject;
	item		*filenameitem, *hasskyitem, *filesizeitem;
	int		n, mode, hasbiggerneighbour, hassmallerneighbour, hasmagicneighbour;
	int		in[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
	int		jn[8] = {-1, 0, 1, 1, 1, 0, -1, -1};
	double		R[2], FS, NU, **ddfs, detddfs, maximum;
	int		calcstats = 1, doddfs, outputimage;
	iostream	*theiostream;
	
	/* defaults */
	nu = DEFAULT_NU;
	rf = DEFAULT_RF;
	mode = FINDPEAKS;
	doddfs = 0;
	outputimage = 0;
	
	if (argc < 2)
		error_exit(usage);
	if (!strcmp(argv[1], "-u")) 
		error_exit(usage);
	while (arg < argc) {
		if (*argv[arg] != '-')
			error_exit(usage);
		switch (*(argv[arg++]+1)) {
			case 'r':
				if (1 != sscanf(argv[arg++], "%f", &rf))
					error_exit(usage);
				break;
			case 'n':
				if (1 != sscanf(argv[arg++], "%f", &nu))
					error_exit(usage);
				break;
			case 'e':
				mode = FINDEXTREMA;
				break;
			case 'm':
				mode = FINDMINIMA;
				break;
			case 's':
				if (1 != sscanf(argv[arg++], "%f", &(srec.sigma)))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%f", &(srec.fmode)))
					error_exit(usage);
				calcstats = 0;
				break;
			case 'd':
				doddfs = 1;
				break;
			case 'o':
				outputimage = 1;
				break;
			default:
				error_exit(usage);
				break;
		}
	}

	fitsfilename = argv[1];	
	theiostream = openiostream(fitsfilename, "r");
	fitsfile = theiostream->f;
	if (!fitsfile)
		error_exit("findpeaks: can't open fits file\n");
	read2Dfloatimage(&f, &N1, &N2, &fits, fitsfile);

	if (outputimage) {
		add_comment(argc, argv, fits);
		writefitsheader(fits);
		writefitsplane((void **) f, fits);
	}

	if (calcstats) {
		fdo_stats(f, N1, N2, margin, &srec);
	}
	if (rf > 0.0) {
		allocFloatArray(&fs, N1, N2);
		magicsubstitute = srec.fmode;
		gaussfilter(f, N1, N2, fs, rf, rf, 0, magicsubstitute);
		/* reset magic pixels */
		for (i = 0; i < N2; i++) {
			for (j = 0; j < N1; j++) {
				if (f[i][j] == MAGIC) {
					fs[i][j] = MAGIC;
				}
			}
		}
	} else {
		fs = f;
	}

	margin = 0;

	if (calcstats)
		fdo_stats(fs, N1, N2, margin, &srec);
	else		
		srec.sigma /= 2 * sqrt(PI) * rf;
	f_thresh = srec.fmode + nu * srec.sigma;

	/* simplest to just create the cathead with 'lc */
	if (theiostream->type == FILE_IOSTREAM_TYPE) {
		sprintf(lcstring, "lc -C -x -H 'fits_name = {%s}' ", fitsfilename);
	} else {
		sprintf(lcstring, "lc -C -x -H 'fits_name = {iostream_input}' ");
	}
	sprintf(tempstring, "-H 'fits_size = %d %d 2 vector' -H 'has_sky = 0' ", N1, N2);
	strcat(lcstring, tempstring);
	if (doddfs) {
		sprintf(tempstring, "-N '1 2 x' -N '1 1 fs' -N '1 1 nu' -N '2 2 2 ddfs' -n detddfs -n maximum < /dev/null");
	} else {
		sprintf(tempstring, "-N '1 2 x' -N '1 1 fs' -N '1 1 nu' -n maximum < /dev/null");
	}
	strcat(lcstring, tempstring);
	lcpipe = popen(lcstring, "r");
	setcatipf(lcpipe);
	thecathead = readcathead();
	pclose(lcpipe);

	/* now add the history */
	addargscomment(argc, argv, thecathead);

	/* and write the cathead */
	setcatopfiletype(BINARY_FILE_TYPE);
	writecathead(thecathead);

	/* now create the object */
	theobject = newobject(thecathead);
	connectcatheadtoobject(theobject);
	setaddress(theobject, getobjectitemindex("x", theobject), R);
	setaddress(theobject, getobjectitemindex("fs", theobject), &FS);
	setaddress(theobject, getobjectitemindex("nu", theobject), &NU);
	setaddress(theobject, getobjectitemindex("maximum", theobject), &maximum);
	if (doddfs) {
		setaddress(theobject, getobjectitemindex("detddfs", theobject), &detddfs);	
		ddfs = (double **) calloc(2, sizeof(double *));
		ddfs[0] = (double *) calloc(2, sizeof(double));
		ddfs[1] = (double *) calloc(2, sizeof(double));
		setaddress(theobject, getobjectitemindex("ddfs", theobject), ddfs);
	}

	/* and now find peaks */
	for (i = 1; i < N2 - 1; i++) {
		for (j = 1; j < N1 - 1; j++) {
		        if (mode == FINDPEAKS && fs[i][j] < f_thresh)
			     	continue;
			if (f[i][j] == MAGIC)
				continue;
			hasbiggerneighbour = hassmallerneighbour = hasmagicneighbour = 0;
			for (n = 0; n < 8; n++) {
				fsn = fs[i + in[n]][j + jn[n]];
				fn = f[i + in[n]][j + jn[n]];
				if (fn == (float) MAGIC) {
					hasmagicneighbour = 1;
					break;
				}
				if (fsn > fs[i][j]) {
					hasbiggerneighbour = 1;
					if (hassmallerneighbour)
						break;
				}
				if (fsn < fs[i][j]) {
					hassmallerneighbour = 1;
					if (hasbiggerneighbour)
						break;
				}
			}
			if (hasmagicneighbour)
				continue;
			if (hasbiggerneighbour && hassmallerneighbour)
				continue;
			if (mode == FINDPEAKS && hasbiggerneighbour)
				continue;
			if (mode == FINDMINIMA && hassmallerneighbour)
				continue;
			R[0] = (double) j + 0.5;
			R[1] = (double) i + 0.5;
			if (doddfs) {
				ddfs[0][0] = 2 * fs[i][j] - fs[i+1][j] - fs[i-1][j];
				ddfs[1][1] = 2 * fs[i][j] - fs[i][j+1] - fs[i][j-1];
				ddfs[0][1] = ddfs[1][0] = 0.25 * (fs[i+1][j+1] - fs[i+1][j-1] - fs[i-1][j+1] + fs[i-1][j-1]);
				detddfs = ddfs[0][0] * ddfs[1][1] - ddfs[0][1] * ddfs[1][0];
			}
			/* now we add the simple cludge to get more precise position */
			Fp  = 0.5 * (fs[i][j+1] - fs[i][j-1]);
			Fpp = fs[i][j+1] - 2 * fs[i][j] + fs[i][j-1];
			if (fabs(Fpp) > TINY)
				R[0] -= Fp / Fpp;
			Fp  = 0.5 * (fs[i+1][j] - fs[i-1][j]);
			Fpp = fs[i+1][j] - 2 * fs[i][j] + fs[i-1][j];
			if (fabs(Fpp) > TINY)
				R[1] -= Fp / Fpp;
			FS = fs[i][j];
			NU = (fs[i][j] - srec.fmode) / srec.sigma;
			maximum = (double) hassmallerneighbour;
			writeobject(theobject);
		}
	}
	exit(0);
}










