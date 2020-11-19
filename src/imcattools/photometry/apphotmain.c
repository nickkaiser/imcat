/*
 *	apphotmain.c - front end largely copies from getsky.c
 */

#define	usage "\n\n\n\
NAME\n\
	apphot --- perform aperture photometry\n\
\n\
SYNOPSIS\n\
	apphot	[options...]\n\
		-m mul		# r_aperture = mul * r (3.0)\n\
		-r rname	# name for aperture radius ('rg')\n\
		-R rap		# fixed aperture\n\
		-f fitsfile	# specify fits file explicitly\n\
		-i 		# ignore local sky parameters\n\
		-z zeropoint	# zeropoint for magnitude scale\n\
		-Z r mul	# ignore pixels within mul * r of other objects\n\
		-M rapmax	# maximum aperture for photometry (50)\n\
\n\
DESCRIPTION\n\
	By default, 'apphot' measures flux and mag for aperture of radius\n\
	mul * rg as determined by hfindpeaks, though you can use the -r\n\
	(and optionally the -m option) to choose an alternative, or use\n\
	the -R option to specify a fixed aperture.\n\
\n\
	It adds to the catalogue objects the following items\n\
\n\
		flux	# flux (sum of pixel values) within the aperture\n\
		mag	# magnitude = zeropoint - 2.5 * log_10(flux)\n\
		rh	# radius within which 1/2 of flux is found\n\
		rp	# Petrosian radius\n\
		rql	# radius within which 1/4 of flux is found\n\
		rqu	# radius within which 3/4 of flux is found\n\
		nbad	# number of magic pixels within the aperture\n\
		fmax	# highest pixel value within aperture.   \n\
\n\
	The 'Petrosian' radius 'rp' is defined to be the first\n\
	maximum in the cumulative enclosed flux divided by the radius.\n\
	In order to avoid assigning an unreasonably small petrosian radii\n\
	to small objects where the centroid happens to lie very close to\n\
	the centre of a pixel we 'soften' the radius: r -> sqrt(r^2 + 0.3^2)\n\
\n\
	Apphot takes all of the pixels whose centres fall within a distance\n\
	less than r_aperture, sorts them by (softened) distance, and computes the\n\
	half-light radius, and also the 'quartile' radii rql, rqu which\n\
	contain 25 and 75 percent of the light respectively.\n\
	It also outputs a count of 'bad pixels' nbad, which are pixels\n\
	whose centres fall within the aperture, but are either MAGIC or\n\
	lie off the image.\n\
	To compute the magnitude it looks for zeropoint in the input catalogue\n\
	header unless you override this with -z option.\n\
\n\
	It looks for fits filename in the header variable 'fits_name'\n\
	unless you specify alternative with -f option.  The argument here\n\
	can be 'somecommand |' to generate the image on the fly.\n\
\n\
	It will use the local sky background parameters fb0 and dfb\n\
	if present, unless you give the -i option.\n\
	With the '-Z' option, we ignore pixels around other objects if\n\
	distance is <= r * mul.\n\
\n\
	Objects with negative fluxes are assigned magnitude -100.0.\n\
\n\
	In order to get good 'total' magnitudes, it is necessary that you\n\
	use a good radius (big enough to get nearly all the light, but small\n\
	enough to avoid counting neighbour object light).  This is a tricky\n\
	problem in general.  I have found that 'rg' computed by 'hfindpeaks'\n\
	provides a good choice (that's why I adopted it as the default), but\n\
	if you have a catalogue of objects with no decent size parameter\n\
	a reasonable alternative is to run apphot to compute the petrosian\n\
	radius using a fixed aperture of say 10 pixels, and then use a suitable\n\
	multiple (say 2.0) of rp as the aperture for a second pass.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n\n\n"		


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../catlib/cat.h"					/* READ THIS! */
#include "../../imlib/fits.h"
#include "../../utils/error.h"
#include "../../utils/arrays.h"
#include "../../utils/iostream.h"
#include "../../utils/stats_stuff.h"
#include "apphot.h"
#include "zap.h"


main(int argc, char *argv[])	
{
	int		arg = 1, needfitsname, needzeropoint, usesky, usefixedap, dozap, pid, ix, iy;
	double		mul, zapmul;
	char		fitsfile[128], *rname, defaultrname[3] = "rg";
	char		*rzapname, tempfilename[L_tmpnam];
	int		N1, N2;
	fitsheader	*fits;
	FILE		*fitsf, *tempf;
	iostream	*fitsiostream;
	short		**nzap;
	float		**f, **fzap;
	cathead		*inputcathead, *outputcathead;		/* catalogue stuff... */
	object		*inputobject, *outputobject;
	double		*has_sky, *zeropoint, zp, rap;
	item		*fluxitem, *magitem, *rhitem, *rpitem;	/* new items we create here... */
	item		*rqlitem, *rquitem, *nbaditem, *fmaxitem;
	int		iflux, imag, irh, irp, irql, irqu, inbad, ifmax;	/* we'll want their indices... */
	double		*flux, *mag, *rh, *rp, *rql, *rqu, *nbad, *fmax;	/* and local handles */
	double		*x, *r, *fb0, *dfb, *rzap;			/* handles for input variables we want */
	double		rapmax;

	/* defaults */
	mul = 3.0;
	needfitsname = 1;
	needzeropoint = 1;
	rname = defaultrname;
	usesky = 1;
	usefixedap = 0;
	dozap = 0;
	rapmax = 50.0;
	
	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		} else {
			switch (argv[arg++][1]) {
				case 'm':
					sscanf(argv[arg++], "%lf", &mul);
					break;
				case 'r':
					rname = argv[arg++];
					break;
				case 'f':
					strcpy(fitsfile, argv[arg++]);
					needfitsname = 0;
					break;
				case 'i':
					usesky = 0;
					break;
				case 'z':
					sscanf(argv[arg++], "%lf", &zp);
					needzeropoint = 0;
					break;
				case 'Z':
					rzapname = argv[arg++];
					sscanf(argv[arg++], "%lf", &zapmul);
					dozap = 1;
					break;
				case 'M':
					sscanf(argv[arg++], "%lf", &rapmax);
					needzeropoint = 0;
					break;
				case 'R':
					sscanf(argv[arg++], "%lf", &rap);
					usefixedap = 1;
					break;
				default:
					error_exit(usage);
					break;
			}
		}
	}

	setcatopfiletype(BINARY_FILE_TYPE);

	inputcathead = readcathead();			/* read the cat head */
	inputobject = newobject(inputcathead);		/* make the input object */
	connectobjecttocathead(inputobject);		/* obj addresses point back to cathead */
	allocobjectcontents(inputobject);		/* and allocate space for obj data */

	outputcathead = (cathead *) calloc(1, sizeof(cathead));			/* new cathead */
	copyheaderinfo(outputcathead, inputcathead);				/* copy header stuff */

	if (!(inputcathead->headeritembase)) {
		usesky = 0;
		zp = 30.0;
	} else {
		has_sky = (double *) getheaderitemaddress("has_sky", inputcathead);	/* get header field */
		if (!has_sky) {
			usesky = 0;
		} else {
			if  (!(*has_sky)) {
				usesky = 0;
			}
		}
		if (needzeropoint) {
			zeropoint = (double *) getheaderitemaddress("zeropoint", inputcathead);
			if (zeropoint) {
				zp = *zeropoint;
			} else {
				zp = 30.0;
			}
		}
	}
	addargscomment(argc, argv, outputcathead);		/* add history */

	if (needfitsname) {		/* get the fits file name from the header */
		strcpy(fitsfile, *((char **) getheaderitemaddress("fits_name", inputcathead)));
	}
	fitsiostream = openiostream(fitsfile, "r");
	fitsf = fitsiostream->f;
	read2Dfloatimage(&f, &N1, &N2, &fits, fitsf);

	copycontentinfo(outputcathead, inputcathead);		/* copy over pre-exisiting object items */
	fluxitem = newitem("flux", NUM_TYPE, 1, 1);		/* add a new object item */
	addobjectitem(fluxitem, outputcathead);
	magitem = newitem("mag", NUM_TYPE, 1, 1);		/* add another new object item */
	addobjectitem(magitem, outputcathead);
	rhitem = newitem("rh", NUM_TYPE, 1, 1);			/* add another new object item */
	addobjectitem(rhitem, outputcathead);
	rpitem = newitem("rp", NUM_TYPE, 1, 1);			/* add another new object item */
	addobjectitem(rpitem, outputcathead);
	rqlitem = newitem("rql", NUM_TYPE, 1, 1);			/* add another new object item */
	addobjectitem(rqlitem, outputcathead);
	rquitem = newitem("rqu", NUM_TYPE, 1, 1);			/* add another new object item */
	addobjectitem(rquitem, outputcathead);
	nbaditem = newitem("nbad", NUM_TYPE, 1, 1);			/* add another new object item */
	addobjectitem(nbaditem, outputcathead);
	fmaxitem = newitem("fmax", NUM_TYPE, 1, 1);			/* add another new object item */
	addobjectitem(fmaxitem, outputcathead);
	writecathead(outputcathead);				/* and write cathead out */			

	outputobject = newobject(outputcathead);		/* make the output object */
	iflux = getobjectitemindex("flux", outputobject);	/* get indices for new items */
	imag = getobjectitemindex("mag", outputobject);
	irh = getobjectitemindex("rh", outputobject);
	irp = getobjectitemindex("rp", outputobject);
	irql = getobjectitemindex("rql", outputobject);
	irqu = getobjectitemindex("rqu", outputobject);
	inbad = getobjectitemindex("nbad", outputobject);
	ifmax = getobjectitemindex("fmax", outputobject);
	inheritcontents(outputobject, inputobject);		/* pre-existing output object addresses point back to inputobject */
	allocitemcontents(fluxitem, &((outputobject->addrlist)[iflux]), 0);	/* allocate space for new data */
	allocitemcontents(magitem, &((outputobject->addrlist)[imag]), 0);
	allocitemcontents(rhitem, &((outputobject->addrlist)[irh]), 0);
	allocitemcontents(rpitem, &((outputobject->addrlist)[irp]), 0);
	allocitemcontents(rqlitem, &((outputobject->addrlist)[irql]), 0);
	allocitemcontents(rquitem, &((outputobject->addrlist)[irqu]), 0);
	allocitemcontents(nbaditem, &((outputobject->addrlist)[inbad]), 0);
	allocitemcontents(fmaxitem, &((outputobject->addrlist)[ifmax]), 0);
	flux = (double *) ((outputobject->addrlist)[iflux]);	/* get local handles for the new items. */
	mag = (double *) ((outputobject->addrlist)[imag]);
	rh = (double *) ((outputobject->addrlist)[irh]);
	rp = (double *) ((outputobject->addrlist)[irp]);
	rql = (double *) ((outputobject->addrlist)[irql]);
	rqu = (double *) ((outputobject->addrlist)[irqu]);
	nbad = (double *) ((outputobject->addrlist)[inbad]);
	fmax = (double *) ((outputobject->addrlist)[ifmax]);

	/* now we get the handles to the input object items we will need */
	x = (double *) ((inputobject->addrlist)[getobjectitemindex("x", inputobject)]);
	if (!usefixedap) {
		r = (double *) ((inputobject->addrlist)[getobjectitemindex(rname, inputobject)]);
	} else {
		r = NULL;
	}
	if (usesky) {
		fb0 = (double *) ((inputobject->addrlist)[getobjectitemindex("fb0", inputobject)]);
		dfb = (double *) ((inputobject->addrlist)[getobjectitemindex("dfb", inputobject)]);
	} else {
		fb0 = dfb = NULL;
	}


	/* if we want to ignore pixels around other objects we need to create */
	/* fzap, and nzap and spool the cat into a temporary file for re-reading later */
	if (dozap) {
		rzap = (double *) ((inputobject->addrlist)[getobjectitemindex(rzapname, inputobject)]);
		allocFloatArray(&fzap, N1, N2);
		allocShortArray(&nzap, N1, N2);
		for (iy = 0; iy < N2; iy++) {
			for (ix = 0; ix < N1; ix++) {
				fzap[iy][ix] = f[iy][ix];
			}
		}
		tmpnam(tempfilename); 
		if (!(tempf = fopen(tempfilename, "w"))) {
			error_exit("apphot: unable to open temporary file for writing\n");
		}
		setcatopf(tempf);
		while (readobject(inputobject)) {
			zap(ZAP_MODE, *rzap * zapmul, x, f, fzap, nzap, N1, N2);
			writeobject(inputobject);
		}
		fclose(tempf);
		if (!(tempf = fopen(tempfilename, "r"))) {
			error_exit("apphot: unable to open temporary file for reading\n");
		}
		setcatipf(tempf);
		setcatopf(stdout);
	} else {
		fzap = f;
	}


	while (readobject(inputobject)) {					/* big loop */
		if (dozap) {
			zap(UNZAP_MODE, *rzap * zapmul, x, f, fzap, nzap, N1, N2);
		}
		if (!usefixedap) {
			rap = *r * mul;
		}
		if (rap <= 0.0) {
			writeobject(outputobject);
			continue;
		}
		if (rap > rapmax) {
			rap = rapmax;
		}
		if (apphot(flux, rh, rql, rqu, rp, nbad, fmax, x , rap, fb0, dfb, fzap, N1, N2)) {
			if (*flux > 0) {
				*mag = zp - 2.5 * log10(*flux);
			} else {
				*mag = -100.0;
			}
		}
		writeobject(outputobject);
		if (dozap) {
			zap(ZAP_MODE, *rzap * zapmul, x, f, fzap, nzap, N1, N2);
		}
	}
	if (dozap) {
		remove (tempfilename);
	}
	exit(0);
}

