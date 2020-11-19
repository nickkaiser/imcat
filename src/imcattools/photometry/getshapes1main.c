/*
 *	Getshapesmain.c - a heavily documented example of front end for catalogue processing 
 */

#define	usage "\n\n\n\
NAME\n\
	getshapes1 --- calculate ellipticities etc. for catalogue of objects\n\
\n\
SYNOPSIS\n\
	getshapes1	[options...]\n\
		-r rname	# name for window radius ('rg')\n\
		-m rmult	# multiplier: r_window = rmult * rname (1.0)\n\
		-f fitsfile	# specify fits file explicitly\n\
		-i 		# ignore the sky background information.\n\
		-Z r mul	# ignore pixels within mul * r of other objects\n\
		-R		# always output object\n\
\n\
DESCRIPTION\n\
	\"getshapes1\" calculates second moments of sky brightness\n\
	(and related quantities) for objects detected by (h)findpeaks,\n\
	though possibly after having been processed by getsky and/or\n\
	apphot.  It uses a gaussian window of size determined by\n\
	flags -r, -m.  Use '-r unity' for r_window = rmult.\n\
	It adds the following items to the catalogue:\n\
		qll		# trace of flux normalised quadrupole moment matrix\n\
		q[2]		# q[a] = M[a][l][m] q[l][m]\n\
		R[2][2]		# response to psf anisotropy\n\
		P[2][2]		# response to (post seeing) shear\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
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
#include "../../utils/stats_stuff.h"
#include "getshape1.h"
#include "zap.h"

#define Q_DIM	3

main(int argc, char *argv[])	
{
	int		arg = 1, mode, needfitsname, useunitradius, usesky, dozap, ix, iy, rambomode;
	double		rwindow, rmult;
	char		fitsfile[128];	/* fits stuff */
	char		*rzapname, tempfilename[L_tmpnam];
	int		N1, N2;
	fitsheader	*fits;
	FILE		*fitsf, *tempf;
	float		**f, **fzap;
	short		**nzap;
	cathead		*inputcathead, *outputcathead;		/* catalogue stuff... */
	object		*inputobject, *outputobject;
	double		*has_sky;				/* will point to header item */
	item		*qitem, *Ritem, *Pitem;			/* new items we create here... */
	int		qindex, Rindex, Pindex;			/* we'll want their indices... */
	double		*q, **R, **P;				/* and local handles */
	double		*x, *r, *fb0, *dfb, unity[1] = {1.0}, *rzap, zapmul, *flux;		/* handles for input variables we want */
	char		*rname, defrname[3] = "rg";

	/* defaults */
	rname = defrname;
	rmult = 1.0;
	needfitsname = 1;
	usesky = 1;
	useunitradius = 0;
	dozap = 0;
	
	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		} else {
			switch (argv[arg++][1]) {
				case 'r':
					if (strcmp(argv[arg], "unity")) {
						rname = argv[arg++];
					} else {
						useunitradius = 1;
						arg++;
					}
					break;
				case 'm':
					sscanf(argv[arg++], "%lf", &rmult);
					break;
				case 'f':
					sscanf(argv[arg++], "%s", fitsfile);
					needfitsname = 0;
					break;
				case 'i':
					usesky = 0;
					break;
				case 'Z':
					rzapname = argv[arg++];
					sscanf(argv[arg++], "%lf", &zapmul);
					dozap = 1;
					break;
				case 'R':
					rambomode = 1;
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

	if (!(outputcathead->headeritembase)) {
		fprintf(stderr, "# getshapes: warning - no headeritems - usuing no sky correction\n");
		usesky = 0;
	} else {
		has_sky = (double *) getheaderitemaddress("has_sky", outputcathead);
		if (!(*has_sky))
			usesky = 0;
	}
	addargscomment(argc, argv, outputcathead);		/* add history */

	if (needfitsname) {		/* get the fits file name from the header */
		strcpy(fitsfile, *((char **) getheaderitemaddress("fits_name", inputcathead)));
	}
	if (!(fitsf = fopen(fitsfile, "r")))
		error_exit("getsky: can't open fits file\n");
	read2Dfloatimage(&f, &N1, &N2, &fits, fitsf);

	copycontentinfo(outputcathead, inputcathead);		/* copy over pre-exisiting object items */
	qitem = newitem("q", NUM_TYPE, 1, Q_DIM);			/* add a new object item */
	addobjectitem(qitem, outputcathead);
	Ritem = newitem("R", NUM_TYPE, 2, Q_DIM, Q_DIM);			/* add a new object item */
	addobjectitem(Ritem, outputcathead);
	Pitem = newitem("P", NUM_TYPE, 2, Q_DIM, Q_DIM);			/* add a new object item */
	addobjectitem(Pitem, outputcathead);
	writecathead(outputcathead);				/* and write cathead out */			

	outputobject = newobject(outputcathead);		/* make the output object */
	qindex = getobjectitemindex("q", outputobject);
	Rindex = getobjectitemindex("R", outputobject);
	Pindex = getobjectitemindex("P", outputobject);
	inheritcontents(outputobject, inputobject);		/* pre-existing output object addresses point back to inputobject */
	allocitemcontents(qitem, &((outputobject->addrlist)[qindex]), 0);
	allocitemcontents(Ritem, &((outputobject->addrlist)[Rindex]), 0);
	allocitemcontents(Pitem, &((outputobject->addrlist)[Pindex]), 0);
	q = (double *) ((outputobject->addrlist)[qindex]);
	R = (double **) ((outputobject->addrlist)[Rindex]);
	P = (double **) ((outputobject->addrlist)[Pindex]);

	/* now we get the handles to the input object items we will need */
	x = (double *) ((inputobject->addrlist)[getobjectitemindex("x", inputobject)]);
	if (useunitradius) {
		r = unity;
	} else {
	  	r = (double *) ((inputobject->addrlist)[getobjectitemindex(rname, inputobject)]);
	}
	if (usesky) {
	  	fb0 = (double *) ((inputobject->addrlist)[getobjectitemindex("fb0", inputobject)]);
	  	dfb = (double *) ((inputobject->addrlist)[getobjectitemindex("dfb", inputobject)]);
	} else {
		fb0 = dfb = NULL;
	}
	flux = (double *) ((inputobject->addrlist)[getobjectitemindex("flux", inputobject)]);

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
		rwindow = *r * rmult;
		if (dozap)
			zap(UNZAP_MODE, *rzap * zapmul, x, f, fzap, nzap, N1, N2);
		Getshape(x, fb0, dfb, rwindow, flux, q, R, P, fzap, N1, N2);
		writeobject(outputobject);
		if (dozap)
			zap(ZAP_MODE, *rzap * zapmul, x, f, fzap, nzap, N1, N2);
	}
	remove (tempfilename);
	
	exit(0);
}

