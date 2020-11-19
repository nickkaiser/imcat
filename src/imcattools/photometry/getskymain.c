/*
 *	getskymain.c - a heavily documented example of front end for catalogue processing 
 */

#define	usage "\n\n\n\
NAME\n\
	getsky --- determine model for local sky background around objects\n\
\n\
SYNOPSIS\n\
	getsky	[options...]\n\
		-m mode		# mode for background sky calculation (1)\n\
		-a a1 a2	# inner and outer radii (default 16 32)\n\
		-f fitsfile	# specify fits file explicitly\n\
		-Z r mul	# ignore pixels within mul * r of other objects\n\
\n\
DESCRIPTION\n\
	\"getsky\" determines the local sky background for objects detected\n\
	by (h)findpeaks. We determine a mean + gradient model using\n\
	modal sky value for each of four quadrants of an\n\
	annulus around the object.  The inner and outer annulus radii\n\
	are a1, and a2 times:\n\
		1 pixel		# for mode = 1\n\
		rg		# for mode = 2\n\
		rh		# for mode = 3\n\
		rp		# for mode = 4   \n\
	The latter two require that we have already run analyse once.\n\
	Getsky adds entries fb0, dfb[2] to the catalogue and sets the\n\
	header value 'has_sky'.\n\
	By default we use the image named in the catalogue, but you can\n\
	specify an alternative (if you want to use one where the objects\n\
	have been 'zapped' with makechart for instance).\n\
	We require at least 16 good pixels in a viable quadrant.  Our strategy is:\n\
		0 good quadrants:	# mean and gradient = 0.0\n\
		1,2 good quadrants:	# only calculate mean, gradient = 0\n\
		3 good quadrants:	# determine mean and gradient\n\
		4 good quadrants:	# remove most extreme quadrant.   \n\
	With the '-Z' option, we ignore pixels around other objects if\n\
	distance is <= r * mul.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
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
#include "getsky.h"
#include "zap.h"

#define	CONST_ANNULUS_MODE	1
#define RG_ANNULUS_MODE		2
#define RH_ANNULUS_MODE		3
#define RP_ANNULUS_MODE		4

main(int argc, char *argv[])	
{
	int		arg = 1, mode, needfitsname;		/* command line stuff */
	double		a1, a2, r1, r2;
	char		fitsfile[128];				/* fits stuff */
	char		*rzapname, tempfilename[L_tmpnam];
	int		N1, N2;
	fitsheader	*fits;
	FILE		*fitsf, *tempf;
	short		**nzap;
	float		**f, **fzap;
	cathead		*inputcathead, *outputcathead;		/* catalogue stuff... */
	object		*inputobject, *outputobject;
	double		*has_sky;				/* will point to header item */
	item		*fb0item, *dfbitem;			/* new items we create here... */
	int		fb0index, dfbindex;			/* we'll want their indices... */
	double		*fb0, *dfb;				/* and local handles */
	double		*x, *r, unity[1] = {1.0}, *rzap;		/* handles for input variables we want */
	int		dozap, ix, iy;
	double		zapmul;

	/* defaults */
	mode = CONST_ANNULUS_MODE;
	a1 = 16;
	a2 = 32;
	needfitsname = 1;
	dozap = 0;
	
	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		} else {
			switch (argv[arg++][1]) {
				case 'm':
					sscanf(argv[arg++], "%d", &mode);
					break;
				case 'a':
					sscanf(argv[arg++], "%lf", &a1);
					sscanf(argv[arg++], "%lf", &a2);
					break;
				case 'f':
					sscanf(argv[arg++], "%s", fitsfile);
					needfitsname = 0;
					break;
				case 'Z':
					rzapname = argv[arg++];
					sscanf(argv[arg++], "%lf", &zapmul);
					dozap = 1;
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

	has_sky = (double *) getheaderitemaddress("has_sky", outputcathead);	/* modify header field */
	*has_sky = 1;
	addargscomment(argc, argv, outputcathead);		/* add history */

	if (needfitsname) {		/* get the fits file name from the header */
		strcpy(fitsfile, *((char **) getheaderitemaddress("fits_name", inputcathead)));
	}
	if (!(fitsf = fopen(fitsfile, "r")))
		error_exit("getsky: can't open fits file\n");
	read2Dfloatimage(&f, &N1, &N2, &fits, fitsf);

	copycontentinfo(outputcathead, inputcathead);		/* copy over pre-exisiting object items */
	fb0item = newitem("fb0", NUM_TYPE, 1, 1);		/* add a new object item */
	addobjectitem(fb0item, outputcathead);
	dfbitem = newitem("dfb", NUM_TYPE, 1, 2);		/* add another new object item */
	addobjectitem(dfbitem, outputcathead);
	writecathead(outputcathead);				/* and write cathead out */			

	outputobject = newobject(outputcathead);		/* make the output object */
	fb0index = getobjectitemindex("fb0", outputobject);	/* get indices for new items */
	dfbindex = getobjectitemindex("dfb", outputobject);
	inheritcontents(outputobject, inputobject);		/* pre-existing output object addresses point back to inputobject */
	allocitemcontents(fb0item, &((outputobject->addrlist)[fb0index]), 0);	/* allocate space for new data */
	allocitemcontents(dfbitem, &((outputobject->addrlist)[dfbindex]), 0);
	fb0 = (double *) ((outputobject->addrlist)[fb0index]);	/* get local handles for the new items. */
	dfb = (double *) ((outputobject->addrlist)[dfbindex]);	/* Note that both 'scalar' fb0 and 'vector' dfb are */
								/* both 'double *'.  A matrix would be 'double **' though. */
								/* We can now assign to '*fb0' and to dfb[0], dfb[1]. */

	/* now we get the handles to the input object items we will need */
	x = (double *) ((inputobject->addrlist)[getobjectitemindex("x", inputobject)]);
	switch (mode) {
	 case CONST_ANNULUS_MODE:
	  r = unity;
	  break;
	 case RG_ANNULUS_MODE:
	  r = (double *) ((inputobject->addrlist)[getobjectitemindex("rg", inputobject)]);
 	  break;
	 case RH_ANNULUS_MODE:
	  r = (double *) ((inputobject->addrlist)[getobjectitemindex("rh", inputobject)]);
 	  break;
	 case RP_ANNULUS_MODE:
	  r = (double *) ((inputobject->addrlist)[getobjectitemindex("rp", inputobject)]);
 	  break;
	 default:
	  error_exit("getsky: bad mode\n");
	  break;
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
		r1 = *r * a1;
		r2 = *r * a2;
		if (dozap)
			zap(UNZAP_MODE, *rzap * zapmul, x, f, fzap, nzap, N1, N2);
		getsky(fb0, dfb, r1, r2, x, fzap, N1, N2);			/* do the biz...*/
		writeobject(outputobject);
		if (dozap)
			zap(ZAP_MODE, *rzap * zapmul, x, f, fzap, nzap, N1, N2);
	}
	remove (tempfilename);	
	exit(0);
}

