/*
 *	getshapes2main.c - new version shape measurer 
 */

#define	usage "\n\n\n\
NAME\n\
	getshapes2 --- calculate ellipticities etc. for catalogue of objects\n\
\n\
SYNOPSIS\n\
	getshapes2 psfimage rf [options...]\n\
		-f fitsfile	# specify fits file explicitly\n\
		-i 		# ignore the sky background information.\n\
		-Z r mul	# ignore pixels within mul * r of other objects\n\
\n\
DESCRIPTION\n\
	\"getshapes2\" is a catalogue filter which calculates polarisation\n\
	and polarisability for objects.\n\
\n\
	It requires that the input catalogue contain at least position\n\
	vector 'x[2]', and a (aperture) flux 'flux'.\n\
\n\
	It uses 'makekernel' to compute the kernels W[i] and K[i][j]\n\
	from the psf supplied in the fits image 'psfimage' and\n\
	then computes\n\
		F 	= sum f[y][x] w[y][x]\n\
		q[l] 	= sum f[y][x] W[l][y][x] / F\n\
		P[l][m] = sum f[y][x] K[l][m][y][x] / F\n\
		R[m]	= sum f[y][x] R[m][y][x] / F\n\
	for l = 0,1,2 and m = 1,2\n\
\n\
	If sky background values (from getsky) are present they will\n\
	be used (unless you give -i flag).\n\
\n\
	Use -Z option to zap circles around neighbouring objects.\n\
\n\
	By default the source image name is taken from the catalogue\n\
	header item 'fits_name', but you can specify alternative\n\
	explicitly with -f option.\n\
\n\
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
#include "getshape2.h"
#include "zap.h"

#define MAGIC FLOAT_MAGIC

main(int argc, char *argv[])	
{
	int		arg = 1, needfitsname, usesky, dozap, ix, iy, l, m;
	char		fitsfilename[128], *psffilename;
	char		*rzapname, tempfilename[L_tmpnam], tmpstring[1024];
	int		N1, N2, M1, M2;
	fitsheader	*fits, *kernel;
	FILE		*fitsf, *tempf, *kernelf;
	float		**f, **fzap, **w, ***W, ***RR, ****K;
	short		**nzap;
	cathead		*inputcathead, *outputcathead;		/* catalogue stuff... */
	object		*inputobject, *outputobject;
	double		*has_sky;				/* will point to header item */
	item		*Fitem, *qitem, *Pitem, *Ritem;			/* new items we create here... */
	int		Findex, qindex, Pindex, Rindex;			/* we'll want their indices... */
	double		*q, **P, *F, *R;				/* and local handles */
	double		*x, *fb0, *dfb, *rzap, zapmul, *flux;		/* handles for input variables we want */
	double		rf;

	/* defaults */
	needfitsname = 1;
	usesky = 1;
	dozap = 0;
	
	if (argc < 3) {
		error_exit(usage);
	}
	psffilename = argv[arg++];
	if (1 != sscanf(argv[arg++], "%lf", &rf)) {
		error_exit(usage);
	}

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		} else {
			switch (argv[arg++][1]) {
				case 'f':
					strcpy(fitsfilename, argv[arg++]);
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
		strcpy(fitsfilename, *((char **) getheaderitemaddress("fits_name", inputcathead)));
	}
	if (!(fitsf = fopen(fitsfilename, "r")))
		error_exit("getsky: can't open fits file\n");
	read2Dfloatimage(&f, &N1, &N2, &fits, fitsf);

	/* read the kernels */
	W = (float ***) calloc(3, sizeof(float **));
	RR = (float ***) calloc(2, sizeof(float **));
	K = (float ****) calloc(3, sizeof(float ***));
	K[0] = (float ***) calloc(2, sizeof(float **));
	K[1] = (float ***) calloc(2, sizeof(float **));
	K[2] = (float ***) calloc(2, sizeof(float **));
	sprintf(tmpstring, "test -f %s", psffilename);
	if (system(tmpstring)) {
		error_exit("getshapes2: psf file does not exist; bailing out!\n");
	}
	sprintf(tmpstring, "makekernel -w %lf < %s", rf, psffilename);
	kernelf = (FILE *) popen(tmpstring, "r");
	if (!kernelf) {
		error_exit("getshapes2: makekernel failed\n");
	}
	read2Dfloatimage(&w, &M1, &M2, &kernel, kernelf);
	pclose(kernelf);
	for (l = 0; l < 3; l++) {
		sprintf(tmpstring, "makekernel -W %d %lf < %s", l, rf, psffilename);
		kernelf = (FILE *) popen(tmpstring, "r");
		if (!kernelf) {
			error_exit("getshapes2: makekernel failed\n");
		}
		read2Dfloatimage(&(W[l]), &M1, &M2, &kernel, kernelf);
		pclose(kernelf);
		if (l < 2) {
			sprintf(tmpstring, "makekernel -R %d %lf < %s", l + 1, rf, psffilename);
			kernelf = (FILE *) popen(tmpstring, "r");
			if (!kernelf) {
				error_exit("getshapes2: makekernel failed\n");
			}
			read2Dfloatimage(&(RR[l]), &M1, &M2, &kernel, kernelf);
			pclose(kernelf);
		}
		for (m = 0; m < 2; m++) {
			sprintf(tmpstring, "test -f %s", psffilename);
			if (system(tmpstring)) {
				error_exit("getshapes2: psf file does not exist; bailing!\n");
			}
			sprintf(tmpstring, "makekernel -K %d %d %lf < %s", l, m + 1, rf, psffilename);
			kernelf = (FILE *) popen(tmpstring, "r");
			if (!kernelf) {
				error_exit("getshapes2: makekernel failed\n");
			}
			read2Dfloatimage(&(K[l][m]), &M1, &M2, &kernel, kernelf);
			pclose(kernelf);
		}
	}

	copycontentinfo(outputcathead, inputcathead);		/* copy over pre-exisiting object items */
	Fitem = newitem("F", NUM_TYPE, 1, 1);			/* add a new object item */
	addobjectitem(Fitem, outputcathead);
	qitem = newitem("q", NUM_TYPE, 1, 3);			/* add a new object item */
	addobjectitem(qitem, outputcathead);
	Pitem = newitem("P", NUM_TYPE, 2, 3, 2);			/* add a new object item */
	addobjectitem(Pitem, outputcathead);
	Ritem = newitem("R", NUM_TYPE, 1, 2);			/* add a new object item */
	addobjectitem(Ritem, outputcathead);
	writecathead(outputcathead);				/* and write cathead out */			

	outputobject = newobject(outputcathead);		/* make the output object */
	Findex = getobjectitemindex("F", outputobject);	/* get indices for new items */
	qindex = getobjectitemindex("q", outputobject);
	Pindex = getobjectitemindex("P", outputobject);
	Rindex = getobjectitemindex("R", outputobject);
	inheritcontents(outputobject, inputobject);		/* pre-existing output object addresses point back to inputobject */
	allocitemcontents(Fitem, &((outputobject->addrlist)[Findex]), 0);	/* allocate space for new data */
	allocitemcontents(qitem, &((outputobject->addrlist)[qindex]), 0);
	allocitemcontents(Pitem, &((outputobject->addrlist)[Pindex]), 0);
	allocitemcontents(Ritem, &((outputobject->addrlist)[Rindex]), 0);
	F = (double *) ((outputobject->addrlist)[Findex]);	/* get local handles for the new items. */
	q = (double *) ((outputobject->addrlist)[qindex]);
	P = (double **) ((outputobject->addrlist)[Pindex]);
	R = (double *) ((outputobject->addrlist)[Rindex]);

	/* now we get the handles to the input object items we will need */
	x = (double *) ((inputobject->addrlist)[getobjectitemindex("x", inputobject)]);
	if (usesky) {
	  	fb0 = (double *) ((inputobject->addrlist)[getobjectitemindex("fb0", inputobject)]);
	  	dfb = (double *) ((inputobject->addrlist)[getobjectitemindex("dfb", inputobject)]);
	} else {
		fb0 = dfb = NULL;
	}
	/* flux = (double *) ((inputobject->addrlist)[getobjectitemindex("flux", inputobject)]); */

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
		if (dozap)
			zap(UNZAP_MODE, *rzap * zapmul, x, f, fzap, nzap, N1, N2);
		getshape2(x, fb0, dfb, F, q, P, R, fzap, N1, N2, w, W, RR, K, M1, M2);
		writeobject(outputobject);
		if (dozap)
			zap(ZAP_MODE, *rzap * zapmul, x, f, fzap, nzap, N1, N2);
	}
	remove (tempfilename);
	
	exit(0);
}

