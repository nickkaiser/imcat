/*
 *	getshapes3main.c - new version shape measurer 
 */

#define	usage "\n\n\n\
NAME\n\
	getshapes3 --- calculate ellipticities etc. for catalogue of objects\n\
\n\
SYNOPSIS\n\
	getshapes3 f_image fs_image psf_image rf [options...]\n\
		-i 		# ignore the sky background information.\n\
\n\
DESCRIPTION\n\
	\"getshapes3\" is a catalogue filter which calculates polarisation\n\
	and polarisability for objects.\n\
\n\
	It requires that the input catalogue contain at least position\n\
	vector 'x[2]', and a (aperture) flux 'flux'.\n\
\n\
	It uses 'makekernel' to compute the kernels W[i], K[i][j], ZZ[i][j]\n\
	from the psf supplied in the fits image 'psf_image' and\n\
	then computes\n\
\n\
		F 	= sum fs[y][x] w[y][x]\t\n\
		q0 	= sum fs[y][x] W[0][y][x] / F\t\n\
		q[l] 	= sum fs[y][x] W[1+l][y][x] / F\t\n\
		P0[m] = sum f[y][x] K[0][1+m][y][x] / F\t\n\
		P[l][m] = sum f[y][x] K[1+l][1+m][y][x] / F\t\n\
		Z[l][m] = sum f[y][x] ZZ[1+l][1+m][y][x] / F\t\n\
		R[m]	= sum f[y][x] R[m][y][x] / F\t\n\
		Psh[l][m] = 2 q0 delta[l][m] - Z[l][m] / 2 rf^2 + 2 q[l] q[m] / rf^2\t\n\
		Psh[l][m] = ((rf^2 - 2 q0) delta[l][m] + Z[l][m] / 4 rf^2 - q[l] q[m] / rf^2)\t\n\
\n\
	for l,m = 0,1\n\
\n\
	If sky background values (from getsky) are present they will\n\
	be used (unless you give -i flag).\n\
\n\
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
#include "getshape3.h"
#include "zap.h"

#define MAGIC FLOAT_MAGIC

main(int argc, char *argv[])	
{
	int		arg = 1, needfitsname, usesky, dozap, ix, iy, l, m;
	char		*ffilename, *fsfilename, *psffilename, tmpstring[256];
	int		N1, N2, M1, M2;
	fitsheader	*ffits, *fsfits, *kernel;
	FILE		*ffitsf, *fsfitsf, *tempf, *kernelf;
	float		**f, **fs, **w, ***W, ***RR, ****K, ****ZZ;
	cathead		*inputcathead, *outputcathead;		/* catalogue stuff... */
	object		*inputobject, *outputobject;
	double		*has_sky;				/* will point to header item */
	item		*Fitem, *q0item, *qitem, *P0item, *Pitem, *Zitem, *Ritem, *Pshitem, *Psmitem;	/* new items we create here... */
	int		Findex, q0index, qindex, P0index, Pindex, Zindex, Rindex, Pshindex, Psmindex;	/* we'll want their indices... */
	double		*q0, *q, *P0, **P, **Z, *F, *R, **Psm, **Psh;	/* and local handles */
	double		*x, *fb0, *dfb, *rzap, zapmul, *flux;		/* handles for input variables we want */
	double		rf;

	/* defaults */
	usesky = 1;
	dozap = 0;
	
	if (argc < 5) {
		error_exit(usage);
	}
	ffilename = argv[arg++];
	fsfilename = argv[arg++];
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
				case 'i':
					usesky = 0;
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

	if (!(ffitsf = fopen(ffilename, "r")))
		error_exit("getshapes3: can't open f_fits file\n");
	read2Dfloatimage(&f, &N1, &N2, &ffits, ffitsf);
	if (!(fsfitsf = fopen(fsfilename, "r")))
                error_exit("getshapes3: can't open fs_fits file\n");
        read2Dfloatimage(&fs, &N1, &N2, &fsfits, fsfitsf);


	/* read the kernels */
	W = (float ***) calloc(3, sizeof(float **));
	RR = (float ***) calloc(2, sizeof(float **));
	K = (float ****) calloc(3, sizeof(float ***));
	K[0] = (float ***) calloc(2, sizeof(float **));
	K[1] = (float ***) calloc(2, sizeof(float **));
	K[2] = (float ***) calloc(2, sizeof(float **));
	ZZ = (float ****) calloc(3, sizeof(float ***));
	ZZ[0] = (float ***) calloc(2, sizeof(float **));
	ZZ[1] = (float ***) calloc(2, sizeof(float **));
	ZZ[2] = (float ***) calloc(2, sizeof(float **));
	sprintf(tmpstring, "test -f %s", psffilename);
	if (system(tmpstring)) {
		error_exit("getshapes3: psf file does not exist; bailing out!\n");
	}
	sprintf(tmpstring, "makekernel -w %lf < %s", rf, psffilename);
	kernelf = (FILE *) popen(tmpstring, "r");
	if (!kernelf) {
		error_exit("getshapes3: makekernel failed\n");
	}
	read2Dfloatimage(&w, &M1, &M2, &kernel, kernelf);
	pclose(kernelf);
	for (l = 0; l < 3; l++) {
		sprintf(tmpstring, "makekernel -W %d %lf < %s", l, rf, psffilename);
		kernelf = (FILE *) popen(tmpstring, "r");
		if (!kernelf) {
			error_exit("getshapes3: makekernel failed\n");
		}
		read2Dfloatimage(&(W[l]), &M1, &M2, &kernel, kernelf);
		pclose(kernelf);
		if (l < 2) {
			sprintf(tmpstring, "makekernel -R %d %lf < %s", l + 1, rf, psffilename);
			kernelf = (FILE *) popen(tmpstring, "r");
			if (!kernelf) {
				error_exit("getshapes3: makekernel failed\n");
			}
			read2Dfloatimage(&(RR[l]), &M1, &M2, &kernel, kernelf);
			pclose(kernelf);
		}
		for (m = 0; m < 2; m++) {
			sprintf(tmpstring, "test -f %s", psffilename);
			if (system(tmpstring)) {
				error_exit("getshapes3: psf file does not exist; bailing!\n");
			}
			sprintf(tmpstring, "makekernel -K %d %d %lf < %s", l, m + 1, rf, psffilename);
			kernelf = (FILE *) popen(tmpstring, "r");
			if (!kernelf) {
				error_exit("getshapes3: makekernel failed\n");
			}
			read2Dfloatimage(&(K[l][m]), &M1, &M2, &kernel, kernelf);
			pclose(kernelf);
			sprintf(tmpstring, "makekernel -Z %d %d %lf < %s", l, m + 1, rf, psffilename);
			kernelf = (FILE *) popen(tmpstring, "r");
			if (!kernelf) {
				error_exit("getshapes3: makekernel failed\n");
			}
			read2Dfloatimage(&(ZZ[l][m]), &M1, &M2, &kernel, kernelf);
			pclose(kernelf);
		}
	}

	copycontentinfo(outputcathead, inputcathead);		/* copy over pre-exisiting object items */
	Fitem = newitem("F", NUM_TYPE, 1, 1);			/* add a new object item */
	addobjectitem(Fitem, outputcathead);
	q0item = newitem("q0", NUM_TYPE, 1, 1);			/* add a new object item */
	addobjectitem(q0item, outputcathead);
	qitem = newitem("q", NUM_TYPE, 1, 2);			/* add a new object item */
	addobjectitem(qitem, outputcathead);
	P0item = newitem("P0", NUM_TYPE, 1, 2);			/* add a new object item */
	addobjectitem(P0item, outputcathead);
	Pitem = newitem("P", NUM_TYPE, 2, 2, 2);			/* add a new object item */
	addobjectitem(Pitem, outputcathead);
	Zitem = newitem("Z", NUM_TYPE, 2, 2, 2);			/* add a new object item */
	addobjectitem(Zitem, outputcathead);
	Ritem = newitem("R", NUM_TYPE, 1, 2);			/* add a new object item */
	addobjectitem(Ritem, outputcathead);
	Pshitem = newitem("Psh", NUM_TYPE, 2, 2, 2);			/* add a new object item */
	addobjectitem(Pshitem, outputcathead);
	Psmitem = newitem("Psm", NUM_TYPE, 2, 2, 2);			/* add a new object item */
	addobjectitem(Psmitem, outputcathead);
	writecathead(outputcathead);				/* and write cathead out */			

	outputobject = newobject(outputcathead);		/* make the output object */
	Findex = getobjectitemindex("F", outputobject);	/* get indices for new items */
	q0index = getobjectitemindex("q0", outputobject);
	qindex = getobjectitemindex("q", outputobject);
	P0index = getobjectitemindex("P0", outputobject);
	Pindex = getobjectitemindex("P", outputobject);
	Zindex = getobjectitemindex("Z", outputobject);
	Rindex = getobjectitemindex("R", outputobject);
	Pshindex = getobjectitemindex("Psh", outputobject);
	Psmindex = getobjectitemindex("Psm", outputobject);
	inheritcontents(outputobject, inputobject);		/* pre-existing output object addresses point back to inputobject */
	allocitemcontents(Fitem, &((outputobject->addrlist)[Findex]), 0);	/* allocate space for new data */
	allocitemcontents(q0item, &((outputobject->addrlist)[q0index]), 0);
	allocitemcontents(qitem, &((outputobject->addrlist)[qindex]), 0);
	allocitemcontents(P0item, &((outputobject->addrlist)[P0index]), 0);
	allocitemcontents(Pitem, &((outputobject->addrlist)[Pindex]), 0);
	allocitemcontents(Zitem, &((outputobject->addrlist)[Zindex]), 0);
	allocitemcontents(Ritem, &((outputobject->addrlist)[Rindex]), 0);
	allocitemcontents(Pshitem, &((outputobject->addrlist)[Pshindex]), 0);
	allocitemcontents(Psmitem, &((outputobject->addrlist)[Psmindex]), 0);
	F = (double *) ((outputobject->addrlist)[Findex]);	/* get local handles for the new items. */
	q0 = (double *) ((outputobject->addrlist)[q0index]);
	q = (double *) ((outputobject->addrlist)[qindex]);
	P0 = (double *) ((outputobject->addrlist)[P0index]);
	P = (double **) ((outputobject->addrlist)[Pindex]);
	Z = (double **) ((outputobject->addrlist)[Zindex]);
	R = (double *) ((outputobject->addrlist)[Rindex]);
	Psh = (double **) ((outputobject->addrlist)[Pshindex]);
	Psm = (double **) ((outputobject->addrlist)[Psmindex]);

	/* now we get the handles to the input object items we will need */
	x = (double *) ((inputobject->addrlist)[getobjectitemindex("x", inputobject)]);
	if (usesky) {
	  	fb0 = (double *) ((inputobject->addrlist)[getobjectitemindex("fb0", inputobject)]);
	  	dfb = (double *) ((inputobject->addrlist)[getobjectitemindex("dfb", inputobject)]);
	} else {
		fb0 = dfb = NULL;
	}
	/* flux = (double *) ((inputobject->addrlist)[getobjectitemindex("flux", inputobject)]); */

	while (readobject(inputobject)) {					/* big loop */
		getshape3(x, fb0, dfb, F, q0, q, P0, P, Z, R, f, fs, N1, N2, w, W, RR, K, ZZ, M1, M2);
		/* compute Psm, Psh */
		for (l = 0; l < 2; l++) {
			for (m = 0; m < 2; m++) {
				Psh[l][m] = (2 * q[l] * q[m] - 0.5 * Z[l][m]) / (rf * rf);
				Psm[l][m] = (0.25 * Z[l][m] - q[l] * q[m]) / (rf * rf);
			}
		}
		for (l = 0; l < 2; l++) {
			Psh[l][l] += 2 * (*q0);
			Psm[l][l] += (rf * rf - 2 * (*q0));
		}
		writeobject(outputobject);
	}
	
	exit(0);
}

