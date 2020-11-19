/*
 * getworldcoords.c
 */

#define usage "\n\
NAME\n\
	getworldcoords --- add world coordinates to a catalogue\n\
\n\
SYNOPSIS\n\
	getworldcoords fitsfile\n\
\n\
DESCRIPTION\n\
	getworldcoords reads a catalogue from stdin and computes\n\
	from the pixel coordinates x[] the celestial coordinates\n\
	RA and DEC according to the transformation model encoded\n\
	in 'fitsfile' (which need only be a header).\n\
	RA and DEC are expressed in decimal degrees.\n\
\n\
	Uses Doug Mink's wcs utilities code.\n\
\n\
AUTHOR\n\
	Nick Kaiser	kaiser@hawaii.edu\n\
\n\n"


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


#include "utils/error.h"
#include "catlib/cat.h"
#include "imlib/fits.h"
#include "wcs.h"

main (int argc, char *argv[])
{
	double		*RA, *DEC, *x, xpix, ypix, PC[2][2], PCdet, PCtol;
	cathead		*ipcat, *opcat;
	object		*ipobj, *opobj;
	item		*xitem, *RAitem, *DECitem;
	int		xindex, RAindex, DECindex;
	struct WorldCoor 	*wcs;
	char 		header[100*2880], comment[1024], tmpstring[1024];
	FILE		*headerf;
	int		i, cattype, hasPC, l, m;
	double		rah, ram, ras, decd, decm, decs;


	if (argc != 2) {
		error_exit(usage);
	}
	if (!strcmp("-u", argv[1])) {
		error_exit(usage);
	}

	/* read the header */
	headerf = fopen(argv[1], "r");
	if (!headerf) {
		error_exit("getworldcoords: cannot open fits header\n");
	}
	for(i = 0; i < 3600; i++) {
		fgets(header + COM_LENGTH * i, COM_LENGTH1, headerf);
		if (!strncmp(header + COM_LENGTH * i, "END", 3)) break;
	}
	close(headerf);

	/* get the PC matrix if present */
	hasPC = 0;
	for (l = 0; l < 2; l++) {
		for (m = 0; m < 2; m++) {
			sprintf(tmpstring, "imhead -v PC00%1d00%1d < %s 2> /dev/null", l + 1, m + 1, argv[1]);
			headerf = popen(tmpstring, "r");
			if (1 == fscanf(headerf, "%lf", &(PC[l][m]))) {
				hasPC++;
			}
		}
	}

	wcs = wcsinit(header);
	if (nowcs(wcs)) {
		error_exit("getworldcoords: header does not contain wcs info\n");
	}
	if (hasPC) {
		/* check that PC matrix looks reasonable
		PCtol = 1.e-6;
		PCdet = PC[0][0] * PC[1][1] - PC[0][1] * PC[1][0];
		if (fabs(fabs(PCdet) - 1) > PCtol || fabs(PC[0][0] - PC[1][1]) > PCtol ||  fabs(PC[0][1] + PC[1][0]) > PCtol) {
			error_exit("getworldcoords: PC matrix does not seem to be a simple rotation matrix\n");
		} */
		wcs->rot = atan2(PC[1][0] * wcs->yinc / wcs->xinc, PC[0][0]);
		wcs->crot = cos(wcs->rot);
		wcs->srot = sin(wcs->rot);
		wcs->rot *= 180.0 / M_PI;
	}

	ipcat = readcathead();
	getcatipfiletype(&cattype);
	ipobj = newobject(ipcat);
	connectobjecttocathead(ipobj);
	allocobjectcontents(ipobj);
	opcat = (cathead *) calloc(1, sizeof(cathead));
	copyheaderinfo(opcat, ipcat);
	setcatopfiletype(cattype);
	copycontentinfo(opcat, ipcat);

	addargscomment(argc, argv, opcat);
	sprintf(comment, "WCS epoch:   %13.8lf", wcs->epoch);
	addcomment(comment, opcat); 
	sprintf(comment, "WCS equinox: %6.1lf", wcs->equinox);
	addcomment(comment, opcat); 
	sprintf(comment, "WCS system:  %s", wcs->sysout);
	addcomment(comment, opcat); 

	RAitem = newitem("RA", NUM_TYPE, 1, 1);
	addobjectitem(RAitem, opcat);
	DECitem = newitem("DEC", NUM_TYPE, 1, 1);
	addobjectitem(DECitem, opcat);

	opobj = newobject(opcat);
	RAindex = getobjectitemindex("RA", opobj);
	DECindex = getobjectitemindex("DEC", opobj);
	inheritcontents(opobj, ipobj);
	allocitemcontents(RAitem, &((opobj->addrlist)[RAindex]), 0);
	allocitemcontents(DECitem, &((opobj->addrlist)[DECindex]), 0);
	RA = (double *) ((opobj->addrlist)[RAindex]);
	DEC = (double *) ((opobj->addrlist)[DECindex]);
	xindex = getobjectitemindex("x", ipobj);
	x = (double *) ((ipobj->addrlist)[xindex]);
	
	writecathead(opcat);
	while (readobject(ipobj)) {
		/* add 1/2 because of dumb fits convention */
		xpix = x[0] + 0.5;
		ypix = x[1] + 0.5;
		pix2wcs (wcs, xpix, ypix, RA, DEC);		
		writeobject(opobj);
	}

	exit(0);
}

