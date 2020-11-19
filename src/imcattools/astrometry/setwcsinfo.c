/*
 * getwcsinfo.c
 */

#define usage "\n\
NAME\n\
	setwcsinfo --- insert world coordinate system info into fits header\n\
\n\
SYNOPSIS\n\
	setwcsinfo CTYPE1 CDELT1 CRVAL1 CRPIX1 CTYPE2 CDELT2 CRVAL2 CRPIX2 CROTA2 [-c | -p] [-e EPOCH EQUINOX]\n\
\n\
DESCRIPTION\n\
	setwcsinfo reads a fits image from stdin, removes\n\
	any existing WCS info; inserts new ones; and writes\n\
	out the data.\n\
\n\
	Rotation angle CROTA2 must be given in degrees.\n\
\n\
	With -c option we encode the rotation and scaling in the CDn_m matrix.  However,\n\
	this practice is deprecated.\n\
\n\
	With -p option we encode the rotation angle in the PCnnnmmm matrix.\n\
\n\
	Refer to Greisen and Calabretta, 96 for description\n\
	of FITS WCS convention.\n\
\n\
AUTHOR\n\
	Nick Kaiser	kaiser@hawaii.edu\n\
\n\n"


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


#include "utils/error.h"
#include "utils/args.h"
#include "imlib/fits.h"

#define	NNAMES	32

#define CROTA2_ENCODING	0
#define CD_ENCODING	1
#define PC_ENCODING	2

main (int argc, char *argv[])
{
	fitsheader	*fits;
	fitscomment	*com;
	char		*name[NNAMES] = {
				"RA", 
				"DEC", 
				"EQUINOX", 
				"EPOCH", 
				"CRPIX1", 
				"CRVAL1", 
				"CDELT1",
				"CTYPE1",
				"CRPIX2", 
				"CRVAL2", 
				"CDELT2",
				"CTYPE2",
				"PLTRAH",
				"PLTRAM",
				"PLTRAS",
				"PLTDECH",
				"PLTDECM",
				"PLTDECS",
				"OBJCTRA",
				"OBJCTDEC",
				"CROTA2",
				"PC001001",
				"PC001002",
				"PC002001",
				"PC002002",
				"CD1_1",
				"CD1_2",
				"CD2_1",
				"CD2_2",
				"LIM1_1",
				"LIM2_2",
				"",
			}, string[9], *epoch, *equinox;
	int		i, rotencoding;
	void		*f;
	double		c, s, theta, cd11, cd12, cd21, cd22, cdelt1, cdelt2, S;
	char		*CTYPE1, *CDELT1, *CRVAL1, *CRPIX1, *CTYPE2, *CDELT2, *CRVAL2, *CRPIX2, *CROTA2, *EPOCH, *EQUINOX;
	char		*flag;

	epoch   = "2000";
	equinox = "2000";
	rotencoding = CROTA2_ENCODING;

	argsinit(argc, argv, usage);
	if (nextargtype() == FLAG_ARG) {
		error_exit(usage);
	}
	CTYPE1 = getargs();
	CDELT1 = getargs();
	CRVAL1 = getargs();
	CRPIX1 = getargs();
	CTYPE2 = getargs();
	CDELT2 = getargs();
	CRVAL2 = getargs();
	CRPIX2 = getargs();
	CROTA2 = getargs();
	while (flag = getflag()) {
		 switch (flag[0]) {
			case 'c':
				rotencoding = CD_ENCODING;
				break;
			case 'p':
				rotencoding = PC_ENCODING;
				break;
			case 'e':
				epoch = getargs();
				equinox = getargs();
				break;
			default:
				error_exit(usage);
		}
	}
	
	fits = readfitsheader(stdin);
	for (i = 0; i < NNAMES; i++) {
		removenamedcomments(name[i], fits);
		sprintf(string, "PPO%d", i);
		removenamedcomments(string, fits);
		sprintf(string, "AMDX%d", i);
		removenamedcomments(string, fits);
		sprintf(string, "AMDY%d", i);
		removenamedcomments(string, fits);
	}
	for (i = 0; i < 10; i++) {
		sprintf(string, "WAT%d_001", i);
		removenamedcomments(string, fits);
	}

	if (1 != sscanf(argv[9], "%lf", &theta)) {
		error_exit("setwcsinfo: failed to decipher CROTA2 arg\n");
	}
	c = cos(M_PI * theta / 180.0);
	s = sin(M_PI * theta / 180.0);

	sscanf(CDELT1, "%lf", &cdelt1);
	sscanf(CDELT2, "%lf", &cdelt2);
	S = cdelt2 / cdelt1;

	switch (rotencoding) {
		case CROTA2_ENCODING:
			com = newnumericcomment("CROTA2", theta, "WCS --- image rotation");
			prependcomment(com, fits);
			break;
		case PC_ENCODING:
			com = newnumericcomment("PC002002",  c, "WCS --- PC matrix");
			prependcomment(com, fits);
			com = newnumericcomment("PC002001",  s / S, "WCS --- PC matrix");
			prependcomment(com, fits);
			com = newnumericcomment("PC001002", -s * S, "WCS --- PC matrix");
			prependcomment(com, fits);
			com = newnumericcomment("PC001001",  c, "WCS --- PC matrix");
			prependcomment(com, fits);
			break;
		case CD_ENCODING:
			cd11 = c * cdelt1;
			cd12 = -s * cdelt2;
			cd21 = s * cdelt1;
			cd22 = c * cdelt2;
			com = newnumericcomment("CD2_2",  cd22, "WCS --- CD matrix");
			prependcomment(com, fits);
			com = newnumericcomment("CD2_1",  cd21, "WCS --- CD matrix");
			prependcomment(com, fits);
			com = newnumericcomment("CD1_2", cd12, "WCS --- CD matrix");
			prependcomment(com, fits);
			com = newnumericcomment("CD1_1",  cd11, "WCS --- CD matrix");
			prependcomment(com, fits);
			break;
	}
	com = newtextcomment("CRPIX2", CRPIX2, "WCS --- y pixel at reference point");
	prependcomment(com, fits);
	com = newtextcomment("CRVAL2", CRVAL2, "WCS --- DEC at reference point");
	prependcomment(com, fits);
	if (rotencoding != CD_ENCODING) {
		com = newtextcomment("CDELT2", CDELT2, "WCS --- y increment per pixel (degrees)");
		prependcomment(com, fits);
	}
	com = newtextcomment("CTYPE2", CTYPE2, "WCS --- y coordinate type");
	prependcomment(com, fits);
	com = newtextcomment("CRPIX1", CRPIX1, "WCS --- x pixel at reference point");
	prependcomment(com, fits);
	com = newtextcomment("CRVAL1", CRVAL1, "WCS --- RA at reference point");
	prependcomment(com, fits);
	if (rotencoding != CD_ENCODING) {
		com = newtextcomment("CDELT1", CDELT1, "WCS --- x increment per pixel (degrees)");
		prependcomment(com, fits);
	}
	com = newtextcomment("CTYPE1", CTYPE1, "WCS --- x coordinate type");
	prependcomment(com, fits);
	com = newtextcomment("EPOCH", epoch, "WCS --- obs epoch");
	prependcomment(com, fits);
	com = newtextcomment("EQUINOX", equinox, "WCS --- RA/DEC equinox");
	prependcomment(com, fits);

	fits->intpixtype = fits->extpixtype;	
	f = calloc(fits->n[0] * pixsize(fits->intpixtype), sizeof(char));

	writefitsheader(fits);
	for (i = 0; i < fits->n[1]; i++) {
		readfitsline(f, fits);
		writefitsline(f, fits);
	}
	writefitstail(fits);
	exit(0);
}
