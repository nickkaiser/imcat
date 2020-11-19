/*
 * lensmap.c
 */

#define	usage "\n\n\n\
NAME\n\
	mapbynumericdef --- map image using numerical deflection func\n\
\n\
SYNOPSIS\n\
	mapbynumericdef [options...] sourcefits deflectionfits\n\
		-i L1 L2	# perform inverse mapping\n\
		-u		# print this message\n\
		-m mapmode	# type of mapping\n\
		-M		# initialise image to magic\n\
		-s sfac		# use scrunched deflection image\n\
		-b bgfits	# background fits image to paint onto\n\
\n\
DESCRIPTION\n\
	Mapbynumericdef maps a source image to a target image using a numerically\n\
	defined deflection function.  It reads a N1 x N2 image fs = 'sourcefits' and\n\
	a M1 x M2 x 2 image d = 'deflectionfits' whose 0th ans 1st planes are\n\
	x and y components of a deflection function.\n\
	By default it generates a M1 * M2 image ft(r) = fs(r + d(r)), i.e.\n\
	the deflection is given as a function of target image coordinates.\n\
\n\
	If either of the source images is given as 'stdin' then the image\n\
	will be read from standard input.\n\
\n\
	With -i option it performs 'inverse' mapping --- i.e. deflection\n\
	is given as a function of source image --- so the resulting image\n\
	satisfies ft(x + d(x)) = fs(x) * M, where M is the magnification.\n\
	With this option, the deflection\n\
	function image size should be matched to the source image so\n\
	M1 = N1 and M2 = N2, and the output image size is L1 x L2.\n\
	Physically, inverse mapping projects a source image through a\n\
	lens.\n\
\n\
	Use -m flag to specify mapping mode:\n\
		mode = 0:	# nearest pixel\n\
		mode = 1:	# linear interpolation\n\
		mode = 2:	# sum over triangles\n\
	the default being linear interpolation.  With inverse mapping\n\
	this flag is ignorred and triangle mapping is used (which\n\
	can be very slow).\n\
\n\
	The -s option allows you to read a miniature deflection image\n\
	which has been scrunched by a factor 2^sfac. This is currently\n\
	not implemented for inverse mapping.\n\
\n\
	The -b option allows you to supply a background image to be\n\
	painted onto. The size of this image will override either the\n\
	dimensions determined from the deflection image size for\n\
	forward mapping or the L1, L2 values set with '-i' option.\n\
\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@ifa.hawaii.edu\n\
\n\n\n"		
	

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../utils/error.h"
#include "../imlib/fits.h"
#include "../imlib/map.h"

static int 	N1, N2, M1, M2, L1, L2, MM1, MM2, oflagset;
static float	**defx, **defy, scalefac;

int	deflection(float ri, float rj, float *di, float *dj);
int	invdeflection(float ri, float rj, float *di, float *dj);
void	magicinit(float **f, int N1, int N2);


main(int argc, char *argv[])	
{
	int		arg = 1; 
	fitsheader	*fitssource, *fitstarget, *fitsdef;
	float		**fsource, **ftarget, **def;
	FILE		*fsourcefile, *deffile, *bgimagef;
	int		doinversemap = 0, mapmode, initmagic;
	int		x, y, is, sfac, backgroundimage;
	int		targetN1, targetN2;
	char		*fsourcefilename, *deffilename;


	/* defaults */	
	mapmode = FAST_MAP_MODE;
	initmagic = 0;
	sfac = 0;
	backgroundimage = 0;
	oflagset = 0;

	/* parse args */
	while (arg < argc) {
		if (*argv[arg] == '-') {
			switch (*(argv[arg++]+1)) {
				case 'i':
					doinversemap = 1;
					if ((argc - arg) < 2)
						error_exit(usage);
					if (1 != sscanf(argv[arg++], "%d", &L1))
						error_exit(usage);
					if (1 != sscanf(argv[arg++], "%d", &L2))
						error_exit(usage);
					break;
                        	case 'm':
 					if ((argc - arg) < 1)
  						error_exit(usage);
                             		if (1 != sscanf(argv[arg++], "%d", &mapmode))
                                                error_exit(usage);
                                	break;
				case 'M':
 					initmagic = 1;
                               		break;
                        	case 's':
 					if ((argc - arg) < 1)
  						error_exit(usage);
                             		if (1 != sscanf(argv[arg++], "%d", &sfac))
                                                error_exit(usage);
                                	break;
                        	case 'b':
 					if ((argc - arg) < 1)
  						error_exit(usage);
					bgimagef = fopen(argv[arg++], "r");
                             		if (!bgimagef)
                                                error_exit("mapbynumericdef: failed to open background image\n");
					backgroundimage = 1;
                                	break;
				case 'o':
					oflagset = 1;
					break;
				default:
					error_exit(usage);
					break;
			}
		} else {
			break;
		}
	}
	if (doinversemap && sfac)
		error_exit("mapbynumericaldef: inverse mapping not implemented for scrunched deflection images\n");

	if ((argc - arg) != 2)
		error_exit(usage);

	/* open the source and deflection files */
	fsourcefilename = argv[arg++];
	deffilename = argv[arg++];
	if (!strcmp(fsourcefilename, deffilename)) {
		error_exit(usage);
	}
	if (strcmp(fsourcefilename, "stdin")) {
		if (!(fsourcefile = fopen(fsourcefilename, "r"))) {
			error_exit("mapbynumericdef: failed to open source image file\n");
		}
	} else {
		fsourcefile = stdin;
	}		
	if (strcmp(deffilename, "stdin")) {
		if (!(deffile = fopen(deffilename, "r"))) {
			error_exit("mapbynumericdef: failed to open deflection image file\n");
		}
	} else {
		deffile = stdin;
	}

	/* read the deflection file */
	fitsdef = readfitsheader(deffile);
	if (fitsdef->ndim != 3) {
		error_exit("mapbynumericdef: deflection file not 3-dimensional\n");
	}
	M1 = fitsdef->n[0];
	M2 = fitsdef->n[1];
	allocFloatArray(&defx, M1, M2);
	allocFloatArray(&defy, M1, M2);
	readfitsplane((void **) defx, fitsdef);
	readfitsplane((void **) defy, fitsdef);

	scalefac = 1;
	MM1 = M1;
	MM2 = M2;
	for (is = 0; is < sfac; is++) {
		scalefac *= 0.5;
		MM1 *= 2;
		MM2 *= 2;
	}

	/* read the source file */
	read2Dfloatimage(&fsource, &N1, &N2, &fitssource, fsourcefile);

	if (doinversemap) {
		if ((M1 != N1) || (M2 != N2))
			error_exit("mapbynumericdef: deflection image size does not match source image size\n");
	}

	if (doinversemap) {
		targetN1 = L1;
		targetN2 = L2;
	} else {
		targetN1 = MM1;
		targetN2 = MM2;
	}

	/* create or read the target image */
	if (backgroundimage) {
		read2Dfloatimage(&ftarget, &targetN1, &targetN2, &fitstarget, bgimagef);
	} else {
		fitstarget = copyfitsheader(fitssource);
		allocFloatArray(&ftarget, targetN1, targetN2);
		if (initmagic) {
			magicinit(ftarget, targetN1, targetN2);
		}
	}
	add_comment(argc, argv, fitstarget);

	/* apply the mapping */
	if (doinversemap) {
		set_triangle_map_mode(INVERSEMAPMODE);
		map(fsource, N1, N2, ftarget, targetN1, targetN2, invdeflection);
	} else {
		switch (mapmode) {
			case  ULTRAFAST_MAP_MODE:
				ultrafastmap(ftarget, targetN1, targetN2, fsource, N1, N2, deflection);
				break;
			case  FAST_MAP_MODE:
				fastmap(ftarget, targetN1, targetN2, fsource, N1, N2, deflection);
				break;
			case  TRIANGLE_MAP_MODE:
				map(ftarget, targetN1, targetN2, fsource, N1, N2, deflection);
				break;
			default:
				error_exit("mapbynumericdef: bad mapping mode\n");
		}
	}
	
	/* write the target image */
	write2Dfloatimage(ftarget, fitstarget);
	
	/* all done */
	exit(0);
}



int	deflection(float ri, float rj, float *di, float *dj)
{
	int	i, j;

	i = (int) floor(0.5 + scalefac * ri);
	i = (i < 0 ? 0 : i);
	i = (i >= M2 ? M2 - 1 : i);
	j = (int) floor(0.5 + scalefac * rj);
	j = (j < 0 ? 0 : j);
	j = (j >= M1 ? M1 - 1 : j);
	if ((defx[i][j] == FLOAT_MAGIC) || (defy[i][j] == FLOAT_MAGIC)) {
		return(0);
	} else {
		*di = defy[i][j];
		*dj = defx[i][j];
		return(1);
	}
}




int	invdeflection(float ri, float rj, float *di, float *dj)
{
	int	i, j;

	i = (int) floor(0.5 + ri);
	i = (i < 0 ? 0 : i);
	i = (i >= N2 ? N2 - 1 : i);
	j = (int) floor(0.5 + rj);
	j = (j < 0 ? 0 : j);
	j = (j >= N1 ? N1 - 1 : j);
	if ((defx[i][j] == FLOAT_MAGIC) || (defy[i][j] == FLOAT_MAGIC)) {
		return(0);
	} else {
		*di = defy[i][j];
		*dj = defx[i][j];
		return(1);
	}
}



void	magicinit(float **f, int N1, int N2)
{
	int	x, y;

	for (y = 0; y < N2; y++) {
		for (x = 0; x < N1; x++) {
			f[y][x] = FLOAT_MAGIC;
		}
	}
}


