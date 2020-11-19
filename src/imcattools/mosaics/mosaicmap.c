/*
 * mosaicmap.c
 *
 * maps a source fits file to a target fits file using transformation defined by
 * parameters alpha Phi dX dY phi dx dy (see mosaicfitting.tex)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../utils/arrays.h"
#include "../../imlib/fits.h"
#include "../../imlib/map.h"

#define usage  "\n\n\
NAME\n\
	mosaicmap --- maps fits file - called by makemosaicstack\n\
\n\
SYNOPSIS\n\
	mosaicmap alpha Phi00 Phi01 Phi10 Phi11 dX dY phi dx dy DX DY X0 Y0 scalefac NX NY target.fits [mapmode]\n\
\n\
DESCRIPTION\n\
	'mosaicmap' maps a source fits file to a target file using transformation defined by\n\
	parameters alpha Phi_ij dX dY phi dx dy (see mosaicfitting.tex for more details)\n\
	DX, DY are the nominal coords of the bottom left corner of the chip,\n\
	X0, Y0 are the coords of the bottom left corner of the sub-image we are making, and\n\
	NX, NY are nominal dimensions of chips (as defined in nominal.db).\n\
	All above dimensions are given in source pixel size units.\n\
	'scalefac' is the target pixel size in units of source pixel size.\n\
	'mapmode' is 0,1,2 for nearest pixel, linear interpolation and triangular\n\
	tesselation mapping modes respctively. Default is linear interpolation.\n\
	Mosiacmap is usually invoked by the script 'makemosaicstack' or somesuch.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n"

/* global parameters - needed by deflection() */
double	alpha, Phi00, Phi01, Phi10, Phi11, dX, dY, phi, dx, dy, DX, DY, X0, Y0, scalefac;
int	NXnom, NYnom;

void	deflect(float ri, float rj, float *di, float *dj);

main(int argc, char *argv[]) {
	int	i, j, NX, NY, NXsrc, NYsrc, comctarg, comcsrc, mapmode; 
	char *comvtarg[MAX_COMMENTS], *comvsrc[MAX_COMMENTS], targetfilename[1024];
	float   **fsource, **ftarget;
	FILE 	*ipf, *opf;

	/* defaults */
	mapmode = FAST_MAP_MODE;

	if (argc < 19 || argc > 20) {
		fprintf(stderr, usage);
		exit(-1);
	}
	sscanf(argv[1], "%lf", &alpha);
	sscanf(argv[2], "%lf", &Phi00);
	sscanf(argv[3], "%lf", &Phi01);
	sscanf(argv[4], "%lf", &Phi10);
	sscanf(argv[5], "%lf", &Phi11);
	sscanf(argv[6], "%lf", &dX);
	sscanf(argv[7], "%lf", &dY);
	sscanf(argv[8], "%lf", &phi);
	sscanf(argv[9], "%lf", &dx);
	sscanf(argv[10], "%lf", &dy);
	sscanf(argv[11], "%lf", &DX);
	sscanf(argv[12], "%lf", &DY);
	sscanf(argv[13], "%lf", &X0);
	sscanf(argv[14], "%lf", &Y0);
	sscanf(argv[15], "%lf", &scalefac);
	sscanf(argv[16], "%d", &NXnom);
	sscanf(argv[17], "%d", &NYnom);
	sscanf(argv[18], "%s", targetfilename);
	if (argc == 20)
		sscanf(argv[19], "%d", &mapmode);

	ipf = fopen(targetfilename, "r");
	if (!ipf) {
		fprintf(stderr, "mosaicmap: can't open target fits file for reading\n");
		exit(-1);
	}
	set_fits_ipf(ipf);

	fread_fits(&ftarget, &NX, &NY, &comctarg, comvtarg);
	fclose(ipf);
	opf = fopen(targetfilename, "w");
	if (!opf) {
		fprintf(stderr, "mosaicmap: can't open target fits file for writing\n");
		exit(-1);
	}
	set_fits_opf(opf);

	set_fits_ipf(stdin);
	fread_fits(&fsource, &NXsrc, &NYsrc, &comcsrc, comvsrc);

	/* apply the mapping */
	switch (mapmode) {
		case ULTRAFAST_MAP_MODE:
			ultrafastmap(ftarget, NX, NY, fsource, NXsrc, NYsrc, deflect);
			break;
		case FAST_MAP_MODE:
			fastmap(ftarget, NX, NY, fsource, NXsrc, NYsrc, deflect);
			break;
		case TRIANGLE_MAP_MODE:
			map(ftarget, NX, NY, fsource, NXsrc, NYsrc, deflect);
			break;
		default:
			fprintf(stderr, "mosaicmap: bad mapping mode\n");
			exit(-1);
			break;
	}

	/* write the target image */
        add_comment(argc, argv, &comctarg, comvtarg);
        fwrite_fits(ftarget, NX, NY, comctarg, comvtarg);

        /* all done */
        exit(0);
	
}


void	deflect(float ri, float rj, float *di, float *dj)
{
	double	X, Y, Xe, Ye, xe, ye, xce, yce, i, j, rr;

	X = scalefac * rj + X0;
	Y = scalefac * ri + Y0;
	X -= dX;
	Y -= dY;
	Xe = X - Phi00 * X - Phi01 * Y;
	Ye = Y - Phi10 * X - Phi11 * Y;
	rr = Xe * Xe + Ye * Ye;
	xe = (1 - alpha * rr) * Xe;
	ye = (1 - alpha * rr) * Ye;
	xe -= dx;
	ye -= dy;
	xce = xe - phi * ye;
	yce = ye + phi * xe;
	j = xce - DX;
	i = yce - DY;
/* no longer implemented
	if (orient < 1) {
		j = NXnom - j;
		i = NYnom - i;
	}
*/
	*di = i - ri;
	*dj = j - rj;
}
