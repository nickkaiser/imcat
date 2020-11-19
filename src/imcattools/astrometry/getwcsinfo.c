/*
 * getwcsinfo.c
 */

#define usage "\n\
NAME\n\
	getwcsinfo --- report on WCS fits header info\n\
\n\
SYNOPSIS\n\
	getwcsinfo\n\
\n\
DESCRIPTION\n\
	getwcsinfo reads a fits header from stdin, runs Doug\n\
	Mink's wcsinit() and reports what, if any, world\n\
	coordinate sysntem information is present.\n\
\n\
AUTHOR\n\
	Nick Kaiser	kaiser@hawaii.edu\n\
\n\n"


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


#include "utils/error.h"
#include "wcs.h"

main (int argc, char *argv[])
{
	struct WorldCoor 	*wcs;
	char 			header[100*2880];
	int			i;


	if (argc != 1) {
		error_exit(usage);
	}

	for(i = 0; i < 3600; i++) {
		fgets(header + 80 * i, 81, stdin);
		if(!strncmp(header + 80 * i, "END", 3)) break;
	}

	wcs = wcsinit(header);
	if (nowcs(wcs)) {
		fprintf(stdout, "# getwcsinfo: header does not contain wcs info\n");
		exit(0);
	}
	fprintf(stdout, "# getwcsinfo:\n");
	fprintf(stdout, "# projection type:    %s\n", wcs->ptype);
	fprintf(stdout, "# x-coordinate type:  %s\n", wcs->c1type);
	fprintf(stdout, "# y-coordinate type:  %s\n", wcs->c2type);
	fprintf(stdout, "# reference frame:    %s\n", wcs->radecsys);
	fprintf(stdout, "# output ref frame:   %s\n", wcs->sysout);
	fprintf(stdout, "# center:             %s\n", wcs->center);
	fprintf(stdout, "# xref:               %13.8lg\n", wcs->xref);
	fprintf(stdout, "# yref:               %13.8lg\n", wcs->yref);
	fprintf(stdout, "# xrefpix:            %13.8lg\n", wcs->xrefpix);
	fprintf(stdout, "# yrefpix:            %13.8lg\n", wcs->yrefpix);
	if (wcs->rotmat) {
		fprintf(stdout, "# cd matrix           {{%13.8lg, %13.8lg},{%13.8lg, %13.8lg}}\n",
			wcs->cd11, wcs->cd12, wcs->cd21, wcs->cd22);
	} else {
		fprintf(stdout, "# rotation:           %13.8lg\n", wcs->rot);
	}
	fprintf(stdout, "# equinox:            %13.8lg\n", wcs->equinox);
	fprintf(stdout, "# epoch:              %13.8lg\n", wcs->epoch);
	fprintf(stdout, "# plate_ra:           %13.8lg\n", wcs->plate_ra);
	fprintf(stdout, "# plate_dec:          %13.8lg\n", wcs->plate_dec);
	fprintf(stdout, "# plate_scale:        %13.8lg\n", wcs->plate_scale);


	exit(0);
}
