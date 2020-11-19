/*
 * iis.c
 *
 * create iis header followed by an image line
 *
 */

#define usage "\n\n\
NAME\n\
	iis - pipe fits image data into saopipe\n\
\n\
SYNOPSIS\n\
	iis [options...]\n\
	where options are:\n\
		-u		# print this message\n\
		-f fmin fmax	# limits for pixel values\n\
		-p pipename	# supply FIFO pipe name\n\
\n\
DESCRIPTION\n\
	\"iis\" reads a fits image from stdin, prepends an iis header\n\
	and writes the output to a FIFO pipe /dev/imt1o so that saoimage can\n\
	display it.\n\
	If the standard FIFO does not exist (it typically won't unless\n\
	you have IRAF installed) then create one with e.g.\n\
		'mknod ~/dev/saopipe p'\n\
	and then fire up saoimage with\n\
		'saoimage -idev ~/dev/saopipe &'\n\
	Beware: some old versions of saoimage crash when you try to\n\
	supply 'idev' argument.\n\
	Many thanks to Karl Glazebrook for refinements\n\
	to original code.\n\
\n\
AUTHOR\n\
	Nick Kaiser (kaiser@cita.utoronto.ca)\n\
\n\n"

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "libiis.h"
#include "../../imlib/fits.h"

#define MAGIC FLOAT_MAGIC
#define MINFLOAT FLT_MIN
#define MAXFLOAT FLT_MAX

main(int argc, char *argv[]) {
	float	**f, fmin, fmax;
	int	x, y, N1, N2, arg = 1, usenamedpipe;
	char	*pipename;
   	int 	fb, fbx, fby, frame, autoscale;
	fitsheader	*fits;

	fb = 1; fbx = 512; fby = 512; frame = 1;

	/* defaults */
	autoscale = 1;
	fmin = 0.0;
	fmax = 0.0;
	usenamedpipe = 0;

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			fprintf(stderr, usage);
			exit(-1);
		}
		switch (argv[arg++][1]) {
			case 'f':
				sscanf(argv[arg++], "%f", &fmin);
				sscanf(argv[arg++], "%f", &fmax);
				autoscale = 0;
				break;
			case 'p':
				pipename = argv[arg++];
				usenamedpipe = 1;
				break;
			case 'u':
			default:
				fprintf(stderr, usage);
				exit(-1);
				break;
		}
	}

	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);
	if (autoscale) {
		fmin = MAXFLOAT;
		fmax = MINFLOAT;
		for (y = 0; y < N2; y++) {
			for (x = 0; x < N1; x++) {
				if (f[y][x] != (float) MAGIC) {
					if (f[y][x] > fmax)
						fmax = f[y][x];
					if (f[y][x] < fmin)
						fmin = f[y][x];
				}
			}
		}
	}
	while (fbx < N1 || fby < N2) {
		fbx *= 2;
		fby *= 2;
		fb += 2;
	}
	if (usenamedpipe) {
		iis_open("", pipename, fb, fbx, fby);
	} else {
		iis_open("","",fb, fbx, fby);
	}
	iis_display(f, N1, N2, fmin, fmax, frame); 
	iis_close();


/*	old version:
	openiispipe();
	iisdisplay(f, N1, N2, fmin, fmax);
*/
	exit(0);
}
