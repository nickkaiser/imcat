#define usage "\n\n\n\
NAME\n\
	fitstomovie --- convert a 3 dimensional FITS file to a movie file\n\
SYNOPSIS\n\
	fitstomovie [-p paletteframe] [-c cmap] [-z zoomfac] [-g | -G dstgif | -M dstmpg] [-x]\n\
\n\
DESCRIPTION\n\
	'fitstomovie' reads a 3 dimensional FITS image from stdin and\n\
	generates a set of individual frame images in tmp/nnnnn.sfx\n\
	using fitstopnm with nnnnn = 00000, 00001, 00002....\n\
\n\
	By default sfx = ppm or pgm depending on color/grayscale image type.\n\
	With -g option we output a set of sfx = gif files\n\
	and with -G dstgif we combine them with multigif to generate\n\
	the multiframe output file dstgif and then deleting the temporary frames.\n\
	Similarly, with -M option we generate a set of .ppm files and then combine\n\
	them with 'mpeg_encode tmp/mpeg.param' where tmp/mpeg.param is generated\n\
	by fitstomovie and then deleted. The -M and -G options require that you\n\
	have installed mpeg_encode and/or multigif as apprioriate.\n\
\n\
	With -x option we don't clear up any temporary files we create.\n\
\n\
	Use -p with -G option to tell multigif to use the paletteframe'th frame for\n\
	the global palette.\n\
\n\
	Images should be scaled to min = 0, max = 255.\n\
\n\
	fitstomovie will first 'rm -f tmp/?????.sfx' so use with care.\n\
\n\
	Options are:\n\
		-c cmap		# pipe fits output through colorize cmap -f 0 255\n\
		-z zoomfac	# paint pixels of size 2^zoomfac (0)\n\
		-u		# print this message\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "utils/error.h"
#include "utils/args.h"
#include "utils/arrays.h"
#include "imlib/fits.h"
#include "mpegtemplate.h"

#define PPM_FRAME_TYPE 0
#define GIF_FRAME_TYPE 1

int		main(int argc, char *argv[])	
{
	int		p, np, nbytes, cmap, z, zoomfac, frametype, cleanup, paletteframe = 1;
	char		*dat, *sfx;
	fitsheader	*fits;
	char		*flag, sysstring[1024], tmpstring[1024], comstring[1024], comstring2[1024], *tmpdirname;
	char		*dstmpgname, *dstgifname;
	FILE		*tmpf;

	/* defaults */
	cmap = -1;
	zoomfac = 0;
	frametype = PPM_FRAME_TYPE;
	sfx = "ppm";
	dstgifname = NULL;
	tmpdirname = "./tmp";
	cleanup = 1;
	dstmpgname = NULL;
 
	/* parse args */
	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'c':
				cmap = getargi();
				break;
			case 'z':
				zoomfac = getargi();
				break;
			case 'M':
				dstmpgname = getargs();
				break;
			case 'G':
				dstgifname = getargs();
			case 'g':
				sfx = "gif";
				frametype = GIF_FRAME_TYPE;
				break;
			case 'x':
				cleanup = 0;
				break;
			case 'p':
				paletteframe = getargi();
				break;
			default:
				error_exit(usage);
		}
	}
	if ((dstgifname && dstmpgname) || (dstmpgname && (frametype == GIF_FRAME_TYPE))) {
		error_exit("fitstomovie: illegal combination of command line parameters\n");
	} 

	/* read the fits header */
	fits = readfitsheader(stdin);
	if (fits->ndim != 3) {
		error_exit("fitstomovie: image must be 3-dimensional\n");
	}
	np = fits->n[2];
	fits->ndim = 2;
	nbytes = fits->n[0] * fits->n[1] * pixsize(fits->extpixtype);

	/* allocate data */
	dat = (char *) checkalloc(nbytes, sizeof(char));

	/* remove any pre-existing frames */
	sprintf(sysstring, "rm -f %s/?????.%s", tmpdirname, sfx);
	fprintf(stderr, "fitstomovie: running %s\n", sysstring);
	if (system(sysstring)) {
		strcat(sysstring, " failed");
		error_exit(sysstring);
	}

	/* generate the output command */ 
	sprintf(comstring, "");
	for (z = 0; z < zoomfac; z++) {
		sprintf(tmpstring, "unscrunch | ");
		strcat(comstring, tmpstring);
	}
	if (cmap >= 0) {
		sprintf(tmpstring, "colorize %d -f 0 255 | ", cmap);
		strcat(comstring, tmpstring);
	} 

	sprintf(comstring2, "fitstopnm %s/tmp.fits", tmpdirname);
	if (cmap < 0) {
		/* black-white option didn't work under solaris 
		strcat(comstring, "| pgmtoppm black-white"); */
		strcat(comstring2, "| pgmtoppm rgbi:0.0/0.0/0.0-rgbi:1.0/1.0/1.0 ");
	}
	switch (frametype) {
		case PPM_FRAME_TYPE:
			break;
		case GIF_FRAME_TYPE:
			strcat(comstring2, "| ppmquant 256 | ppmtogif ");
			break;
		default:
			error_exit("fitstomovie: illegal frame type\n");
	}
	fprintf(stderr, "fitstomovie: command: %s\n", comstring);
	fprintf(stderr, "fitstomovie: command2: %s\n", comstring2);

	/* process the frames */
	for (p = 0; p < np; p++) {
		sprintf(sysstring, "%s cat > %s/tmp.fits", comstring, tmpdirname, p, sfx);
		fprintf(stderr, "fitstomovie: opening %s\n", sysstring);
		tmpf = popen(sysstring, "w");
		if (!tmpf) {
			strcat(sysstring, " failed");
			error_exit(sysstring);
		}
		fits->opstream = tmpf;
		writefitsheader(fits);
		fread(dat, sizeof(char), nbytes, stdin);
		fwrite(dat, sizeof(char), nbytes, tmpf);
		pclose(tmpf);
		sprintf(sysstring, "%s > %s/%.5d.%s", comstring2, tmpdirname, p, sfx);
		if (system(sysstring)) {
			strcat(sysstring, " failed");
			error_exit(sysstring);
		}
	}

	if (dstgifname) {
		/* make multigif */
		sprintf(sysstring, "multigif -l 1-%d -n 1000 -o %s -g %d %s/?????.gif", np, dstgifname, paletteframe, tmpdirname);
		fprintf(stderr, "fitstomovie: running %s\n", sysstring);
		if (system(sysstring)) {
			strcat(sysstring, " failed");
			error_exit(sysstring);
		}
	}
	if (dstmpgname) {
		/* make mpg */
		sprintf(tmpstring, "sed 's:tmp.mpg:%s:g' | sed 's:99999:%.5d:g' > %s/mpeg.param", dstmpgname, np - 1, tmpdirname);
		tmpf = popen(tmpstring, "w");
		if (!tmpf) {
			fprintf(stderr, "fitstomovie: failed to open pipe '%s'\n", tmpstring);
		}
		fprintf(tmpf, mpegparamtext);
		pclose(tmpf);
		sprintf(sysstring, "mpeg_encode %s/mpeg.param", tmpdirname);
		fprintf(stderr, "fitstomovie: running %s\n", sysstring);
		if (system(sysstring)) {
			strcat(sysstring, " failed");
			error_exit(sysstring);
		}
	}
	if (cleanup && (dstgifname || dstmpgname)) {
		sprintf(sysstring, "rm -f %s/?????.%s", tmpdirname, sfx);
		fprintf(stderr, "fitstomovie: running %s\n", sysstring);
		if (system(sysstring)) {
			strcat(sysstring, " failed");
			error_exit(sysstring);
		}
	}
	fprintf(stderr, "fitstomovie: all done\n");
	exit(0);
}

