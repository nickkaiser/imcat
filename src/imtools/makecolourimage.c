#define usage "\n\n\n\
NAME\n\
	makecolourimage --- combines 2 or 3 images into a color fits image\n\
\n\
SYNOPSIS\n\
	makecolourimage [options....] \n\
		-o offsetsfile		# file containing positional offsets\n\
		-f fits1 colorname1	# name of 1st image and its colour\n\
		-f fits2 colorname2	# name of 2nd image and its colour....\n\
\n\
DESCRIPTION\n\
	\"makecolourimage\" combines 2 or 3 images into a single color fits\n\
	file using offsets in offsets.out if supplied.\n\
	Names are arbitrary.  Files should be supplied in increasing\n\
	order of bandpass frequency.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


#include "../imlib/fits.h"
#include "../utils/error.h"

#define	MAX_COLOURS	3

int		main(int argc, char *argv[])	
{
	int	arg = 1, image, nimages, i, j, ii, jj, getoffsets;
	FILE	*ipf, *offsetf;
	int	N1, N2, globalioff, globaljoff;
	int	ioffset[MAX_COLOURS], joffset[MAX_COLOURS];
	fitsheader	*fitsin, *fitsout, *fits;
	fitscomment	*com;
	char	argstring[COM_LENGTH], errorstring[1024], line[1024];
	float	*fin, **fout;
        char    filename[MAX_COLOURS][COM_LENGTH], colorname[MAX_COLOURS][COM_LENGTH];
	char	*offsetsfilename;

	if (argc < 7)
		error_exit(usage);

	/* defaults */
	globalioff = globaljoff = 0;
	getoffsets = 0;

	/* process arguments */
	nimages = 0;
        while (arg < argc) {
                if (argv[arg][0] == '-') {
                        switch (argv[arg++][1]) {
                                case 'f':
                                        if (1 != sscanf(argv[arg++], "%s", filename[nimages]))
                                                error_exit(usage);
                                        if (1 != sscanf(argv[arg++], "%s", colorname[nimages]))
                                                error_exit(usage);
                                        nimages ++;
                                        break;
                                 case 'o':
                                        getoffsets = 1;
					offsetsfilename = argv[arg++];
                                        break;
                                default:
                                        error_exit(usage);
                                        break;
                        }
		}
        }

	/* read the offsets */
	for (image = 0; image < nimages; image++) {
		ioffset[image] = 0;
		joffset[image] = 0;
	}
	if (getoffsets) {
		offsetf = fopen(offsetsfilename, "r");
		if (!offsetf) 
			error_exit("makecolourimage: unable to open offsets.out for input\n");
		fgets (line, 1024, offsetf);
		for (image = 0; image < nimages; image++) {
			if (!fgets (line, 1024, offsetf))
				error_exit("makecolourimage: problem reading from offsets.out\n");
			if (2 != sscanf(line, "%d %d", &(ioffset[image]), &(joffset[image])))
				error_exit("makecolourimage: problem reading from offsets.out\n");
			ioffset[image] += globalioff;
			joffset[image] += globaljoff;
		}
	}

	/* read the first header to get image size, comments */
	if (!(ipf = fopen(filename[0], "r"))) {
		fprintf(stderr, "makecolourimage: failed to open fits file %s\n", filename[0]);
		exit(1);
	}
	fitsin = readfitsheader(ipf);
	N1 = fitsin->n[0];
	N2 = fitsin->n[1];
	fitsout = copyfitsheader(fitsin);
	fitsout->ndim = 3;
	fitsout->n[2] = nimages;
	fin = (float *) calloc(N1, sizeof(float));
	allocFloatArray(&fout, N1, N2);

	/* write the output file header */
	com = newtextcomment("HISTORY", "", "");
	switch (nimages) {
		case 2:
			sprintf(com->value, "color: %s %s", colorname[0], colorname[1]);
			break;
		case 3:
			sprintf(com->value, "color: %s %s %s", colorname[0], colorname[1], colorname[2]);
			break;
		default:
			break;
	}

	appendcomment(com, fitsout);
	writefitsheader(fitsout);

	/* now write the first file */
	image = 0;
	for (i = 0; i < N2; i++) {
		ii = i + ioffset[image];
		readfitsline(fin, fitsin);
		for (j = 0; j < N1; j++) {
			jj = j + joffset[image];
			if (ii >= 0 && ii < N2 && jj >= 0 && jj < N1)
				fout[ii][jj] = fin[j];
		}
	}
	for (i = 0; i < N2; i++)
		writefitsline(fout[i], fitsout);
	fclose(ipf);
	/* read and write all the rest of the files */
	for (image = 1; image < nimages; image++) {
		if (!(ipf = fopen(filename[image], "r"))) {
			fprintf(stderr, "makecolourimage: failed to open fits file %s\n", filename[image]);
			exit(1);
		}
		fits = readfitsheader(ipf);
		if (fits->n[0] != fitsin->n[0] || fits->n[1] != fitsin->n[1])
			error_exit("makecolourimage: input images must have same size\n");
		for (i = 0; i < N2; i++) {
			ii = i + ioffset[image];
			readfitsline(fin, fits);
			for (j = 0; j < N1; j++) {
				jj = j + joffset[image];
				if (ii >= 0 && ii < N2 && jj >= 0 && jj < N1)
					fout[ii][jj] = fin[j];
			}
		}
		for (i = 0; i < N2; i++)
			writefitsline(fout[i], fitsout);
		fclose(ipf);
	}
	writefitstail(fitsout);
	exit(0);
}
