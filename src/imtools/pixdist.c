#define usage "\n\n\n\
NAME\n\
	pixdist - compute histogram of pixel values in a FITS image\n\
\n\
SYNOPSIS\n\
	pixdist [options...]\n\
		-d fmin fmax	# histogram range\n\
		-n nbins	# number of bins (40)\n\
\n\
DESCRIPTION\n\
	\"pixdist\" generates a histogram of pixel values from a fits image.\n\
	Use -d option to specify range of pixvalues (default = fmin fmax)\n\
	Use -n option to specify number of bins.\n\
\n\
	Output is a lc-format catalogue with object items f, n\n\
	giving the bin center and count\n\
	and with header items overcount, undercount \n\
	giving the count of pixels outside the histogram range.\n\
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
#include "../utils/arrays.h"

int		main(int argc, char *argv[])	
{
	int		N1, N2, *count, val, arg = 1, autorange;
	fitsheader	*fits;
	float		**f;
	float		fmin, fmax, df;
	int		bin, nbins, overcount, undercount;
	int		i, j;
	FILE		*opf;
	char		lccom[1024], argstring[128];

	/* defaults */
	autorange = 1;
	nbins = 40;

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-')
			error_exit(usage);
		switch(argv[arg++][1]) {
			case 'd':
				autorange = 0;
				if (argc < arg + 2)
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%f", &fmin) ||
				1 != sscanf(argv[arg++], "%f", &fmax))
						error_exit(usage);
				break;
			case 'n':
				sscanf(argv[arg++], "%d", &nbins);
				break;
			default:
				error_exit(usage);
				break;
		}
	}

	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);

	/* find min and max */
	if (autorange) {
		fmin =  1.e20;
		fmax = -1.e20;
		for (i = 0; i < N2; i++) {
			for (j = 0; j < N1; j++) {
				if (f[i][j] == FLOAT_MAGIC)
					continue;
				fmax = (f[i][j] > fmax ? f[i][j] : fmax);
				fmin = (f[i][j] < fmin ? f[i][j] : fmin);
			}
		}
	}

	df = (fmax - fmin) / (double) nbins;
	count = (int *) calloc(nbins, sizeof(int));
	overcount = undercount = 0;
	for (i = 0; i < N2; i++) {
		for (j = 0; j < N1; j++) {
			bin = floor((f[i][j] - fmin)/ (double) df);
			if (bin < 0) {
				undercount++;
				continue;
			}
			if (bin >= nbins) {
				overcount++;
				continue;
			}
			count[bin]++;
		}
	}
	
	argsToString(argc, argv, argstring);
	sprintf(lccom, "lc -C -n f -n n -x -a '%s' -H 'overcount = %d' -H 'undercount = %d'", 
		argstring, overcount, undercount);
	opf = popen(lccom, "w");
	if (!opf) {
		error_exit("pixdist: failed to open lc-pipe for output\n");
	}
/*
	fprintf(stdout, "# pixdist output\n");
	fprintf(stdout, "# %10d pixels below %g\n", undercount, fmin);
	fprintf(stdout, "# %10d pixels above %g\n", overcount, fmax);
	fprintf(stdout, "#   pixval      count\n");
*/
	for (bin = 0; bin < nbins; bin++) {
		fprintf(opf, "%13.8g %10d\n", fmin + (bin + 0.5) * df, count[bin]);
	}
	pclose(opf);
	exit(0);
}



