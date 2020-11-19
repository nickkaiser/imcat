/*
 *	mova.c
 *
 *	calculate moving average of a table
 */

#define usage "\n\
NAME\n\
	mova - calculate moving average of tabulated data\n\
\n\
SYNOPSIS\n\
	mova xmin xmax nbins binwidth\n\n\
\n\
DESCRIPTION\n\
	mova calculates average of y-values (column 2) binned according\n\
	to x-value (column 1) with nbins bins of width binwidth between xmin and xmax.\n\
	It outputs a table with <y>,x,n for non empty bins.\n\
\n\
AUTHOR\n\
	Nick Kaiser - kaiser@hawaii.edu\n\
\n"

#include <stdio.h>
#include <math.h>


main (int argc, char *argv[])
{
	int	thebin, bin, nbins, dbin, *sumn;
	double	x, y, xmin, xmax, dx, width, halfwidth, *sumy;
	char	line[1024];
 
	if (argc != 5) {
		fprintf(stderr, usage);
		exit(-1);
	}

	sscanf(argv[1], "%lf", &xmin);
	sscanf(argv[2], "%lf", &xmax);
	sscanf(argv[3], "%d", &nbins);
	sscanf(argv[4], "%lf", &width);
	dx = (xmax - xmin) / nbins;
	dbin = ceil(0.5 * width / dx);

	sumy = (double *) calloc(nbins, sizeof(double));
	sumn = (int *) calloc(nbins, sizeof(int));

	while (fgets(line, 1024, stdin)) {
		if (line[0] == '#')
			continue;
		sscanf(line, "%lf %lf", &x, &y);
		thebin = floor(0.5 + (x - xmin) / dx);
		for (bin = thebin - dbin; bin <= thebin + dbin; bin++) {
			if (bin < 0 || bin >= nbins)
				continue;
			sumn[bin]++;
			sumy[bin] += y;
		}
	}
	fprintf(stdout, "# mova output\n#        x        <y>          N\n");
	for (bin = 0; bin < nbins; bin++) {
		if (sumn[bin])
			fprintf(stdout, "%10.5le %10.5le %5d\n",
				xmin + bin * dx, sumy[bin] / sumn[bin], sumn[bin]);
	}
	exit(0);
}
