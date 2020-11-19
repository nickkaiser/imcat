/*
 * beta.c
 *
 * calculates beta parameter
 */

#include <stdio.h>
#include <math.h>

#define usage "\n\
NAME\n\
	rebin - rebin x-y table to logarithmically spaced intervals in x\n\
\n\
SYNOPSIS\n\
	rebin xname yname xmin xmax dlnx\n\
\n\
DESCRIPTION\n\
	Rebin reads a from stdin a lc-format catalogue containing\n\
	x-y pairs, with x values in increasing order. It then\n\
	outputs to stdout a table of x and linearly interpolated y values where\n\
	the x values are now uniformly distributed between xmin and xmax\n\
	with log spacing dlnx.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n"

main(int argc, char *argv[]) 
{
	FILE	*inpipe, *outpipe;
	char	*xname, *yname, lccom[128];
	double	xmin, xmax, dlnx, xlo, ylo, xhi, yhi, x, y;
	int	i;
	
	/* parse args */
	if (argc != 6) {
		fprintf(stderr, usage);
		exit(-1);
	}
	xname = argv[1];
	yname = argv[2];
	sscanf(argv[3], "%lf", &xmin);
	sscanf(argv[4], "%lf", &xmax);
	sscanf(argv[5], "%lf", &dlnx);

	/* open input and output streams */
	sprintf(lccom, "lc -o %s %s", xname, yname);
	inpipe = popen(lccom, "r");
	if (!inpipe) {
		fprintf(stderr, "rebin: failed to open input pipe\n");
		exit(-1);
	}
	sprintf(lccom, "lc -C -n %s -n %s", xname, yname);
	outpipe = popen(lccom, "w");
	if (!outpipe) {
		fprintf(stderr, "rebin: failed to open output pipe\n");
		exit(-1);
	}
	
	/* read the first input line */
	fscanf(inpipe, "%lf %lf", &xlo, &ylo);

	/* loop over output values till we exceed first xlo outputing y = 0 */
	i = 1;
	for (x = xmin; x < xlo; x = xmin * exp(i++ * dlnx)) {
		fprintf(outpipe, "%13.8lg %13.8lg\n", x, 0.0);
	}

	/* read the rest of the catalogue */
	while (2 == fscanf(inpipe, "%lf %lf", &xhi, &yhi)) {
		for (; x <= xhi; x = xmin * exp(i++ * dlnx)) {
			if (x > xmax) {
				pclose(inpipe);
				pclose(outpipe);
				exit(0);
			}					
			if (xhi > xlo) {
				y = ylo + (yhi - ylo) * (x - xlo) / (xhi - xlo);
			} else {
				y = ylo;
			}
			fprintf(outpipe, "%13.8lg %13.8lg\n", x, y);
		}
		xlo = xhi;
		ylo = yhi;
	}

	/* loop over rest of output values outputing y = 0 */
	for (; x <= xmax; x= xmin * exp(i++ * dlnx)) {
		fprintf(outpipe, "%13.8lg %13.8lg\n", x, 0.0);
	}

	pclose(inpipe);
	pclose(outpipe);
	exit(0);	
}
