#define usage "\n\
NAME\n\
	applytrans --- apply transformations from fitstack\n\
\n\
SYNOPSIS\n\
	applytrans parfile\n\
\n\
DESCRIPTION\n\
	'applytrans' reads from stdin a catalogue containing the result of\n\
	merging all pairs of cats for a stack of 'nexp' images (as created\n\
	by 'mergestacks1') and applies the transformations \n\
	defined in 'parfile' (as created by fitstack)\n\
	to generate sky coords r, so we can then filter to remove\n\
	discrepant pairs.\n\
	See also fitstack.tex.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n"

#include <stdio.h>
#include <math.h>
#include "fmode.h"

#define SCALE 1.0

main(int argc, char *argv[])
{
	int	i, j, m, n, e, ep;
	double	x[2], xp[2], r[2], rp[2], *d[2], *phi[2][2], *a[2];
	char	line[1024], *parfilename;
	int	nexp, nM, M, l, lmax, *ll, *mm;
	FILE	*lcinpipe, *lcoutpipe, *paramf;
	double	mag, magp, Mag, Magp, *dm;
	double	inputbuff[8], outputbuff[12];

	/* parse args */
	if (argc != 2 || argv[1][0] == '-')
		error_exit(usage);
	parfilename = argv[1];

	if (!(paramf = fopen(parfilename, "r"))) {
		error_exit("applytrans: unable to open file for transformation parameters\n");
	}
	fgets(line, 1024, paramf);	/* comment line */
	fgets(line, 1024, paramf);	/* nexp line */
	sscanf(line, "%d", &nexp);
	fgets(line, 1024, paramf);	/* lmax line */
	sscanf(line, "%d", &lmax);
	fgets(line, 1024, paramf);	/* comment line */
	/* allocate space for the model parameters */
	nM = lmax * (lmax + 1);		/* generous guess for number of modes */
	for (i = 0; i < 2; i++) {
		d[i] = (double *) calloc(nexp, sizeof(double));
		a[i] = (double *) calloc(nM, sizeof(double));
		for (j = 0; j < 2; j++) {
			phi[i][j] = (double *) calloc(nexp, sizeof(double));
		}
	}
	dm = (double *) calloc(nexp, sizeof(double));
	/* set up arrays of l, m values */
	ll = (int *) calloc(nM, sizeof(int));
	mm = (int *) calloc(nM, sizeof(int));
	for (e = 0; e < nexp; e++) {
		fgets(line, 1024, paramf);
		sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", d[0]+e, d[1]+e, 
			phi[0][0]+e, phi[0][1]+e, phi[1][0]+e, phi[1][1]+e, dm+e);
	}
	fgets(line, 1024, paramf);	/* comment line */
	M = 0;
	while (fgets(line, 1024, paramf)) {
		sscanf(line, "%d %d %lf %lf", ll+M, mm+M, a[0]+M, a[1]+M);
		M++;
	}
	nM = M;
	fclose(paramf);

	if (!(lcinpipe = popen("lc -o x exp mag", "r"))) {
		fprintf(stderr, "mosaicfit: unable to open lc-pipe for input\n");
		exit(-1);
	}
	if (!(lcoutpipe = popen("lc -C -N '2 2 2 x' -N '1 2 exp' -N '1 2 mag' -N '2 2 2 r' -N '1 2 Mag'", "w"))) {
		fprintf(stderr, "mosaicfit: unable to open lc-pipe for output\n");
		exit(-1);
	}
	while (fgets(line, 1024, lcinpipe)) {
		sscanf(line, "%lf %lf %lf %lf %d %d %lf %lf", x, x+1, xp, xp+1, &e, &ep, &mag, &magp);
		for (i = 0; i < 2; i++) {
			r[i] = x[i] + d[i][e];
			rp[i] = xp[i] + d[i][ep];
			for (j = 0; j < 2; j++) {
				r[i] += phi[i][j][e] * x[j];
				rp[i] += phi[i][j][ep] * xp[j];
			}
			for (M = 0; M < nM; M++) {
				r[i] += a[i][M] * f(ll[M], mm[M], x);
				rp[i] += a[i][M] * f(ll[M], mm[M], xp);
			}
		}
		fprintf(lcoutpipe, "%13.8lg %13.8lg %13.8lg %13.8lg %13d %13d %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg\n",
			x[0], x[1], xp[0], xp[1], e, ep, mag, magp, r[0], r[1], rp[0], rp[1], mag, magp);
	}
	exit(pclose(lcoutpipe));
}
