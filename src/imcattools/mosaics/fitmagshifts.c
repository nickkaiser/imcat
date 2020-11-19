/*
 * fitmagshifts.c
 *
 */

#define usage "\n\
NAME\n\
	fitmagshifts --- fit for extinction/gain coefficients for a set of exposures\n\
\n\
SYNOPSIS\n\
	fitmagshifts nc ne [options...]\n\
\n\
DESCRIPTION\n\
	'fitmagshifts' reads a catalogue containing (at least) pairs\n\
	of magnitudes mag[2]; chip numbers c[2] and exposure numbers e[2]\n\
	for a set of reference stars observed on a mosaic of\n\
	nc chips with ne exposures,\n\
	and solves for any gain variations between chips and any\n\
	differential extinction between exposures.\n\
\n\
	More explicitly, we model the magnitude of a star as measured\n\
	in the c'th chip and e'th exposure as:\n\
		m_ce = m + m_c + M_e\n\
	where m is the true magnitude, and we solve for the gain\n\
	variations m_c and the extinctions M_e by least squares\n\
	minimisation (these being measured relative to the 0th\n\
	chip and 0th exposure respectively).\n\
\n\
	We output the coefficients as a pair of lc-format catalogues.\n\
	By default these are concatenated to stdout, but can be sent to\n\
	named files by using the -c and -e options.\n\
	Options are\n\
		-c chipcoefftfile\n\
		-e expcoefftfile\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n"

#include <stdio.h>
#include <math.h>
#include "../../utils/lu.h"
#include "../../utils/error.h"

main(int argc, char *argv[])
{
	int	arg = 1, m, n, *indx, mcoefft, nobjects, c, cp, e, ep;
	double	**Am, *Bm, *Cm, d;
	double	mag, magp;
	double	*dM, *dm;
	char	line[1024], *chipcoefftfilename, *expcoefftfilename, lcstring[1024];
	int	Nc, Ne, dMbase, dmbase;
	FILE	*lcpipe, *dmpipe, *dMpipe;

	/* defaults */
	chipcoefftfilename = expcoefftfilename = NULL;

	/* get ne, nc */
	if (argc < 3) {
		error_exit(usage);
	}
	if (1 != sscanf(argv[arg++], "%d", &Nc)) {
		error_exit(usage);
	}
	if (1 != sscanf(argv[arg++], "%d", &Ne)) {
		error_exit(usage);
	}
	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			fprintf(stderr, usage);
			exit(-1);
		}
		switch (argv[arg++][1]) {
			case 'c':
				chipcoefftfilename = argv[arg++];
				break;
			case 'e':
				expcoefftfilename = argv[arg++];
				break;
			default:
				fprintf(stderr, usage);
				exit(-1);
				break;
		}
	}

	mcoefft = (Ne - 1) + (Nc - 1);
	indx = (int *) calloc(mcoefft, sizeof(int));
	Bm = (double *) calloc(mcoefft, sizeof(double));
	Cm = (double *) calloc(mcoefft, sizeof(double));
	Am = (double **) calloc(mcoefft, sizeof(double));
	for (m = 0; m < mcoefft; m++) {
		Am[m] =  (double *) calloc(mcoefft, sizeof(double));
	}
	dm = (double *) calloc(Nc, sizeof(double));
	dM = (double *) calloc(Ne, sizeof(double));
	dMbase = -1;
	dmbase = Ne - 2;

	if (!(lcpipe = popen("lc -o c e mag", "r"))) {
		fprintf(stderr, "mosaicfit: unable to open lc-pipe for input\n");
		exit(-1);
	}
	nobjects = 0;
	while (fgets(line, 1024, lcpipe)) {
		if (line[0] == '#')
			continue;
		sscanf(line, "%d %d %d %d %lf %lf", &c, &cp, &e, &ep, &mag, &magp);
		if (c < 0 || c > (Nc - 1) || cp < 0 || cp > (Nc - 1)) {
			fprintf(stderr, "mosaicfit: chip number out of allowed range\n");
			exit(-1);
		}
		if (e < 0 || e > (Ne - 1) || ep < 0 || ep > (Ne - 1)) {
			fprintf(stderr, "mosaicfit: exposure number out of allowed range\n");
			exit(-1);
		}
		nobjects++;
		/* calculate C-vectors */
		/* first the exposure terms...*/
		for (m = 1; m < Ne; m++) {
			Cm[m + dMbase] = 0.0;
			if (e == m) {
				Cm[m + dMbase] += 1.0;
			}
			if (ep == m) {
				Cm[m + dMbase] -= 1.0;
			}
		}
		/* now the chip terms */
		for (n = 1; n < Nc; n++) {
			Cm[n + dmbase] = 0.0;
			if (c == n) {
				Cm[n + dmbase] += 1.0;
			}
			if (cp == n) {
				Cm[n + dmbase] -= 1.0;
			}
		}
		/* accumulate Am matrix, Bm-vector */
		for (m = 0; m < mcoefft; m++) {
			Bm[m] += (mag - magp) * Cm[m];
			for (n = 0; n < mcoefft; n++) {
				Am[m][n] += Cm[m] * Cm[n];
			}
		}
	}

	if (nobjects <= mcoefft) {
		fprintf(stderr, "mosaicfit: too few objects\n");
		exit(-1);
	}


	/* solve linear equations */
	myludcmp(Am, mcoefft, indx, &d);
	mylubksb(Am, mcoefft, indx, Bm);

	/* extract coefficients */
	for (e = 1; e < Ne; e++) {
		dM[e] = Bm[e + dMbase];
	}
	for (c = 1; c < Nc; c++) {
		dm[c] = Bm[c + dmbase];
	}

	/* and output */
	sprintf(lcstring, "lc -C -n e -n dM");
	if (expcoefftfilename) {
		strcat(lcstring, " > ");
		strcat(lcstring, expcoefftfilename);
	}
	lcpipe = popen(lcstring, "w");
	if (!lcpipe) {
		error_exit("fitmagshifts: failed to open lc-pipe for output\n");
	}
	for (e = 0; e < Ne; e++) {
		fprintf(lcpipe, "%d %13.8lg\n", e, dM[e]);
	}
	pclose(lcpipe);

	sprintf(lcstring, "lc -C -n c -n dm");
	if (chipcoefftfilename) {
		strcat(lcstring, " > ");
		strcat(lcstring, chipcoefftfilename);
	}
	lcpipe = popen(lcstring, "w");
	if (!lcpipe) {
		error_exit("fitmagshifts: failed to open lc-pipe for output\n");
	}
	for (c = 0; c < Nc; c++) {
		fprintf(lcpipe, "%d %13.8lg\n", c, dm[c]);
	}
	pclose(lcpipe);

	exit(0);
}
