/*
 * fitmagshifts.c
 *
 */

#define usage "\n\
NAME\n\
	fitmagshifts --- fit for extinction/gain coefficients for a set of exposures\n\
\n\
SYNOPSIS\n\
	fitmagshifts nc\n\
\n\
DESCRIPTION\n\
	'fitmagshifts' reads a catalogue containing (at least) pairs\n\
	of magnitudes mag[2], chip numbers c[2], and positions x[2][2]\n\
	for a set of reference stars observed on a mosaic of nc chips\n\
	and solves for any gain variations between chips and gradients\n\
	thereof.\n\
\n\
	More explicitly, we model the magnitude of a star as measured\n\
	at position x, y on the c'th chip as:\n\
		m_ce = m + m_c + x * m_cx + y * m_cy\n\
	where m is the true magnitude, and we solve for the coefficients.\n\
	We set M_0 and m_0 to be zero.\n\
\n\
	We output the coefficients as a lc-format catalog.\n\
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
	int	arg = 1, m, n, *indx, mcoefft, nobjects, c, cp;
	double	**Am, *Bm, *Cm, d;
	double	mag, magp, x, y, xp, yp;
	double	*dM, *dm, *dmx, *dmy;
	char	line[1024], lcstring[1024];
	int	Nc, dmbase, dmxbase, dmybase;
	FILE	*lcpipe, *dmpipe, *dMpipe;

	/* get Nc */
	if (argc != 2) {
		error_exit(usage);
	}
	if (1 != sscanf(argv[arg++], "%d", &Nc)) {
		error_exit(usage);
	}

	mcoefft = 3 * Nc - 1;
	indx = (int *) calloc(mcoefft, sizeof(int));
	Bm = (double *) calloc(mcoefft, sizeof(double));
	Cm = (double *) calloc(mcoefft, sizeof(double));
	Am = (double **) calloc(mcoefft, sizeof(double));
	for (m = 0; m < mcoefft; m++) {
		Am[m] =  (double *) calloc(mcoefft, sizeof(double));
	}
	dm = (double *) calloc(Nc, sizeof(double));
	dmx = (double *) calloc(Nc, sizeof(double));
	dmy = (double *) calloc(Nc, sizeof(double));
	dmbase = -1;
	dmxbase = dmbase + Nc;
	dmybase = dmbase + 2 * Nc;

	if (!(lcpipe = popen("lc -o c mag x", "r"))) {
		fprintf(stderr, "mosaicfit: unable to open lc-pipe for input\n");
		exit(-1);
	}
	nobjects = 0;
	while (fgets(line, 1024, lcpipe)) {
		if (line[0] == '#')
			continue;
		sscanf(line, "%d %d %lf %lf %lf %lf %lf %lf", 
			&c, &cp, &mag, &magp, &x, &y, &xp, &yp);
		if (c < 0 || c > (Nc - 1) || cp < 0 || cp > (Nc - 1)) {
			fprintf(stderr, "mosaicfit: chip number out of allowed range\n");
			exit(-1);
		}
		nobjects++;
		/* calculate C-vectors */
		/* constant chip terms */
		for (n = 1; n < Nc; n++) {
			Cm[n + dmbase] = 0.0;
			if (c == n) {
				Cm[n + dmbase] += 1.0;
			}
			if (cp == n) {
				Cm[n + dmbase] -= 1.0;
			}
		}
		/* gradient terms */
		for (n = 0; n < Nc; n++) {
			Cm[n + dmxbase] = 0.0;
			if (c == n) {
				Cm[n + dmxbase] += x;
			}
			if (cp == n) {
				Cm[n + dmxbase] -= xp;
			}
		}
		for (n = 0; n < Nc; n++) {
			Cm[n + dmybase] = 0.0;
			if (c == n) {
				Cm[n + dmybase] += y;
			}
			if (cp == n) {
				Cm[n + dmybase] -= yp;
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

	if (nobjects < mcoefft) {
		fprintf(stderr, "mosaicfit: too few objects\n");
		exit(-1);
	}


	/* solve linear equations */
	myludcmp(Am, mcoefft, indx, &d);
	mylubksb(Am, mcoefft, indx, Bm);

	/* extract coefficients */
	for (c = 1; c < Nc; c++) {
		dm[c] = Bm[c + dmbase];
	}
	for (c = 0; c < Nc; c++) {
		dmx[c] = Bm[c + dmxbase];
	}
	for (c = 0; c < Nc; c++) {
		dmy[c] = Bm[c + dmybase];
	}

	/* and output */
	sprintf(lcstring, "lc -C -n c -n dm -n dmx -n dmy");
	lcpipe = popen(lcstring, "w");
	if (!lcpipe) {
		error_exit("fitmagshifts: failed to open lc-pipe for output\n");
	}
	for (c = 0; c < Nc; c++) {
		fprintf(lcpipe, "%d %13.8lg %13.8lg %13.8lg\n", c, dm[c], dmx[c], dmy[c]);
	}
	pclose(lcpipe);

	exit(0);
}
