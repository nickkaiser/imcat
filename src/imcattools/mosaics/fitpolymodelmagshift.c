/*
 * fitpolymodelmagshift.c
 *
 */

#define usage "\n\
NAME\n\
	fitpolymodelmagshift --- fit for extinction/gain coefficients for a set of exposures\n\
\n\
SYNOPSIS\n\
	fitpolymodelmagshift ne fitorder opname [-outputarray] [-noextinction]\n\
\n\
DESCRIPTION\n\
	'fitpolymodelmagshift' reads a catalogue containing (at least) pairs\n\
	of magnitudes mag[2]; detector coords xdet[2][2] and exposure numbers e[2]\n\
	for a set of reference stars observed in ne exposures.\n\
\n\
	It solves for differential extinction between exposures and for a\n\
	spatial polynomial magnitude shift.\n\
\n\
	More explicitly, we model the magnitude of the i'th star as measured\n\
	at position xdet and e'th exposure as:\n\
		m_ei = m + m_e + sum_j m_j f_j(xdet)\n\
	where m is the true magnitude and the functions f_j are\n\
	polynomials (no DC term).\n\
\n\
	We solve for the coefficients by least squares.\n\
	minimisation (these being measured relative to the 0th\n\
	chip and 0th exposure respectively).\n\
\n\
	We output the m_e coefficients as a lc-format catalogue opname.cat.\n\
	The coefficents m_j are output as a l-model par file opname.par.\n\
\n\
	With -outputarray option we output the A-matrix etc to stdout.\n\
\n\
	With -noextinction option we don't solve for the extinction terms.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n"

#include <stdio.h>
#include <math.h>
#include "../../utils/lu.h"
#include "../../utils/error.h"
#include "../../utils/args.h"
#include "../../catlib/cat.h"
#include "../../utils/lmodel.h"

double f(int m, int *e, double **xdet);

int	ne, fitorder, fitextinction;
lmodel	*thelmodel;

main(int argc, char *argv[])
{
	int	arg = 1, m, mm, n, *indx, mcoefft, e[2], ee, nitems, nobjects;
	double	**A, *B, det, *dme, *s, *modefunc;
	double	**xdet, *mag, *ipbuff;
	char	*flag, *opname, opcatname[1024], opparname[1024], lcstring[1024];
	FILE	*lcpipe, *oppar;
	item	*xdetitem, *aitem;
	int	outputarray;

	/* defaults */
	outputarray = 0;
	fitextinction = 1;

	/* get args */
	argsinit(argc, argv, usage);
	if (FLAG_ARG == nextargtype()) {
		error_exit(usage);
	}
	ne = getargi();
	fitorder = getargi();
	opname = getargs();
	while (flag = getflag()) {
		if (!strcmp(flag, "outputarray")) {
			outputarray = 1;
		} else if (!strcmp(flag, "noextinction")) {
			fitextinction = 0;
		} else {
			error_exit(usage);
		}
	}
	if (fitextinction) {
		if (fitorder < 2) {
			error_exit("fitpolymodelmagshift: fitorder must exceed 2\n");
		}
	} else {
		if (fitorder < 1) {
			error_exit("fitpolymodelmagshift: fitorder must exceed 1\n");
		}
	}

	/* create op file names */
	sprintf(opcatname, "%s.cat", opname);
	sprintf(opparname, "%s.par", opname);		

	/* create the lmodel structure */
	xdetitem = newitem("xdet", NUM_TYPE, 1, 2);
	aitem = newitem("dm", NUM_TYPE, 1, 1);
	thelmodel = newpolylmodel(xdetitem, aitem, 1 + fitextinction, fitorder);
	fprintf(stderr, "# created model with %d modes\n", thelmodel->nmodes);

	/* create arrays, matrices */
	nitems = 8;
	ipbuff = (double *) calloc(nitems, sizeof(double));
	mcoefft = thelmodel->nmodes + fitextinction * (ne - 1);
	indx = (int *) calloc(mcoefft, sizeof(int));
	B = (double *) calloc(mcoefft, sizeof(double));
	s = (double *) calloc(mcoefft, sizeof(double));
	modefunc = (double *) calloc(mcoefft, sizeof(double));
	A = (double **) calloc(mcoefft, sizeof(double *));
	for (m = 0; m < mcoefft; m++) {
		A[m] =  (double *) calloc(mcoefft, sizeof(double));
	}
	dme = (double *) calloc(ne, sizeof(double));
	xdet = (double **) calloc(2, sizeof(double *));

	if (!(lcpipe = popen("lc -b -o e xdet mag", "r"))) {
		fprintf(stderr, "fitpolymodelmagshift: unable to open lc-pipe for input\n");
		exit(-1);
	}
	nobjects = 0;
	while (nitems == fread(ipbuff, sizeof(double), nitems, lcpipe)) {
		e[0] = (int) floor(ipbuff[0]);
		e[1] = (int) floor(ipbuff[1]);
		xdet[0] = ipbuff + 2;
		xdet[1] = ipbuff + 4;
		mag = ipbuff + 6;
		for (m = 0; m < mcoefft; m++) {
			modefunc[m] = f(m, e, xdet);
		}
		for (m = 0; m < mcoefft; m++) {
			B[m] += modefunc[m] * (mag[0] - mag[1]);
			for (n = 0; n < mcoefft; n++) {
				A[m][n] += modefunc[m] * modefunc[n];
			}
		}
		nobjects++;
	}
	fprintf(stderr, "# %d pairs read\n", nobjects);

	/* scale the array */
	for (m = 0; m < mcoefft; m++) {
		s[m] = sqrt(A[m][m]);
	}
	for (m = 0; m < mcoefft; m++) {
		B[m] /= s[m];
		for (n = 0; n < mcoefft; n++) {
			A[m][n] /= s[m] * s[n];
		}
	}

	if (outputarray) {
		fprintf(stdout, "A[][]:\n");
		for (m = 0; m < mcoefft; m++) {
			for (n = 0; n < mcoefft; n++) {
				fprintf(stdout, "%14.8lg ", A[m][n]);
			}
			fprintf(stdout, "\n");
		}
		fprintf(stdout, "\nB[]:\n");
		for (m = 0; m < mcoefft; m++) {
			fprintf(stdout, "%14.8lg ", B[m]);
		}
		fprintf(stdout, "\n\ns[]:\n");
		for (m = 0; m < mcoefft; m++) {
			fprintf(stdout, "%14.8lg ", s[m]);
		}
		fprintf(stdout, "\n");
	}

	/* solve linear equations */
	myludcmp(A, mcoefft, indx, &det);
	for (m = 0; m < mcoefft; m++) {
		det *= A[m][m];
	}
	mylubksb(A, mcoefft, indx, B);

	if (outputarray) {
		fprintf(stdout, "\n\ndet = %14.8lg\n", det);
		fprintf(stdout, "\n\nx[]:\n");
		for (m = 0; m < mcoefft; m++) {
			fprintf(stdout, "%14.8lg ", B[m]);
		}
		fprintf(stdout, "\n");
	}

	/* extract coefficients */
	if (fitextinction) {
		for (ee = 1; ee < ne; ee++) {
			dme[ee] = B[ee - 1] / s[ee - 1];
		}
	}

	for (m = 0; m < thelmodel->nmodes; m++) {
		mm = m + fitextinction * (ne - 1);
		((double **) thelmodel->a)[m][0] = B[mm]  / s[mm];
	}

	/* and output */
	if (fitextinction) {
		sprintf(lcstring, "lc -C -n e -n dM > %s", opcatname);
		lcpipe = popen(lcstring, "w");
		if (!lcpipe) {
			error_exit("fitpolymodelmagshift: failed to open lc-pipe for output\n");
		}
		for (ee = 0; ee < ne; ee++) {
			fprintf(lcpipe, "%d %13.8lg\n", ee, dme[ee]);
		}
		pclose(lcpipe);
	}

	oppar = fopen(opparname, "w");
	writelmodel(thelmodel, oppar);

	exit(0);
}


double	f(int m, int *e, double **xdet) {
	double	val = 0;
	int	mm = m - fitextinction * (ne - 1);

	if (fitextinction) {
		if (e[0] == m + 1) {
			val += 1.0;
		}
		if (e[1] == m + 1) {
			val -= 1.0;
		}
	}
	if (mm >= 0) {
		val = lmodelfunc(thelmodel, mm, xdet[0]) - lmodelfunc(thelmodel, mm, xdet[1]);
	}
	return (val);
}
