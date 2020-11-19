#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"
#include "deltam.h"

#define usage "\nNAME\n\
	getdet - keeps track of the detections of a set of objects\n\
\n\
SYNOPSIS\n\
	getdet np\n\
\n\
DESCRIPTION\n\
	getdet reads from stdin an lc catalog containing variables\n\
	'pnum' (particle number) and 'fnum' (frame number) and outputs\n\
	a sequence of fnum's and the count 'ndet' of distinct particles that\n\
	have appeared in the input stream.\n\
\n\
	It also outputs 'risksum' - the sum of the risk values for\n\
	the detected particles.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"


main (int argc, char *argv[])
{
	double	oprec[3], iprec[3], *risk, therisk, risksum;
	int	ip, i, np, *det, ndet, fnum, pnum, lastfnum;
	FILE	*ipf;

	/* parse args */
	if (argc != 2) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[1], "%d", &np) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}

	/* allocate space for det[np] */
	det = (int *) calloc(np, sizeof(int));
	/* and for risk[np] */
	risk = (double *) calloc(np, sizeof(double));

	/* open lc pipe for input */
	ipf = popen("lc -b -o fnum pnum risk", "r");
	if (!ipf) {
		fprintf(stderr, "getdet : failed to open lc-pipe for input\n");
		exit(-1);
	}

	/* send the cat header to stdout */
	system("lc -C -b -n fnum -n ndet -n risksum < /dev/null");

	/* initialize */
	lastfnum = -1;

	/* process the input */
	while (3 == fread(iprec, sizeof(double), 3, ipf)) {
		fnum = (int) iprec[0];
		pnum = (int) iprec[1];
		therisk = iprec[2];
		if (fnum != lastfnum) {
			ndet = 0;
			risksum = 0.0;
			for (ip = 0; ip < np; ip++) {
				ndet += det[ip];
				risksum += risk[ip];
			}
			oprec[0] = (double) fnum;
			oprec[1] = (double) ndet;
			oprec[2] = risksum;
			fwrite(oprec, sizeof(double), 3, stdout);
			lastfnum = fnum;
		}
		if (pnum < 0 || pnum >= np) {
			fprintf(stderr, "getdet: pnum value out of allowed bounds\n");
			exit(-1);
		}
		det[pnum] = 1;
		risk[pnum] = therisk;
	}

	exit(0);
}

