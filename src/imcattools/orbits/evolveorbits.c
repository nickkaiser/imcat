#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vectors.h"
#include "tcl.h"
#include "utils/ipbuff.h"
#include "kepler.h"
#include "extravars.h"

#define usage "\nNAME\n\
	evolveorbits - evolve orbits by two-body or time-centered leapfrog\n\
\n\
SYNOPSIS\n\
	evolveorbits dt nsteps [-writekepler] [-tcl nsubsteps] [-extravars vardefs] [-jupiter phase]\n\
\n\
DESCRIPTION\n\
	evolveorbits reads an lc catalog containing a set of positions r[3]\n\
	and velocities v[3] from stdin, evolves the positions and\n\
	velocities through nsteps steps of length and writes the results\n\
	to stdout.\n\
\n\
	By default, it does this by converting from phase-space to\n\
	Kepler orbital elements (a, e, i, omega, Omega, M) and incrementing\n\
	M and then performing the inverse transformation.\n\
\n\
	With -writekepler option, we also output the orbital elements.\n\
\n\
	With -extravars option we carry defined variables along.  For example, use\n\
		-extrvars myscalar:1:myvector:3\n\
	to carry along myscalar and myvector[3]\n\
\n\
	With option -tcl we do the evolution using time-centered\n\
	leapfrog with nsubsteps steps of length dt / nsubsteps\n\
\n\
	Times are given in units of earths dynamical time (approx 58 days).\n\
	Distances are in AU.\n\
\n\
SEE ALSO\n\
	maketestpparfile.pl makeobs_circ laplace3 evolve\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"


main (int argc, char *argv[])
{
	double	*r, *v, *risk, **ipbuff, *keplerbuff, t, dt, jupiterphase;
	keplerorbit	*theorbit;
	int	i, it, nsteps, nsubsteps, np, dotcl, writekepler, arg;
	FILE	*ipf;
	char	*iplist, *opdefs, tmpcom[1024];
	int	extravarssize, ipbuffsize, opbuffsize;

	/* defaults */
	dotcl = 0;
	writekepler = 0;
	iplist = opdefs = NULL;
	ipbuffsize = 6;

	/* parse args */
	if (argc < 3) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[1], "%lf", &dt) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (sscanf(argv[2], "%d", &nsteps) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}
	arg = 3;
	while (arg < argc) {
		if (!strncmp("-writekepler", argv[arg], 12)) {
			writekepler = 1;
		}
		if (!strncmp("-tcl", argv[arg], 4)) {
			dotcl = 1;
			arg++;
			if (arg >= argc) {
				fprintf(stderr, usage);
				exit(-1);
			}
			if (sscanf(argv[arg], "%d", &nsubsteps) != 1) {
				fprintf(stderr, usage);
				exit(-1);
			}
		}
		if (!strncmp("-jupiter", argv[arg], 8)) {
			arg++;
			if (arg >= argc) {
				fprintf(stderr, usage);
				exit(-1);
			}
			if (sscanf(argv[arg], "%lf", &jupiterphase) != 1) {
				fprintf(stderr, usage);
				exit(-1);
			}
			jupiteron(jupiterphase);
		}
		if (!strncmp("-extravars", argv[arg], 10)) {
			arg++;
			if (arg >= argc) {
				fprintf(stderr, usage);
				exit(-1);
			}
			parseextravars(argv[arg], &extravarssize, &iplist, &opdefs);
			ipbuffsize += extravarssize;
		}
		arg++;
	}

	/* open lc pipe for input */
	strcpy(tmpcom, "lc -b -o r v ");
	if (iplist) {
		strcat(tmpcom, iplist);
	}
	ipf = popen(tmpcom, "r");
	if (!ipf) {
		fprintf(stderr, "evolve : failed to open lc-pipe for input\n");
		exit(-1);
	}

	/* send header to stdout */
	strcpy(tmpcom, "lc -C -b -n t -N '1 3 r' -N '1 3 v' ");
	if (opdefs) {
		strcat(tmpcom, opdefs);
	}
	if (writekepler) {
		strcat(tmpcom, "-n a -n e -n i -n omega -n Omega -n M ");
	}
	strcat(tmpcom, " < /dev/null");
	system(tmpcom);

	/* read the data */
	ipbuff = readdoublebuff(ipbuffsize, ipf, &np);

	/* allocate space for the Kepler elements */
	keplerbuff = (double *) calloc(6, sizeof(double));
	theorbit = (keplerorbit *) calloc(1, sizeof(keplerorbit));

	/* fprintf(stderr, "# evolveorbits : read %d particle coords\n", np); */

	for (it = 0; it < nsteps; it++) {
		t = it * dt;
		for (i = 0; i < np; i++) {
			r = ipbuff[i];
			v = ipbuff[i] + 3;
			if (dotcl) {
				tcl(dt / nsubsteps, nsubsteps, r, v);
			} else {
				cartesiantokepler(r, v, theorbit);
				theorbit->M += dt / pow(theorbit->a, 1.5);
				keplertocartesian(theorbit, r, v);
			}
			fwrite(&t, sizeof(double), 1, stdout);
			fwrite(ipbuff[i], sizeof(double), ipbuffsize, stdout);
			if (writekepler) {
				fillbufferfromkeplerorbit(theorbit, keplerbuff);
				fwrite(keplerbuff, sizeof(double), 6, stdout);
			}
		}
	}

	exit(0);
}

