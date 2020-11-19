#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"
#include "tcl.h"
#include "utils/ipbuff.h"

#define usage "\nNAME\n\
	tcl_evolveN\n\
\n\
SYNOPSIS\n\
	tcl_evolveN dt nsteps nrep\n\
\n\
DESCRIPTION\n\
	tcl_evolve reads an lc catalog containing a set of positions r[3]\n\
	and velocities v[3] from stdin, evolves the positions and\n\
	velocities through nrep repetitions of nsteps steps of length dt\n\
	using time centered leapfrog, and writes the results to stdout.\n\
\n\
	Times are given in units of earths dynamical time (approx 58 days).\n\
	Distances are in AU.\n\
\n\
SEE ALSO\n\
	maketestpparfile.pl makeobs_circ laplace3 tcl_evolve\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"


main (int argc, char *argv[])
{
	double	*r, *v, *risk, **buff;
	double	dt;
	int	i, it, nsteps, nrep, np;
	FILE	*ipf;

	/* parse args */
	if (argc != 4) {
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
	if (sscanf(argv[3], "%d", &nrep) != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}

	/* open lc pipe for input */
	ipf = popen("lc -b -o r v risk v_inf", "r");
	if (!ipf) {
		fprintf(stderr, "tcl_evolve : failed to open lc-pipe for input\n");
		exit(-1);
	}

	/* send header to stdout */
	system("lc -C -b -N '1 3 r' -N '1 3 v' -n risk -n v_inf < /dev/null");

	/* read the data */
	buff = readdoublebuff(8, ipf, &np);

	fprintf(stderr, "# tcl_evolveN : read %d particle coords\n", np);

	while (nrep--) {
		for (i = 0; i < np; i++) {
			r = buff[i];
			v = buff[i] + 3;
			tcl(dt, nsteps, r, v);
			fwrite(buff[i], sizeof(double), 8, stdout);
		}
	}

	exit(0);
}

