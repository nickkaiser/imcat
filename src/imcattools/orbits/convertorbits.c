/*
 * convertorbits - convert between Kepler elements and cartesian coords
 *
 * Nick Kaiser 11/25/02
 */

#include <stdio.h>
#include <math.h>
#include "vectors.h"
#include "kepler.h"
#include "extravars.h"

#define usage "\nNAME\n\
	convertorbits - convert between Kepler elements and phase-space coordinates\n\
\n\
SYNOPSIS\n\
	convertorbits k2c | c2k [-extravars vardefs]\n\
\n\
DESCRIPTION\n\
	convertorbits converts between Keplerian elements\n\
	a, e, i, omega, Omega, M and cartesian r[3], v[3].\n\
\n\
	With -extravars option we carry defined variables along.  For example, use\n\
		-extrvars myscalar:1:myvector:3\n\
	to carry along myscalar and myvector[3]\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"

#define K2C_MODE 0
#define C2K_MODE 1

main (int argc, char *argv[])
{
	double	*buff, *r, *v, *keplerbuff;
	FILE	*ipf, *opf;
	int	opmode, extravarssize, buffsize;
	keplerorbit	*theorbit;
	char	*iplist, *opdefs, tmpcom[1024];

	/* defaults */
	extravarssize = 0;

	if (argc < 2) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (!strcmp(argv[1], "k2c")) {
		opmode = K2C_MODE;
	} else {
		if (!strcmp(argv[1], "c2k")) {
			opmode = C2K_MODE;
		} else {
			fprintf(stderr, usage);
			exit(-1);
		}
	}
	if (argc > 2) {
		if (strncmp("-extravars", argv[2], 10) || (argc < 4)) {
			fprintf(stderr, usage);
			exit(-1);
		}
		parseextravars(argv[3], &extravarssize, &iplist, &opdefs);
	}

	/* allocate buffers */
	buffsize = 6 + extravarssize;
	buff = (double *) calloc(buffsize, sizeof(double));
	theorbit = (keplerorbit *) calloc(sizeof(keplerorbit));
	if (opmode == K2C_MODE) {
		r = (double *) calloc(3, sizeof(double));
		v = (double *) calloc(3, sizeof(double));
	} else {
		keplerbuff = (double *) calloc(6, sizeof(double));
		r = buff;
		v = buff + 3;
	}

	/* open pipe to input catalog */
	if (opmode == K2C_MODE) {
		sprintf(tmpcom, "lc -b -o a e i omega Omega M ");
	} else {
		sprintf(tmpcom, "lc -b -o r v ");
	}
	if (extravarssize) {
		strcat(tmpcom, iplist);
	}
	ipf = popen(tmpcom, "r");
	if (!ipf) {
		fprintf(stderr, "convertorbits: failed to open pipe for input\n");
		exit(-1);
	}

	/* generate output header */
	if (opmode == K2C_MODE) {
		sprintf(tmpcom, "lc -C -b -N '1 3 r' -N '1 3 v' -n a -n e -n i -n omega -n Omega -n M ");
	} else {
		sprintf(tmpcom, "lc -C -b -n a -n e -n i -n omega -n Omega -n M -N '1 3 r' -N '1 3 v' ");
	}
	if (extravarssize) {
		strcat(tmpcom, opdefs);
	}
	strcat(tmpcom, " < /dev/null");
	system(tmpcom);

	while(buffsize == fread(buff, sizeof(double), buffsize, ipf)) {
		if (opmode == K2C_MODE) {
			assignkeplerorbitfrombuffer(theorbit, buff);
			keplertocartesian(theorbit, r, v);
			fwrite(r, sizeof(double), 3, stdout);
			fwrite(v, sizeof(double), 3, stdout);
		} else {
			cartesiantokepler(r, v, theorbit);
			fillbufferfromkeplerorbit(theorbit, keplerbuff);
			fwrite(keplerbuff, sizeof(double), 6, stdout);
		}
		fwrite(buff, sizeof(double), buffsize, stdout);
	}

	pclose(ipf);
	exit(0);
}

