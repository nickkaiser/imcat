/*
 * getmoid - get velocity distribution at MOID
 *
 * Nick Kaiser 11/25/02
 */

#include <stdio.h>
#include <math.h>
#include "vectors.h"
#include "kepler.h"
#include "planetdata.h"
#include "extravars.h"

#define usage "\nNAME\n\
	getmoid - compute the MOID from Kepler orbit elements\n\
\n\
SYNOPSIS\n\
	getmoid [-extravars vardefs]\n\
\n\
DESCRIPTION\n\
	getmoid reads a catalog containing Keplerian elements\n\
	a, e, i, omega, Omega, M and cartesian r[3], v[3].\n\
\n\
	It then varies M to find the moid and outputs the\n\
	elements as well as the velocity of the particles at MOID.\n\
\n\
	With -extravars option we carry defined variables along.  For example, use\n\
		-extrvars myscalar:1:myvector:3\n\
	to carry along myscalar and myvector[3]\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"

keplerorbit	*theorbit;
double		*r, *v, phi;
double 		oid(double M);
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double));
double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);


main (int argc, char *argv[])
{
	double	*buff, M, dM, Mmin, d, dmin, tol;
	double	ax, bx, cx, fa, fb, fc;
	double	vesc, v1, dv, vperp, torb, Re, sigma, Rc;
	int	NM, buffsize, extravarssize, i;
	FILE	*ipf, *opf;
	char	*iplist, *opdefs, tmpcom[1024];

	/* hard wired parameters */
	NM = 128;
	dM = 2 * M_PI / NM;
	vesc = sqrt(2.0 * G_NEWTON * M_EARTH / R_EARTH);
	v1 = sqrt(G_NEWTON * M_SUN / AU);
	vesc /= v1;

	Re = R_EARTH / AU;

	extravarssize = 0;

	if (argc < 1) {
		fprintf(stderr, usage);
		exit(-1);
	}
	if (argc > 1) {
		if (strncmp("-extravars", argv[1], 10) || (argc < 3)) {
			fprintf(stderr, usage);
			exit(-1);
		}
		parseextravars(argv[2], &extravarssize, &iplist, &opdefs);
	}

	buffsize = 6 + extravarssize;
	buff = (double *) calloc(buffsize, sizeof(double));
	r = (double *) calloc(3, sizeof(double));
	v = (double *) calloc(3, sizeof(double));
	theorbit = (keplerorbit *) calloc(sizeof(keplerorbit));

	/* open pipe to input catalog */
	sprintf(tmpcom, "lc -b -o a e i omega Omega M ");
	if (extravarssize) {
		strcat(tmpcom, iplist);
	}
	ipf = popen(tmpcom, "r");
	if (!ipf) {
		fprintf(stderr, "getmoid: failed to open pipe for input\n");
		exit(-1);
	}

	/* generate output header */
	sprintf(tmpcom, "lc -C -b -N '1 3 r' -N '1 3 v' -n vperp -n dv -n dmin -n R -n a -n e -n i -n omega -n Omega -n M ");
	/* sprintf(tmpcom, "lc -C -b -n a -n e -n i -n omega -n Omega -n M ");*/
	if (extravarssize) {
		strcat(tmpcom, opdefs);
	}
	strcat(tmpcom, " < /dev/null");
	system(tmpcom);

	while(buffsize == fread(buff, sizeof(double), buffsize, ipf)) {
		dmin = 1.e10;
		assignkeplerorbitfrombuffer(theorbit, buff);
		/* find a crude estimate of the minimum */
		dmin = 1.e10;
		for (M = 0; M < 2 * (M_PI + dM); M += dM) {
			d = oid(M);
			if (d < dmin) {
				dmin = d;
				Mmin = M;
			}
		}
		/* now get a refined minimum */
		ax = Mmin - dM;
		bx = Mmin + dM;
		tol = 1.e-10;
		mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, oid);
		dmin = brent(ax, bx, cx, oid, tol, &Mmin);
		theorbit->M = Mmin;
		keplertocartesian(theorbit, r, v);
		phi = atan2(r[1], r[0]);
		rotz(r, -phi);
		rotz(v, -phi);
		v[1] -= 1;
		dv = length(v);	
		vperp = sqrt(v[0] * v[0] + v[2] * v[2]);
		torb = pow(theorbit->a, 1.5);
		sigma = M_PI * Re * Re * (1 + vesc * vesc / (dv * dv));
		Rc = sigma * dv / (4 * M_PI * vperp * torb * 0.05);
		fillbufferfromkeplerorbit(theorbit, buff);
		fwrite(r, sizeof(double), 3, stdout);
		fwrite(v, sizeof(double), 3, stdout);
		fwrite(&vperp, sizeof(double), 1, stdout);
		fwrite(&dv, sizeof(double), 1, stdout);
		fwrite(&dmin, sizeof(double), 1, stdout);
		fwrite(&Rc, sizeof(double), 1, stdout);
/*
fprintf(stderr, "%14.8lg %14.8lg %14.8lg %14.8lg\n", vperp, dv, dmin, Rc);
fprintf(stderr, "%14.8lg %14.8lg %14.8lg %14.8lg %14.8lg %14.8lg\n", r[0], r[1], r[2], v[0], v[1], v[2]);
fprintf(stderr, "%14.8lg %14.8lg %14.8lg %14.8lg %14.8lg %14.8lg %d \n", buff[0], buff[1], buff[2], buff[3], buff[4], buff[5], buffsize);
*/
		fwrite(buff, sizeof(double), buffsize, stdout);
	}

	pclose(ipf);
	exit(0);
}


double oid(double M)
{
	double	phi, d;

	theorbit->M = M;
	keplertocartesian(theorbit, r, v);
	phi = atan2(r[1], r[0]);
	rotz(r, -phi);
	d = sqrt((r[0] - 1.0) * (r[0] - 1.0) + r[2] * r[2]);
	return(d);
}
