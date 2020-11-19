/*
 * nfw.c - compute nfw density profile
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils/args.h"


#define usage "\nNAME\n\
	nfw - compute Navarro-Frenk-White profile\n\
\n\
SYNOPSIS\n\
	nfw [-u] [-v v200] [-x x1 x2 nx] [-z z0] [-c c] [-l]\n\
\n\
DESCRIPTION\n\
	nfw computes nfw density profile. It outputs a table\n\
	(in lc format if IMCAT is installed) containing:\n\
		x	# dimensionless radius = r/r200\n\
		v	# circular velocity\n\
		r	# physical radius in cm\n\
		dr	# physical delta-radius in cm\n\
		m	# mass in M_solar\n\
		rho	# physical density in g/cm^3\n\
		sigma	# physical projected density in g/cm^2\n\
		theta	# angle in arcseconds\n\
		kappa	# dimensionless surface density\n\
\n\
	Options are:\n\
		-u	# print this message\n\
		-v v200	# circular velocity at r_200 in km/s (200)\n\
		-x x1 x2	# min and max dimensionless radii (0.01, 100)\n\
		-z z0	# lens redshift (0.2)\n\
		-c c	# compactness parameter (10.0)\n\
		-l	# use pure lambda cosmology (otherwise use EdeS)\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"

#define	EDS	0
#define LAMBDA	1

/* cgs units */
#define	G	6.67e-8
#define	PC	3.086e18
#define	KPC	3.086e21
#define H0	(10.0/PC)
#define	C	3.e10
#define ARCSEC	(M_PI/(180*60*60))
#define KM	1.e5
#define	MSUN	1.99e33	

main (int argc, char *argv[])
{
	char	*flag, comstring[512];
	double	c, *v, v200, *x, x1, x2, dlnx, z0, *r, *dr, r200, *rho, *sigma, b, *kappa;
	double	*kappasum, *areasum, *kappabar, *gamma, *m;
	double	c1, cx, cx1, deltac, omegafac, rhofac, Da, sigmacrit, dA;
	int	i, j, nx, cosmo;
	FILE	*opf;

	/* defaults */
	v200 = 220.0;
	x1 = 0.01;
	x2 = 10.0;
	nx = 100;
	z0 = 0.2;
	c  = 10.0;
	cosmo = EDS;

	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'v':
				v200 = getargd();
				break;
			case 'x':
				x1 = getargd();
				x2 = getargd();
				nx = getargi();
				break;
			case 'z':
				z0 = getargd();
				break;
			case 'c':
				c = getargd();
				break;
			case 'l':
				cosmo = LAMBDA;
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}

	/* some constants */
	dlnx = log(x2 / x1) / nx;
	c1 = c + 1;
	deltac = 200 * c * c * c / (3 * (log(c1) - c / c1));

	/* omegafac = Omega_m(0) (1 + z0)^3 / Omega_m(z_0) */
	switch (cosmo) {
		case EDS:
			omegafac = pow(1 + z0, 3.0);
			Da = 2 * C * (1.0 - 1.0 / sqrt(1 + z0)) / (H0 * (1 + z0));
			break;
		case LAMBDA:
			omegafac = 1.0;
			Da = C * z0 / (H0 * (1 + z0));
			break;
		default:
			error_exit("frw: illegal cosmology\n");
	}
	r200 	= v200 * KPC / sqrt(omegafac);
	rhofac	= 3 * H0 * H0 * omegafac / (8 * M_PI * G);
	sigmacrit = C * C / (4 * M_PI * G * Da);	

	/* open opf */
	if (getenv("IMCATDIR")) {
		sprintf(comstring, "lc -C -n x -n v -n r -n dr -n m -n rho -n sigma -n theta -n kappa -n kappabar -n gamma -H 'c = %14.8lg' -H 'v200 = %14.8lg' -H 'z0 = %14.8lg'", c, v200, z0);
		opf = popen(comstring, "w");
	} else {
		opf = stdout;
		fprintf(opf, "# nfw output\n#\tc\t= %14.8lg\n#\tv200\t= %14.8lg\n#\tz0\t= %14.8lg\n", c, v200, z0);
		fprintf(opf, "#             x              v              r             dr              m            rho");
		fprintf(opf, "          sigma          theta          kappa       kappabar          gamma\n");
	}

	/* allocate arrays */
	x 	= (double *) calloc(nx, sizeof(double));
	m 	= (double *) calloc(nx, sizeof(double));
	r 	= (double *) calloc(nx, sizeof(double));
	dr 	= (double *) calloc(nx, sizeof(double));
	rho 	= (double *) calloc(nx, sizeof(double));
	v 	= (double *) calloc(nx, sizeof(double));
	sigma 	= (double *) calloc(nx, sizeof(double));
	kappa 	= (double *) calloc(nx, sizeof(double));
	kappasum 	= (double *) calloc(nx, sizeof(double));
	areasum 	= (double *) calloc(nx, sizeof(double));
	kappabar	= (double *) calloc(nx, sizeof(double));
	gamma	= (double *) calloc(nx, sizeof(double));

	/* compute x, r, dr, v, rho, m */
	for (i = 0; i < nx; i++) {
		x[i] = x1 * exp(i * dlnx);
		r[i] = x[i] * r200;
		if (i) {
			dr[i] = r[i] - r[i - 1];
		}
		cx = c * x[i];
		cx1 = cx + 1.0;
		v[i] = v200 * sqrt((log(cx1) - cx / cx1) / (x[i] * (log(c1) - c / c1)));
		m[i] = v[i] * v[i] * KM * KM * r[i] / (MSUN * G);
		rho[i] = rhofac * deltac / (cx * cx1 * cx1);
	}

	/* compute surface density sigma */
	for (i = 0; i < nx; i++) {
		b = r[i];
		for (j = i + 1; j < nx; j++) {
			sigma[i] += dr[i] * rho[i] / sqrt(1 - b * b / (r[j] * r[j]));
		}
	}
	sigma[0] = sigma[1];
	sigma[nx - 1] = sigma[nx - 2] * sigma[nx - 2] / sigma[nx - 3];

	/* compute dimensionless surface density kappa and kappa, area integrals */
	for (i = 0; i < nx; i++) {
		kappa[i] = sigma[i] / sigmacrit;
		if (i) {
			dA = 2 * M_PI * dr[i] * r[i];
			kappasum[i] = kappasum[i-1] + kappa[i] * dA;
			areasum[i] = areasum[i-1] + dA;
		}
	}

	/* compute kappabar */
	kappabar[0] = kappa[0];
	for (i = 1; i < nx; i++) {
		kappabar[i] = kappasum[i] / areasum[i];
	}

	/* compute gamma */
	for (i = 2; i < nx; i++) {
		gamma[i] = - (kappabar[i] - kappabar[i-1]) / log(areasum[i] / areasum[i-1]);
	}
	gamma[0] = gamma[1] = gamma[2];

	/* output */
	for (i = 0; i < nx; i++) {
		fprintf(opf, " %14.8lg %14.8lg %14.8lg %14.8lg %14.8lg %14.8lg %14.8lg %14.8lg %14.8lg %14.8lg %14.8lg\n", 
			x[i], v[i], r[i], dr[i], m[i], rho[i], sigma[i], r[i] / (Da * ARCSEC), kappa[i], kappabar[i], gamma[i]); 
	}

	pclose(opf);
	exit(0);
}
