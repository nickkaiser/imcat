/*
 * ppsi.c --- compute the shear power spectrum kernel
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define usage "\n\
NAME\n\
	frw - compute physical quantities in FRW cosmology\n\
\n\
SYNOPSIS\n\
	frw [options....]\n\
	where options are:\n\
		-m Omega_matter	# present Omega in matter (0.99)\n\
		-l Omega_lambda	# present Omega in lambda (0.0)\n\
		-Z Z0		# starting redshift (1.e5)\n\
		-s Zs		# max source redshift (9.0)\n\
		-d dlna		# logarithmic step size (0.01)\n\
		-u		# print this message\n\
\n\
DESCRIPTION\n\
	'frw' numerically integrates Freidmann equation\n\
	to obtain physical quantities:\n\
		z 	= redshift\n\
		y 	= a / a0 = 1 / (1 + z)\n\
		h 	= H / H0 = sqrt(Omegam / y^3 + Omegal + (1 - Omega0) / y^2)\n\
		eta	= conformal time\n\
		zeta 	= eta0 - eta = conformal lookback time\n\
		d 	= solution of ddotdot + 2 h ddot - (3/2) Omegam d / y^3 = 0\n\
		dphi 	= d / y = potential growth factor;\n\
		r	= radial comoving distance\n\
		D	= comoving angular diameter distance\n\
		a	= scale factor\n\
	where\n\
		Omega0 = Omega_matter + Omega_lambda and\n\
	We start integrating at very high Z0 in order to get into pure\n\
	growing mode, but only output for Z < Zs.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoroto.ca\n\
\n\n"


main (int argc, char *argv[]) 
{
	double	Omegam, Omegal, Omega0, Z0, Zmax, dlna;
	double	y1, y2, y, z, zeta, eta, eta0, deta, h, d, ddot, ddotdot, dphi, a0;
	double	*D, *Dphi, *Y, *H, *Eta, drdz, dr, r, dz, s, t, *T;
	int	arg = 1, nalloc, nreal, i;
	FILE	*opf;
	char	opcom[512];

	/* defaults */
	Omegam = 0.99;
	Omegal = 0.0;
	Z0 = 1.e5;
	Zmax = 9.0;
	dlna = 0.01;
	
	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			fprintf(stderr, usage);
			exit(-1);
		}
		switch (argv[arg++][1]) {
			case 'm':
				sscanf(argv[arg++], "%lf", &Omegam);
				break;
			case 'l':
				sscanf(argv[arg++], "%lf", &Omegal);
				break;
			case 'Z':
				sscanf(argv[arg++], "%lf", &Z0);
				break;
			case 's':
				sscanf(argv[arg++], "%lf", &Zmax);
				break;
			case 'd':
				sscanf(argv[arg++], "%lf", &dlna);
				break;
			case 'u':
			default:
				fprintf(stderr, usage);
				exit(-1);
				break;
		}
	}
	Omega0 = Omegam + Omegal;

	nalloc = 1.2 * log(1 + Z0) / dlna;
	D 	= (double *) calloc(nalloc, sizeof(double));
	Dphi 	= (double *) calloc(nalloc, sizeof(double));
	Y 	= (double *) calloc(nalloc, sizeof(double));
	H 	= (double *) calloc(nalloc, sizeof(double));
	Eta 	= (double *) calloc(nalloc, sizeof(double));
	T 	= (double *) calloc(nalloc, sizeof(double));
	if (!D || !Dphi || !Y || !H || !Eta || !T) {
		fprintf(stderr, "memory allocation failed!\n");
		exit(-1);
	}

	/* initial conditions... */
	y1 = 1.0 / (1 + Z0);
	y2 = 1.0;
	eta = 0.0;
	ddot = 0.0;
	d = 1.0;
	t = 0.0;
	i = 0;
	for (y = y1; y < y2; y *= (1 + dlna)) {
		h = Omegam / (y * y * y) + Omegal + (1 - Omega0) / (y * y);
		if (h < 0.0) {
			fprintf(stderr, "bad Omega values!\n");
			exit(-1);
		}
		h = sqrt(h);
		deta = sqrt(1 - Omega0) * dlna / (h * y);
		eta += deta;
		t += y * deta;
		ddotdot = 1.5 * d * Omegam / (y * y * y) - 2 * h * ddot;
		d += ddot * dlna / h;
		ddot += ddotdot * dlna / h;
		dphi = d / y;
		if (i >= nalloc) {
			fprintf(stderr, "overstepping allocated memory!\n");
			exit(-1);
		}
		D[i] 	= d;
		Dphi[i] = dphi;
		Y[i] 	= y;
		H[i]	= h;
		Eta[i]	= eta;
		T[i]	= t;
		i++;
	}
	eta0 = eta;
	a0 = 3000.0 / sqrt(1 - Omega0);
	s = 1.0 / a0;

	nreal = i;
	/* now output */
	sprintf(opcom, "lc -C -n z -n y -n dphi -n eta -n zeta -n r -n dz -n drdz -n a -n D -n ct -x -a 'frw' \
		-H 'Omega_matter = %14.8lg' -H 'Omega_lambda = %14.8lg' -H 'a0 = %14.8lg' -H 'eta0 = %14.8lg'", 
		Omegam, Omegal, a0, eta0);
	opf = popen(opcom, "w");
	for (i = nreal - 1; i >= 0; i--) {
		if (Y[i] < 1.0 / (1 + Zmax)) {
			exit(0);
		}
		zeta = eta0 - Eta[i];
		r = zeta / s;
		if (i == nreal - 1) {
			dr = 0.0;
			dz = 0.0;
			drdz = 0.0;
		} else {
			dr = (Eta[i+1] - Eta[i]) / s;
			dz = 1.0 / Y[i] - 1.0 / Y[i+1];
			drdz = dr / dz;
		}
		fprintf(opf, "%14lg %14lg %14lg %14lg %14lg %14lg %14lg %14lg %14lg %14lg %14lg\n", 
			1.0 / Y[i] - 1.0, Y[i], Dphi[i] / dphi, Eta[i], zeta, r, dz, drdz, a0 * Y[i], sinh(zeta) / s, a0 * T[i]);
	}
	pclose(opf);

	exit(0);
}


