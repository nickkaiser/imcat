#define usage "\n\n\n\
NAME\n\
	simulate --- generate mock fits image\n\
\n\
SYNOPSIS\n\
	simulate [option...]\n\
		-n N1 N2		# size of image (256, 256)\n\
		-S As gs		# amp and slope of stellar lum fun (As = 3e-4; gs = -1.5)\n\
		-G Ag gg		# amp and slope of galaxy lum fun (Ag = 1e-2; gg = -2)\n\
		-R rbar lnrsig		# av and SD log-normal distn of galaxy sizes (rbar = 2.0; lnrsig = 0.3)\n\
		-s sigma		# rms sky fluctuation (0)\n\
		-i seed			# for ran num generator (1)\n\
		-l lmin			# min luminosity (50)\n\
		-e e			# global ellipticity (0.0)\n\
		-m mumin		# min cos theta for galaxy inclination (0.5)\n\
		-c catfile		# output object positions to an lc-format catalogue\n\
\n\
DESCRIPTION\n\
	\"simulate\" generates a mock fits image\n\
	Galaxies and stars have power law lum funs\n\
	Galaxies lave log normal distn of sizes.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "simulate.h"
#include "../imlib/fits.h"
#include "../utils/error.h"
#include "../utils/arrays.h"
#include "../imlib/filters.h"
#include "../utils/ran1.h"


#ifndef PI
#define PI 3.14159
#endif

double	drand48();

FILE	*catfile;
int	outputcat;

int		main(int argc, char *argv[])	
{
	int	arg = 1, N1, N2, ngals, nstars, donoise;
	float	As, gs, Ag, gg, rbar, lnrbar, lnrsig, sigma, lmin, e, mumin;
	fitsheader	*fits;
	float	**f;
	long 	seed;
	char	*catfilename, lcstring[512];
	
	/* defaults */
	N1 = 256;
	N2 = 256;
	As 	= 3.0e-4;
	gs 	= -1.5;
	Ag 	= 1.0e-2;
	gg 	= -2;
	rbar = 2.0; 
	lnrbar = log(rbar);
	lnrsig = 0.3;
	lmin = 50;
	seed = 1;
	e = 0;
	mumin = 0.5;
	donoise = 0;
	outputcat = 0;
	
	while (arg < argc) {
		if (argv[arg][0] != '-')
			error_exit(usage);
		switch(argv[arg++][1]) {
			case 'n':
				if (1 != sscanf(argv[arg++], "%d", &N1))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%d", &N2))
					error_exit(usage);
				break;
			case 'S':
				if (1 != sscanf(argv[arg++], "%f", &As))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%f", &gs))
					error_exit(usage);
				break;
			case 'G':
				if (1 != sscanf(argv[arg++], "%f", &Ag))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%f", &gg))
					error_exit(usage);
				break;
			case 'R':
				if (1 != sscanf(argv[arg++], "%f", &rbar))
					error_exit(usage);
				lnrbar = log(rbar);
				if (1 != sscanf(argv[arg++], "%f", &lnrsig))
					error_exit(usage);
				break;
			case 's':
				if (1 != sscanf(argv[arg++], "%f", &sigma))
					error_exit(usage);
				donoise = 1;
				break;
			case 'i':
				if (1 != sscanf(argv[arg++], "%ld", &seed))
					error_exit(usage);
				break;
			case 'l':
				if (1 != sscanf(argv[arg++], "%f", &lmin))
					error_exit(usage);
				break;
			case 'e':
				if (1 != sscanf(argv[arg++], "%f", &e))
					error_exit(usage);
				break;
			case 'm':
				if (1 != sscanf(argv[arg++], "%f", &mumin))
					error_exit(usage);
				break;
			case 'c':
				outputcat = 1;
				catfilename = argv[arg++];
				break;
			default:
				error_exit(usage);
				break;
		}
	}
	srand48(seed);
	if (outputcat) {
		sprintf(lcstring, "lc -C -N '1 2 x' -n r -n logl > %s", catfilename);
		catfile = popen(lcstring, "w");
	}
	allocFloatArray(&f, N1, N2);	
	ngals = make_pop(Ag, gg, lmin, lnrbar, lnrsig, e, mumin, N1, N2, f);
	nstars = make_pop(As, gs, lmin, lnrbar, -1, e, mumin, N1, N2, f);
	fclose(catfile);
	if (donoise)
		add_noise(f, N1, N2, sigma);
	fits = new2Dfitsheader(N1, N2, FLOAT_PIXTYPE);
	add_comment(argc, argv, fits);
	write2Dfloatimage(f, fits);
	exit(pclose(catfile));
}


void		add_noise(float **f, int N1, int N2, float sigma)
{
	int 		i, j, idum;

	for (i = 0; i < N2; i++)						
		for (j = 0; j < N1; j++)
			f[i][j] += sigma * gasdev();
}





float gasdev(void)
{
	static int iset=0;
	static float gset;
	float fac,r,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*drand48()-1.0;
			v2=2.0*drand48()-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

#define LN_L_MAX	16


/* make_pop() generates a set of exponential disks with dn / dl = A l^g
 * and with log normal distribution of sizes.
 * galaxies are randomly oriented but with effective thickness of 0.2
 * if lnrsig < 0 then stars are generated
 * x-axis is globally stretched by factor 1 + e
 */
int	make_pop(float A, float g, float lmin, float lnrbar, float lnrsig, 
		float e, float mumin, int N1, int N2, float **f)
{
	float 	r, l, r0, f0, x, y, mu, phi;
	int 	i, j, di, dj, dijmax, iobj, n;
	
	n = ceil(- N1 * N2 * A * pow(lmin, 1 + g) / (1 + g));
	for (iobj = 0; iobj < n; iobj++) {
		i = floor(N2 * drand48());
		j = floor(N1 * drand48());
		if (i == N2 || j == N1)
			continue;
		l = lmin * pow(drand48(), 1 / (1 + g));
		if (log(l) > LN_L_MAX)
			continue;
		if (lnrsig > 0) {						/* a galaxy */
			r = exp(lnrbar + lnrsig * gasdev());
			r0 = 1.44 * r;
			mu = mumin + (1 - mumin) * drand48();
			f0 = l / (mu * PI * r0 * r0);
			dijmax = ceil(10 * r0);
			phi = 2 * PI * drand48();
			for (di = -dijmax; di <= dijmax; di++) {
				for (dj = -dijmax; dj <= dijmax; dj++) {
					if (i+di < 0 || i+di >= N2 || j+dj < 0 || j+dj >= N1)
						continue;
					x = di * (1 + e) * cos(phi) + dj * sin(phi);
					y = - di * (1 + e) * sin(phi) + dj * cos(phi);
					f[i+di][j+dj] += f0 * exp(-sqrt(x * x + y * y / (mu * mu)) / r0);
				}
			}
		} else {								/* a star */
			r = 1;
			f[i][j] += l;
		}
		if (outputcat) {		
			fprintf(catfile, "%5d %5d %8.3lf %8.3lf\n", j, i, r, log(l));
		}
	}
	return (n);
}

#undef LN_L_MAX


