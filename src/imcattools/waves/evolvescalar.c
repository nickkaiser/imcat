/* to evolve an interacting scalar field */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils/args.h"
#include "utils/arrays.h"
#include "utils/error.h"
#include "imlib/fits.h"

/* defaults */
#define NSTEPS 	8
#define KSTAR	0.3
#define DT  	0.25
#define LAMBDA  0.0
#define	GAMMA	0.0

#define swap(a,b) tmp=(a);(b)=(a);(a)=tmp
#define MAX_FRAMES 99999

#define USAGE "\nNAME\n\
	evolvescalar - evolve a 2-D scalar field\n\
\n\
SYNOPSIS\n\
	evolvescalar [-nsteps nsteps (8)] [-nframes nframes (99999)]\n\
		[-kstar kstar (0.3)] [-autoscale] [-dt dt (0.25)]\n\
		[-lambda lambda (0.0)] [-gamma gamma (0.0)] [-u]\n\
		[-gravity V0 sigma]\n\
\n\
DESCRIPTION\n\
	Evolvescalar reads from stdin a 3-D FITS file f[2][Ny][Nx] consisting\n\
	of 2 planes containing the initial field d[Ny][Nx] and the initial\n\
	field velocity v[Ny][Nx].\n\
	It then evolves the coupled equations\n\
		dv/dt = laplacian(d) - kstar^2 d - 4 lambda d^3 - gamma v\n\
		dd/dt = v\n\
	With gamma = 0 these are equivalent to the Klein-Gordon equation\n\
	with mass ~ kstar and a lambda phi^4 self-interaction.\n\
\n\
	With -gravity option we replace kstar by\n\
	kstar (1 - V0 exp(-0.5 r**2 / sigma**2)).  This simulates the effect of\n\
	a fixed external gravitational potential.\n\
\n\
	With non-zero damping coefficient gamma, the waves evolve as\n\
	in an expanding universe with H = gamma / 3.\n\
\n\
	If kstar is negative, it evolves\n\
		dv/dt = laplacian(d) + kstar^2 d - 4 lambda d^3 - gamma v\n\
		dd/dt = v\n\
	which has a negative mass term and gives a 'w'-shaped potential.\n\
\n\
	The laplacian function is computed as\n\
\n\
	laplacian = d[y][x-1] + d[y][x+1] + d[y-1][x] + d[y+1][x] - 4 * d[y][x].\n\
\n\
	The evolution scheme is simple time-centred leapfrog.  I.e. we\n\
	update v -> v + vdot * dt and then update d -> d + v * dt\n\
\n\
	It outputs a 3-D FITS file f[nframes][Ny][Nx] containing the\n\
	evolved field.\n\
\n\
OPTIONS\n\
		-nsteps nsteps		# number of steps between output frames\n\
		-nframes nframes	# total number of output frames\n\
		-kstar kstar		# Compton wave-number\n\
		-autoscale		# scale each output frame to 0-256\n\
		-dt dt			# time-step\n\
		-lambda lambda		# interaction strength\n\
		-gamma gamma		# damping rate\n\
		-u			# print man-page\n\
\n\
SEE ALSO\n\
	edw, generate_dw\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"



static	fitsheader	*fits;

void	outputframe(float **d, float **dout, int Nx, int Ny, int autoscale, int i);

int	main (int argc, char* argv[]) 
{
	int		arg, nsteps, nframes, Nx, Ny, x, y, xp, yp, xm, ym, f, i, autoscale;
	char		*flag;
	double		dt, dx, kstar, kstarsquared, laplacian, lambda, gamma;
	float		**d, **v, **dout, **V;
	int		dogravity;
	double		V0, sigma;

	/* set defaults */
	autoscale = 0;
	nsteps 	= NSTEPS;
	kstar 	= KSTAR;
	dx 	= 1.0;
	dt 	= DT;
	lambda 	= LAMBDA;
	nframes = MAX_FRAMES;
	gamma	= GAMMA;
	dogravity = 0;

	/* parse args */
	argsinit(argc, argv, USAGE);
	while ((flag = getflag())) {
		if (!strncmp(flag, "nsteps", 6)) {
			nsteps = getargi();
		} else {
			if (!strncmp(flag, "kstar", 5)) {
				kstar = getargd();
			} else {
				if (!strncmp(flag, "autoscale", 9)) {
					autoscale = 1;
				} else {
					if (!strncmp(flag, "dt", 2)) {
						dt = getargd();
					} else {
						if (!strncmp(flag, "lambda", 6)) {
							lambda = getargd();
						} else {
							if (!strncmp(flag, "nframes", 7)) {
								nframes = getargi();
							} else {
								if (!strncmp(flag, "gamma", 5)) {
									gamma = getargd();
								} else {
									if (!strncmp(flag, "gravity", 7)) {
										dogravity = 1;
										V0 =  getargd();
										sigma =  getargd();
									} else {
										error_exit(USAGE);
									}
								}
							}
						}
					}
				}
			}
		}
	}

	/* compute mass-squared coefficient */
	kstarsquared = (kstar > 0 ? kstar * kstar : - kstar * kstar);
	
	/* read the data header */
	fits = readfitsheader(stdin);
	if (fits->ndim != 3 || fits->n[2] != 2) {
		error_exit("evolvescalar: source data must be a Nx x Ny x 2 FITS image\n");
	}
	Nx = fits->n[0];
	Ny = fits->n[1];
	if (nframes > 1) {
		fits->n[2] = nframes + 1;
	} else {
		fits->ndim = 2;
	}

	/* allocate the memory */
	allocFloatArray(&d, Nx, Ny);
	allocFloatArray(&v, Nx, Ny);
	if (autoscale) {
		allocFloatArray(&dout, Nx, Ny);
	} else {
		dout = d;
	}
	/* read the initial data */
	readfitsplane((void **) d, fits);
	readfitsplane((void **) v, fits);

	/* output header */
	add_comment(argc, argv, fits);
	writefitsheader(fits);

	/* output initial frame */
	if (nframes > 1) {
		outputframe(d, dout, Nx, Ny, autoscale, 0);
	}

	/* create the gravitational potential */
	if (dogravity) {
		allocFloatArray(&V, Nx, Ny);
		for (y = 0; y < Ny; y++) {
			for (x = 0; x < Nx; x++) {
				V[y][x] = V0 * exp(-0.5 * ((x - Nx / 2) * (x - Nx / 2) + (y - Ny / 2) * (y - Ny / 2))/ ( sigma * sigma));
			}
		}
	}

	/* evolve */
	for (f = 0; f < nframes; f++) {
		for (i = 0; i < nsteps; i++) {
			for (y = 0; y < Ny; y++) {
				for (x = 0; x < Nx; x++) {
					d[y][x] += v[y][x] * dt;
				}
			}							
			for (y = 0; y < Ny; y++) {
				ym = (Ny + y - 1) % Ny;
				yp = (y + 1) % Ny;
				for (x = 0; x < Nx; x++) {
					xm = (Nx + x - 1) % Nx;
					xp = (x + 1) % Nx;
					laplacian = d[y][xm] + d[y][xp] + d[ym][x] + d[yp][x] - 4 * d[y][x];
					v[y][x] += (laplacian - kstarsquared * d[y][x] - 4 * lambda * pow(d[y][x], 3.0) - gamma * v[y][x]) * dt;
				}
			}
		}
		outputframe(d, dout, Nx, Ny, autoscale, f + 1);
	}

	writefitstail(fits);

	exit(0);
}


void	outputframe(float **d, float **dout, int Nx, int Ny, int autoscale, int i)
{
	int	x, y;
	double	dmin, dmax, ddsum, dsum, dbar, sigmad;

	fprintf(stderr, "# i = %d ", i);
	if (autoscale) {
		dmin = dmax = d[0][0];
		ddsum = dsum = 0.0;
		/* get the scaling */
		for (y = 0; y < Ny; y++) {
			for (x = 0; x < Nx; x++) {
				dmin 	= (d[y][x] < dmin ? d[y][x] : dmin);
				dmax 	= (d[y][x] > dmax ? d[y][x] : dmax);
				dsum 	+= d[y][x];
				ddsum 	+= d[y][x] * d[y][x];
			}
		}
		dbar 	= dsum / (Nx * Ny);
		sigmad 	= sqrt(ddsum / (Nx * Ny) - dbar * dbar);
		fprintf(stderr, "\tdmin=%14.8lg\tdmax=%14.8lg\tdbar=%14.8lg\tsigmad=%14.8lg ", 
			dmin, dmax, dbar, sigmad);
		dmin = dbar - 1.5 * sigmad;
		dmax = dbar + 1.5 * sigmad;
		for (y = 0; y < Ny; y++) {
			for (x = 0; x < Nx; x++) {
				dout[y][x] = 256 * (d[y][x] - dmin) / (dmax - dmin);
			}
		}
	}
	fprintf(stderr, "\n");
	writefitsplane((void **) dout, fits);
}
