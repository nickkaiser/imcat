/* to evolve an interacting scalar field */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils/args.h"
#include "utils/arrays.h"
#include "utils/error.h"
#include "imlib/fits.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>


/* choose precision for floating point numbers
#define FPN_FLOAT 1 */
#if FPN_FLOAT
#define FPN	float
#define FITSPIXTYPE -32
#else
#define FPN	double
#define FITSPIXTYPE -64
#endif

#define LEAPFROGMODE 	1
#define RUNGEKUTTAMODE	2

#define swap(a,b) tmp=(a);(b)=(a);(a)=tmp
#define MAX_FRAMES 99999

#define USAGE "\nNAME\n\
	evolverealscalar - evolve a 2-D scalar field\n\
\n\
SYNOPSIS\n\
	evolverealscalar [-nsteps nsteps (8)] [-nframes nframes (99999)]\n\
		[-m m (0.3)] [-autoscale] [-dt dt (0.25)]\n\
		[-lambda lambda d0 (0.0,0.0)] [-gamma gamma (0.0)] [-u]\n\
		[-gravity V0 sigma] [-RK] [-e epsilon (1e-6)]\n\
\n\
DESCRIPTION\n\
	Evolverealscalar reads from stdin a 3-D FITS file f[2][Ny][Nx] consisting\n\
	of 2 planes containing the initial field d[Ny][Nx] and the initial\n\
	field velocity v[Ny][Nx].\n\
\n\
	It then evolves the coupled equations\n\
\n\
		dv/dt = laplacian(d) - m^2 d - 4 lambda (d - d0)^2 d - gamma v\n\
		dd/dt = v\n\
\n\
	where the laplacian function is computed as\n\
\n\
	laplacian = d[y][x-1] + d[y][x+1] + d[y-1][x] + d[y+1][x] - 4 * d[y][x].\n\
\n\
	With gamma = 0 these are the Klein-Gordon equation with mass m\n\
	and optional lambda phi^4 self-interaction or w-shaped potential with d0.\n\
\n\
	With -gravity option we replace m by\n\
	m (1 - V0 exp(-0.5 r**2 / sigma**2)).  This simulates the effect of\n\
	a fixed external gravitational potential.\n\
\n\
	With non-zero damping coefficient gamma, the waves evolve as\n\
	in an expanding universe with H = gamma / 3.\n\
\n\
	The evolution scheme is simple time-centred leapfrog, in which we update\n\
	v -> v + vdot * dt and then update d -> d + v * dt.\n\
\n\
	Use -RK option to use Runge-Kutta (gamma, lambda and gravity not\n\
	implemented yet).\n\
\n\
	It outputs a 3-D FITS file f[nframes][Ny][Nx] containing the\n\
	evolved field.\n\
\n\
OPTIONS\n\
		-nsteps nsteps		# number of steps between output frames\n\
		-nframes nframes	# total number of output frames\n\
		-m m			# mass (Compton wave-number)\n\
		-dt dt			# time-step\n\
		-lambda lambda		# self interaction strength\n\
		-gamma Gamma		# damping rate\n\
		-RK			# use Runge-Kutta integration\n\
		-e epsilon		# Runge-Kutta tolerance parameter\n\
		-u			# print man-page\n\
\n\
SEE ALSO\n\
	edw, generate_dw, evolvescalar, evovecomplexscalar\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"



static	fitsheader	*fits;
static	FPN		**d, **v, **V;
static	int		Nx, Ny, dogravity, dolambda, dogamma;
static	double		Gamma, lambda, d0, msquared;

/* function for Runge Kutta integration */
int	func (double t, const double dvdata[], double dvdot[], void *params)
{
	int	x, xp, xm, y, yp, ym;
	double	laplacian;

  	(void)(t); 		/* to avoid unused parameter warning */
	(void)(params);

	for (y = 0; y < Ny; y++) {
		ym = (Ny + y - 1) % Ny;
		yp = (y + 1) % Ny;
		for (x = 0; x < Nx; x++) {
			xm = (Nx + x - 1) % Nx;
			xp = (x + 1) % Nx;
			/* rate of change of d */
			dvdot[y * Nx + x] = dvdata[(Ny + y) * Nx + x];
			/* rate of change of v */
			laplacian = d[y][xm] + d[y][xp] + d[ym][x] + d[yp][x] - 4 * d[y][x];
			dvdot[(Ny + y) * Nx + x] = laplacian - msquared * d[y][x];
		}
	}							
 	return GSL_SUCCESS;
}

void    allocRealArray(FPN ***f, int N1, int N2);

int	main (int argc, char* argv[]) 
{
	int		nsteps, nframes, x, y, xp, yp, xm, ym, f, i, intmode;
	char		*flag;
	double		dt, dx, m, laplacian, params[10];
	FPN		*dvdata, *dvdot;
	double		V0, sigma;
	/* GSL stuff */
	gsl_odeiv2_driver 	*driver;
	double			t, ti, epsilon;

	/* set defaults */
	nsteps 	= 8;
	m 	= 0.3;
	dx 	= 1.0;
	dt 	= 0.25;
	nframes = 99999;
	Gamma	= 0.0;
	lambda 	= 0.0;
	d0 	= 0.0;
	dogravity = 0;
	dolambda = 0;
	dogamma = 0;
	epsilon	= 1e-6;
	intmode = LEAPFROGMODE;

	/* parse args */
	argsinit(argc, argv, USAGE);
	while ((flag = getflag())) {
		if (!strncmp(flag, "nsteps", 6)) {
			nsteps = getargi();
		} else {
			if (!strncmp(flag, "m", 5)) {
				m = getargd();
			} else {
				if (!strncmp(flag, "dt", 2)) {
					dt = getargd();
				} else {
					if (!strncmp(flag, "lambda", 6)) {
						lambda = getargd();
						d0 = getargd();
						dolambda = 1;
					} else {
						if (!strncmp(flag, "nframes", 7)) {
							nframes = getargi();
						} else {
							if (!strncmp(flag, "gamma", 5)) {
								Gamma = getargd();
								dogamma = 1;
							} else {
								if (!strncmp(flag, "gravity", 7)) {
									dogravity = 1;
									V0 =  getargd();
									sigma =  getargd();
								} else {
									if (!strncmp(flag, "e", 1)) {
										epsilon = getargd();
									} else {
										if (!strncmp(flag, "RK", 2)) {
											intmode = RUNGEKUTTAMODE;
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
	}

	/* compute mass-squared coefficient */
	msquared = m * m;
	
	/* read the data header */
	fits = readfitsheader(stdin);
	if (fits->ndim != 3 || fits->n[2] != 2) {
		error_exit("evolvescalar: source data must be a f[2][Ny][Nx] FITS image\n");
	}
	Nx = fits->n[0];
	Ny = fits->n[1];
	fits->intpixtype = FITSPIXTYPE;
	fits->n[2] = nframes;

	/* allocate the memory */
	if (intmode == LEAPFROGMODE) {
		allocRealArray(&d, Nx, Ny);
		allocRealArray(&v, Nx, Ny); 
	} else {
		/* initialisation for Runge Kutta */
		gsl_odeiv2_system sys = {func, NULL, 2 * Nx * Ny, params};
  		driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, epsilon, epsilon, 0.0);
		/* global pointers to the columns we can dereference as 2D arrays */
        	d = (FPN **) calloc(Ny, sizeof(FPN *));
		v = (FPN **) calloc(Ny, sizeof(FPN *));
        	if ((!d) || (!v)) 
        	        error_exit("evolvescalarRK: memory allocation failure\n");
		/* pointer to the actual memory for d and v */
		dvdata = (FPN *) calloc(2 * Ny * Nx, sizeof(FPN));
		dvdot  = (FPN *) calloc(2 * Ny * Nx, sizeof(FPN));
         	if ((!dvdata) || (!dvdot))
        	        error_exit("evolvescalarRK: memory allocation failure\n");
        	for (i = 0; i < Ny; i++) {
        	        d[i] = dvdata + i * Nx;
			v[i] = dvdata + (Ny + i) * Nx;
		}
	}

	/* read the initial data */
	readfitsplane((void **) d, fits);
	readfitsplane((void **) v, fits);

	/* output header */
	add_comment(argc, argv, fits);
	writefitsheader(fits);

	/* create the gravitational potential */
	if (dogravity) {
		allocRealArray(&V, Nx, Ny);
		for (y = 0; y < Ny; y++) {
			for (x = 0; x < Nx; x++) {
				V[y][x] = V0 * exp(-0.5 * ((x - Nx / 2) * (x - Nx / 2) + (y - Ny / 2) * (y - Ny / 2))/ ( sigma * sigma));
			}
		}
	}

	/* evolve */
	for (f = 0; f < nframes; f++) {
		for (i = 0; i < nsteps; i++) {
		    if (intmode == RUNGEKUTTAMODE) {
			ti = (f * nsteps + (i + 1)) * dt;
			/* fprintf(stderr, "# t = %10.5f ti = %10.5f\n", t, ti); */
	      		int status = gsl_odeiv2_driver_apply (driver, &t, ti, dvdata);
			if (status != GSL_SUCCESS){
          			printf ("error, GSL return value=%d\n", status);
          			break;
        		}
		    } else {
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
					v[y][x] += (laplacian 
						- (dogravity ? msquared * (1 - V[y][x]) * (1 - V[y][x]) : msquared) * d[y][x] 
						- (dolambda ? 4 * lambda * d[y][x] * (d[y][x] * d[y][x] - d0 * d0) : 0) 
						- (dogamma ? Gamma * v[y][x] : 0)) * dt;
				}
			}
		    }
		}
		fprintf(stderr, "# frame %4d t = %10.3f\n", f, t);
		writefitsplane((void **) d, fits);
	}

	writefitstail(fits);

	exit(0);
}


void    allocRealArray(FPN ***f, int N1, int N2)
{
        int             i;
        
        (*f) = (FPN **) calloc(N2, sizeof(FPN *));
        if (!*f)
                error_exit("allocRealArray: memory allocation failure\n");
        (*f)[0] = (FPN *) calloc(N1 * N2, sizeof(FPN));
        if (!(*f)[0])
                error_exit("allocRealArray: memory allocation failure\n");
        for (i = 1; i < N2; i++)
                (*f)[i] = (*f)[i - 1] + N1;
}

