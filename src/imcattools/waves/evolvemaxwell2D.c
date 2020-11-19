/* to evolve Maxwell's equations in 2D */

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
#define DT  	0.25

#define swap(a,b) tmp=(a);(b)=(a);(a)=tmp
#define MAX_FRAMES 99999

#define USAGE "\nNAME\n\
	evolvemaxwell2D - evolve a 2-D scalar field\n\
\n\
SYNOPSIS\n\
	evolvemaxwell2D [-nsteps nsteps (8)] [-nframes nframes (99999)] [-autoscale]\n\
		[-dt dt (0.25) [-N N (256)] [-output fieldname (Fxy)] [-u] [-j srcfits]\n\
\n\
DESCRIPTION\n\
	Evolvemaxwell2D evolves\n\
		dot-Ax = - Ftx\n\
		dot-Ay = - Fty\n\
		dot-Ftx = -Fxy,y - jx\n\
		dot-Fty = +Fxy,x - jy\n\
		dot-Fxy = - Ftx,y + Fty,x\n\
\n\
	where comma is spatial derivative computed as:\n\
		G,x[y][x] = 0.5 (G[y][x+1] - G[y][x-1])\n\
	and\n\
		G,x[y][x] = 0.5 (G[y+1][x] - G[y-1][x])\n\
\n\
	The evolution scheme is\n\
		- update Ftx, Fty 1/2 timestep\n\
		- update Fxy (and Ax, Ay) by full time step\n\
		- update Ftx, Fty 1/2 timestep\n\
\n\
	It outputs a 3-D FITS file f[nframes][Ny][Nx] containing the\n\
	evolving field.\n\
\n\
OPTIONS\n\
		-nsteps nsteps		# number of steps between output frames\n\
		-nframes nframes	# total number of output frames\n\
		-dt dt			# time-step\n\
		-N N			# grid-size\n\
		-output opfieldindex	# what to output? options: 0=Ax, 1=Ay, 2=Ftx, 3=Fty, 4=Fxy (4)\n\
		-j srcfits		# read current from srcfits\n\
		-u			# print man-page\n\
\n\
SEE ALSO\n\
	edw, generate_dw, evolvescalar, evolvecomplexscalar\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"



static	fitsheader	*fits;

void    outputframe(float **d, float **dout, int Nx, int Ny, int autoscale, int i);


int	main (int argc, char* argv[]) 
{
	int		arg, nsteps, nframes, N, Nx, Ny, x, y, xp, yp, xm, ym, x0, y0, f, i;
	int		autoscale, opfieldindex, intmode;
	char		*flag, *srcfits;
	double		dt, dx;
	float		**Ax, **Ay, **Ftx, **Fty, **Fxy, **jx, **jy, **fout[5], **dout;
	/* for the current */
	double		j0, sigma, T = 64;
	int		jstep, jstepmax = 128;
	fitsheader	*jfits;
	FILE		*jstream;			

	/* set defaults */
	autoscale = 0;
	nsteps 	= NSTEPS;
	dx 	= 1.0;
	dt 	= DT;
	nframes = MAX_FRAMES;
	N	= 256;
	opfieldindex	= 4;
	sigma	= 8;

	intmode = 1;
	

	/* parse args */
	argsinit(argc, argv, USAGE);
	while ((flag = getflag())) {
		if (!strncmp(flag, "nsteps", 6)) {
			nsteps = getargi();
		} else {
			if (!strncmp(flag, "dt", 2)) {
				dt = getargd();
			} else {
				if (!strncmp(flag, "nframes", 7)) {
					nframes = getargi();
				} else {
					if (!strncmp(flag, "N", 1)) {
						N = getargi();
					} else {
						if (!strncmp(flag, "output", 6)) {
							opfieldindex = getargi();
						} else {
							if (!strncmp(flag, "autoscale", 9)) {
								autoscale = 1;
							} else {
								if (!strncmp(flag, "j", 1)) {
									srcfits = getargs();
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

	Nx = Ny = N;

	/* the current goes here */
	x0 = Nx / 2;
	y0 = Ny / 2;

	if (srcfits) {
		/* we will read the current from an image */
		jstream = fopen(srcfits, "r");
		if (!jstream) {
			error_exit("evolvemaxwell2D: failed to open image for current\n");
		}
		jfits = readfitsheader(jstream);
		if (jfits->ndim != 3) {
			error_exit("evolvemaxwell2D: current fits file must be 3-dimensional\n");
		}
		Ny = jfits->n[1] / 2;
		jfits->n[1] = Ny;
		Nx = jfits->n[0];
		nframes = jfits->n[2];
		fprintf(stderr, "image size = Ny=%d Nx%d  nframes=%d\n", Ny, Nx, nframes);
	} 	
	
	/* create the FITS header */
	fits = new2Dfitsheader(Nx, Ny, FLOAT_PIXTYPE);
	fits->ndim = 3;
	fits->n[2] = 2;
	if (nframes > 1) {
		fits->n[2] = nframes + 1;
	} else {
		fits->ndim = 2;
	}

	/* allocate the memory */
	allocFloatArray(&Ax, Nx, Ny);
	allocFloatArray(&Ay, Nx, Ny);
	allocFloatArray(&Ftx, Nx, Ny);
	allocFloatArray(&Fty, Nx, Ny);
	allocFloatArray(&Fxy, Nx, Ny);
	fout[0] = Ax;
	fout[1] = Ay;
	fout[2] = Ftx;
	fout[3] = Fty;
	fout[4] = Fxy;
	if (autoscale) {
		allocFloatArray(&dout, Nx, Ny);
	} else {
		dout = fout[opfieldindex];
	}

	/* output header */
	add_comment(argc, argv, fits);
	writefitsheader(fits);


        /* output initial frame
        if (nframes > 1) {
                outputframe(fout[opfieldindex], dout, Nx, Ny, autoscale, 0);
        } */


	/* create the jx[][] and jy[][] arrays */
	allocFloatArray(&jx, Nx, Ny);
	allocFloatArray(&jy, Nx, Ny);
	for (y = 0; y < Ny; y++) {
		for (x = 0; x < Nx; x++) {
			jx[y][x] = exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / (2 * sigma * sigma));
		}
	}

	/* evolve */
	for (f = 0; f < nframes; f++) {
		if (srcfits) {
			/* we read the current from the file */
			j0 = 1.0;
			readfitsplane((void **) jx, jfits);
			readfitsplane((void **) jy, jfits);
		}
		for (i = 0; i < nsteps; i++) {
			jstep = nsteps * f + i;
			/* oscillating current at the origin */
			j0 = cos(2 * M_PI * jstep * dt / T);
			/* j0 = 1.0; */
			for (y = 0; y < Ny; y++) {
				ym = (Ny + y - 1) % Ny;
				yp = (y + 1) % Ny;
				for (x = 0; x < Nx; x++) {
					xm = (Nx + x - 1) % Nx;
					xp = (x + 1) % Nx;
					/* first half time step for Ftx, Fty */
					Ftx[y][x] -= 0.5 * dt * j0 * jx[y][x];
					Fty[y][x] -= 0.5 * dt * j0 * jy[y][x];
					if (intmode) {
						Ftx[y][x] += 0.25 * dt * (Fxy[ym][x] - Fxy[yp][x]);
						Fty[y][x] += 0.25 * dt * (Fxy[y][xp] - Fxy[y][xm]);
					} else {
						Ftx[y][x] -= 0.500 * dt * (Ax[yp][x] - 2 * Ax[y][x] + Ax[ym][x]);			/* -0.5 dt Ax,yy */
						Ftx[y][x] += 0.125 * dt * (Ay[yp][xp] - Ay[yp][xm] - Ay[ym][xp] + Ay[ym][xm]);		/* +0.5 dt Ay,xy */
						Fty[y][x] -= 0.500 * dt * (Ay[y][xp] - 2 * Ay[y][x] + Ay[y][xm]);			/* -0.5 dt Ay,xx */
						Fty[y][x] += 0.125 * dt * (Ax[yp][xp] - Ax[yp][xm] - Ax[ym][xp] + Ax[ym][xm]);		/* +0.5 dt Ax,yx */
					}
				}
			}
			for (y = 0; y < Ny; y++) {
				ym = (Ny + y - 1) % Ny;
				yp = (y + 1) % Ny;
				for (x = 0; x < Nx; x++) {
					xm = (Nx + x - 1) % Nx;
					xp = (x + 1) % Nx;
					/* full time step for A and B */
					Ax[y][x] -= dt * Ftx[y][x];	
					Ay[y][x] -= dt * Fty[y][x];
					Fxy[y][x] += 0.5 * dt * (Fty[y][xp] - Fty[y][xm] - Ftx[yp][x] + Ftx[ym][x]);			/* Bdot */
				}
			}
			for (y = 0; y < Ny; y++) {
				ym = (Ny + y - 1) % Ny;
				yp = (y + 1) % Ny;
				for (x = 0; x < Nx; x++) {
					xm = (Nx + x - 1) % Nx;
					xp = (x + 1) % Nx;
					/* second half time step for Ftx, Fty */
					Ftx[y][x] -= 0.5 * dt * j0 * jx[y][x];
					Fty[y][x] -= 0.5 * dt * j0 * jy[y][x];
					if (intmode) {
						Ftx[y][x] += 0.25 * dt * (Fxy[ym][x] - Fxy[yp][x]);
						Fty[y][x] += 0.25 * dt * (Fxy[y][xp] - Fxy[y][xm]);
					} else {
						Ftx[y][x] -= 0.500 * dt * (Ax[yp][x] - 2 * Ax[y][x] + Ax[ym][x]);			/* -0.5 dt Ax,yy */
						Ftx[y][x] += 0.125 * dt * (Ay[yp][xp] - Ay[yp][xm] - Ay[ym][xp] + Ay[ym][xm]);		/* +0.5 dt Ay,xy */
						Fty[y][x] -= 0.500 * dt * (Ay[y][xp] - 2 * Ay[y][x] + Ay[y][xm]);			/* -0.5 dt Ay,xx */
						Fty[y][x] += 0.125 * dt * (Ax[yp][xp] - Ax[yp][xm] - Ax[ym][xp] + Ax[ym][xm]);		/* +0.5 dt Ax,yx */
					}
				}
			}
		}
		outputframe(fout[opfieldindex], dout, Nx, Ny, autoscale, f + 1);
	}

	writefitstail(fits);

	exit(0);
}

void    outputframe(float **d, float **dout, int Nx, int Ny, int autoscale, int i)
{
        int     x, y;
        double  dmin, dmax, ddsum, dsum, dbar, sigmad;

        fprintf(stderr, "# i = %d ", i);
        if (autoscale) {
                dmin = dmax = d[0][0];
                ddsum = dsum = 0.0;
                /* get the scaling */
                for (y = 0; y < Ny; y++) {
                        for (x = 0; x < Nx; x++) {
                                dmin    = (d[y][x] < dmin ? d[y][x] : dmin);
                                dmax    = (d[y][x] > dmax ? d[y][x] : dmax);
                                dsum    += d[y][x];
                                ddsum   += d[y][x] * d[y][x];
                        }
                }
                dbar    = dsum / (Nx * Ny);
                sigmad  = sqrt(ddsum / (Nx * Ny) - dbar * dbar);
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


