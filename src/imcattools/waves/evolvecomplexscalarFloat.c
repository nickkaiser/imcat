/* to evolve an interacting scalar field */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "utils/args.h"
#include "utils/arrays.h"
#include "utils/error.h"
#include "imlib/fits.h"

/* defaults */
#define NSTEPS 	8
#define KSTAR	0.3
#define DT  	0.25
#define NOP	6

#define swap(a,b) tmp=(a);(b)=(a);(a)=tmp
#define MAX_FRAMES 99999

#define USAGE "\nNAME\n\
	evolvecomplexscalar - evolve complex scalar field in 2-dimensions\n\
\n\
SYNOPSIS\n\
	evolvecomplexscalar [-nsteps nsteps (8)] [-nframes nframes (99999)] [-q q (0.1)]\n\
		[-m m (0.3)] [-autoscale] [-dt dt (0.25)] [-op opbase (tmp/ECS)] [-u]\n\
\n\
DESCRIPTION\n\
	Evolvecomplexscalar reads from stdin a 4-D FITS file f[2][2][Ny][Nx] where\n\
		f[0][0] and f[0][1] are the real and imaginary parts of a scalar field d\n\
	and where\n\
		f[1][0] and f[1][1] are the real and imaginary parts of its velocity v.\n\
\n\
	It then evolves the coupled equations\n\
\n\
		dv/dt = laplacian(d) - m^2 d\n\
	and\n\
		dd/dt = v\n\
\n\
	which are the Klein-Gordon equation for a complex scalar d with mass = m,\n\
	in the the Hamitonian gauge: A = (0,Ax,Ay) (so D/dt = d/dt)\n\
\n\
	It also evolves the Ampere/Maxwell equation\n\
		dFtx/dt = -dFxy/dy - jx\n\
		dFty/dt = +dFxy/dx - jy\n\
	and Faraday's law\n\
		dFxy/dt = dFty/dx - dFyx/dy\n\
	with vector potential determined from\n\
		dAx/dt = -Ftx\n\
		dAy/dt = -Fty\n\
	and the current determined from\n\
		jx = -0.5 q Im[d * conj((d/dx - i q Ax)d) - conj(d) * (d/dx - i q Ax)d]\n\
	and similarly for jy.\n\
\n\
	The (covariant) Laplacian function is\n\
\n\
		laplacian = div.grad d - 2 i q (A.grad) d - i q d div.A - q^2 A.A d\n\
\n\
	with the ordinary Laplacian computed as\n\
\n\
		div.grad d = d[y][x-1] + d[y][x+1] + d[y-1][x] + d[y+1][x] - 4 * d[y][x].\n\
\n\
	The evolution scheme is:\n\
		- advance d, Fxy, Ax and Ay by dt and calculate jx, jy\n\
		- advance v, Ftx and Fty by dt\n\
\n\
	It outputs 6 3-D FITS files f[nframes][][] containing\n\
		d:	real and imaginary parts\n\
		j:	jx and jy\n\
		A:	Ax and Ay\n\
		E:	Ftx and Fty\n\
		B:	Fxy\n\
		r:	charge densities from (d v^* - d^* v) and from Fti,i\n\
	which by default are written to ./tmp/ECS.[djAEBr].fits\n\
	but you can modify the base of the file names with the -op option.\n\
\n\
	For d, j, A and E, the 2 components are concatenated in a Nx by 2 * Ny image.\n\
\n\
OPTIONS\n\
		-nsteps nsteps		# number of steps between output frames\n\
		-nframes nframes	# total number of output frames\n\
		-m m			# mass (Compton frequency)\n\
		-q q			# the charge parameter\n\
		-autoscale		# scale each output frame to facilitate visualisation\n\
		-dt dt			# time-step\n\
		-u			# print man-page\n\
\n\
SEE ALSO\n\
	edw, generate_dw, evolvemaxwell2D, evolvescalar\n\
\n\
AUTHOR\n\
	Nick Kaiser\n\n"



static	fitsheader	*fits;

void    allocComplexFloatArray(float complex ***f, int N1, int N2);

int	main (int argc, char* argv[]) 
{
	int		arg, nsteps, nframes, Nx, Ny, x, y, xp, yp, xm, ym, iframe, istep, i, autoscale;
	char		*flag, opbase[128], opfile[128];
	double		dt, dx, m, msquared, q, epsilon;
	float		**dreal, **dimag, **vreal, **vimag;
	float complex	**d, **v, laplacian, Ddx, Ddy;
	float		**Ax, **Ay, **Ftx, **Fty, **Fxy, **jx, **jy;
	FILE		*opstream[NOP];
	fitsheader	*opfits[NOP];
	char		opchar[7][2] = {"d", "j", "A", "E", "B", "r"};

	/* set defaults */
	autoscale = 0;
	nsteps 	= NSTEPS;
	m 	= KSTAR;
	dx 	= 1.0;
	dt 	= DT;
	nframes = MAX_FRAMES;
	q	= 0.1;
	strcpy(opbase, "tmp/ECS");

	/* parse args */
	argsinit(argc, argv, USAGE);
	while ((flag = getflag())) {
		if (!strncmp(flag, "nsteps", 6)) {
			nsteps = getargi();
		} else {
			if (!strncmp(flag, "m", 1)) {
				m = getargd();
			} else {
				if (!strncmp(flag, "autoscale", 9)) {
					autoscale = 1;
				} else {
					if (!strncmp(flag, "dt", 2)) {
						dt = getargd();
					} else {
						if (!strncmp(flag, "nframes", 7)) {
							nframes = getargi();
						} else {
							if (!strncmp(flag, "op", 2)) {
								strcpy(opbase, getargs());
							} else {
								if (!strncmp(flag, "q", 1)) {
									q = getargd();
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

	/* compute mass-squared coefficient */
	msquared = m * m;
	
	/* read the data header */
	fits = readfitsheader(stdin);
	if (fits->ndim != 4 || fits->n[3] != 2 || fits->n[2] != 2) {
		error_exit("evolvescalar: source data must be a f[2][2][Ny][Nx] FITS file\n");
	}
	Nx = fits->n[0];
	Ny = fits->n[1];

	/* allocate the memory for the input arrays*/
	allocFloatArray(&dreal, Nx, Ny);
	allocFloatArray(&vreal, Nx, Ny);
	allocFloatArray(&dimag, Nx, Ny);
	allocFloatArray(&vimag, Nx, Ny);

	/* and for the complex arrays */
	allocComplexFloatArray(&d, Nx, Ny);
	allocComplexFloatArray(&v, Nx, Ny);

	/* and for the EM field */
	allocFloatArray(&Ax, Nx, Ny);
	allocFloatArray(&Ay, Nx, Ny);
	allocFloatArray(&Ftx, Nx, Ny);
	allocFloatArray(&Fty, Nx, Ny);
	allocFloatArray(&Fxy, Nx, Ny);

	/* and for the current */
	allocFloatArray(&jx, Nx, Ny);
	allocFloatArray(&jy, Nx, Ny);

	/* read the initial data and load into the complex arrays */
	readfitsplane((void **) dreal, fits);
	readfitsplane((void **) vreal, fits);
	readfitsplane((void **) dimag, fits);
	readfitsplane((void **) vimag, fits);
	for (y = 0; y < Ny; y++) {
		for (x = 0; x < Nx; x++) {
			d[y][x] = dreal[y][x] + I * dimag[y][x];
			v[y][x] = vreal[y][x] + I * vimag[y][x];
		}
	}

	/* open output streams */
	add_comment(argc, argv, fits);
	fits->ndim = 3;
	fits->n[0] = Nx;
	fits->n[1] = 2 * Ny;		/* for all except B */
	fits->n[2] = nframes;
	for (i = 0; i < NOP; i++) {
		sprintf(opfile, "%s.%s.fits", opbase, opchar[i]);
		fprintf(stderr, "opfile = %s\n", opfile);
		opstream[i] = fopen(opfile, "w");
		if (!opstream[i]) { 
			error_exit("evolvecomplexscalar: failed to open opstream\n");
		}
		opfits[i] = copyfitsheader(fits);
		(opfits[i])->opstream = opstream[i];
		if (i == 4) {
			(opfits[4])->n[1] = Ny;
		}
		writefitsheader(opfits[i]);
		(opfits[i])->n[1] = Ny;
	}

	/* evolve */
	for (iframe = 0; iframe < nframes; iframe++) {
		fprintf(stderr, "# frame %d\n", iframe);
		for (istep = 0; istep < nsteps; istep++) {
			/* advance d, Fxy, Ax and Ay by dt and calculate jx, jy */
			for (y = 0; y < Ny; y++) {
				ym = (Ny + y - 1) % Ny;
				yp = (y + 1) % Ny;
				for (x = 0; x < Nx; x++) {
					xm = (Nx + x - 1) % Nx;
					xp = (x + 1) % Nx;
					d[y][x] += dt * v[y][x];
					Ax[y][x] -= dt * Ftx[y][x];	
					Ay[y][x] -= dt * Fty[y][x];
					Fxy[y][x] += 0.5 * dt * (Fty[y][xp] - Fty[y][xm] - Ftx[yp][x] + Ftx[ym][x]);
					Ddx = 0.5 * (d[y][xp] - d[y][xm]) - q * I * Ax[y][x] * d[y][x];
					Ddy = 0.5 * (d[yp][x] - d[ym][x]) - q * I * Ay[y][x] * d[y][x];
					jx[y][x] = -0.5 * q * cimagf(d[y][x] * conjf(Ddx) - conjf(d[y][x]) * Ddx);
					jy[y][x] = -0.5 * q * cimagf(d[y][x] * conjf(Ddy) - conjf(d[y][x]) * Ddy);
				}
			}								
			/* advance v, Ftx, Ftx */
			for (y = 0; y < Ny; y++) {
				ym = (Ny + y - 1) % Ny;
				yp = (y + 1) % Ny;
				for (x = 0; x < Nx; x++) {
					xm = (Nx + x - 1) % Nx;
					xp = (x + 1) % Nx;
					laplacian = d[y][xm] + d[y][xp] + d[ym][x] + d[yp][x] - 4 * d[y][x]
						- q * I * (Ax[y][x] * (d[y][xp] - d[y][xm]) + Ay[y][x] * (d[yp][x] - d[ym][x])
						+ 0.5 * d[y][x] * (Ax[y][xp] - Ax[y][xm] + Ay[yp][x] - Ay[ym][x]))
						- q * q * (Ax[y][x] * Ax[y][x] + Ay[y][x] * Ay[y][x]) * d[y][x];
					v[y][x] += dt * (laplacian - msquared * d[y][x]);
					Ftx[y][x] += 0.5 * dt * (Fxy[ym][x] - Fxy[yp][x] - 2 * jx[y][x]);
					Fty[y][x] += 0.5 * dt * (Fxy[y][xp] - Fxy[y][xm] - 2 * jy[y][x]);
				}
			}
		}
		/* output */
		/* the scalar field */
		for (y = 0; y < Ny; y++) {
			for (x = 0; x < Nx; x++) {
				dreal[y][x] = crealf(d[y][x]);
				dimag[y][x] = cimagf(d[y][x]);
			}
		}
		writefitsplane((void **) dreal, opfits[0]);
		writefitsplane((void **) dimag, opfits[0]);
		/* the current */
		writefitsplane((void **) jx, opfits[1]);
		writefitsplane((void **) jy, opfits[1]);
		/* the potential A */
		writefitsplane((void **) Ax, opfits[2]);
		writefitsplane((void **) Ay, opfits[2]);
		/* the E-field */
		writefitsplane((void **) Ftx, opfits[3]);
		writefitsplane((void **) Fty, opfits[3]);
		/* the B-field */
		writefitsplane((void **) Fxy, opfits[4]);
		/* the charge density - from Fti,i and from (d v^* - d^* v) */
		for (y = 0; y < Ny; y++) {
			ym = (Ny + y - 1) % Ny;
			yp = (y + 1) % Ny;
			for (x = 0; x < Nx; x++) {
				xm = (Nx + x - 1) % Nx;
				xp = (x + 1) % Nx;
				dreal[y][x] = 0.5 * q * cimagf(d[y][x] * conjf(v[y][x]) - conjf(d[y][x]) * v[y][x]);
				dimag[y][x] = 0.5 * (Ftx[y][xp] - Ftx[y][xm] + Fty[ym][x] - Fty[ym][x]);
			}
		}
		writefitsplane((void **) dreal, opfits[5]);
		writefitsplane((void **) dimag, opfits[5]);
	}
					
	for (i = 0; i < NOP; i++) {
		if (i != 4) {
			(opfits[i])->n[1] = 2 * Ny;
		}
		writefitstail(opfits[i]);
		fclose(opstream[i]);
	}

	exit(0);
}


void    allocComplexFloatArray(float complex ***f, int N1, int N2)
{
        int             i;
        
        (*f) = (float complex **) calloc(N2, sizeof(float complex *));
        if (!*f)
                error_exit("allocComplexFloatArray: memory allocation failure\n");
        (*f)[0] = (float complex *) calloc(N1 * N2, sizeof(float complex));
        if (!(*f)[0])
                error_exit("allocComplexFloatArray: memory allocation failure\n");
        for (i = 1; i < N2; i++)
                (*f)[i] = (*f)[i - 1] + N1;
}

