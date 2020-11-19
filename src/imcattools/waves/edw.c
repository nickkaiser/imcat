#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "imlib/fits.h"
#include "utils/error.h"
#include "utils/args.h"
#include "utils/arrays.h"
#include "fftlib/myfft.h"

#define usage "\nNAME\n\
	edw - evolve dispersive waves\n\
\n\
SYNOPSIS\n\
	edw nframes dt [-ocean | -debroglie | -scalar kstar] [-verbose] [-novideo]\n\
\n\
DESCRIPTION\n\
	Edw reads a 3D FITS image f[2][Ny][Nx] from standard input,\n\
	the two planes of which are the initial field f = f[0] and\n\
	the initial velocity fdot = f[1].\n\
\n\
	It computes the transforms fk of f and fdotk of fdot and forms\n\
	the positive frequency component:\n\
		fp(k) 	= (fk + fdotk / i omega)\n\
	It then computes and sends to stdout a stream of nframes images\n\
	containing the real part of the inverse transform of\n\
		fp(k) * exp(i omega(k) * t).\n\
	for t = integer multiples of dt.\n\
\n\
	By default, the dispersion relation is omega = k.\n\
	With option -ocean we use omega = sqrt(k).\n\
\n\
	With the -novideo option we do not output a stream\n\
	of images, rather we evolve for one step\n\
	and then output an image f[2][Ny][Nx] containing the\n\
	final field and velocity.\n\
\n\
	With -verbose option we tell stderr what we are doing.\n\
\n\
SEE ALSO\n\
	generate_dw, evolvescalar, xfv\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"

float	omega(int ikx, int iky);
float	t, dkx, dky, kstar;
int	dr;

#define DR_DEFAULT	0
#define DR_OCEAN	1
#define DR_DEBROGLIE	2
#define DR_SCALAR	3

main (int argc, char *argv[])
{
	double		dt;
	float		**f, **fdot, **fk, **fdotk, **fdotkr, **fdotki, **fkr, **fki, ftmp, omeg, cosw, sinw;
	int		i, x, y, kx, ky, Nx, Ny, nsteps, verbose, novideo;
	fitsheader	*fits;
	char		*flag;

	/* defaults */
	nsteps 	= 99999;
	dr 	= DR_DEFAULT;	
	verbose = 0;
	novideo = 0;

	/* parse args */
	argsinit(argc, argv, usage);
	if (nextargtype() == FLAG_ARG) {
		error_exit(usage);
	}
	nsteps = getargd();
	dt = getargd();
	while(flag = getflag()) {
		if (!strncmp(flag, "ocean", 5)) {
			dr = DR_OCEAN;
		} else {
			if (!strncmp(flag, "debroglie", 9)) {
				dr = DR_DEBROGLIE;
			} else {
				if (!strncmp(flag, "scalar", 6)) {
					dr = DR_SCALAR;
					kstar = getargf();
				} else { 	
					if (!strncmp(flag, "verbose", 7)) {
						verbose = 1;
					} else {
						if (!strncmp(flag, "novideo", 7)) {
							novideo = 1;
							nsteps = 1;
						} else {
							error_exit(usage);
						}
					}
				}
			}
		}
	}

	if (verbose) {
		switch (dr) {
			case DR_DEFAULT :
				fprintf(stderr, "# edw : non-dispersive waves\n");
				break;
			 case DR_OCEAN :
				fprintf(stderr, "# edw : deep ocean waves\n");
				break;
			 case DR_DEBROGLIE :
				fprintf(stderr, "# edw : de Broglie waves\n");
				break;
			 case DR_SCALAR :
				fprintf(stderr, "# edw : scalar field\n");
				break;
		}
	} 

	/* read the fits header */
	fits = readfitsheader(stdin);
	if (fits->ndim != 3) {
		error_exit("edw: input image must be three dimensional\n");
	}
	if (fits->n[2] != 2) {
		error_exit("edw: NAXIS3 must be 2\n");
	}
	if (!novideo) {
		fits->n[2] = nsteps;
	}
	Nx = fits->n[0];
	Ny = fits->n[1];
	add_comment(argc, argv, fits);
	writefitsheader(fits);

	/* set globals dkx, dky */
	dkx = 2 * M_PI / Nx;
	dky = 2 * M_PI / Ny;

	/* allocate and read the input data */
	allocFloatArray(&f, Nx, 2 * Ny);
	allocFloatArray(&fdot, Nx, 2 * Ny);
	readfitsplane((void *) f, fits);
	readfitsplane((void *) fdot, fits);

	/* compute the transforms */
	allocFloatArray(&fk, Nx, 2 * Ny);
	allocFloatArray(&fdotk, Nx, 2 * Ny);
	forward_cfft(f, Nx, Ny, fk);
	forward_cfft(fdot, Nx, Ny, fdotk);

	/* get pointers to real and imaginary parts of fk and fdotk */
	fdotkr = fdotk;
	fdotki = fdotk + Ny;
	fkr = fk;
	fki = fk + Ny;

	/* multiply fdotk by 1 / i omega */
	for (ky = 0; ky < Ny; ky++) {
		for (kx = 0; kx < Nx; kx++) {
			omeg = omega(kx - Nx / 2, ky - Ny / 2);
			ftmp = fdotkr[ky][kx];
			fdotkr[ky][kx] = fdotki[ky][kx] / omeg;
			fdotki[ky][kx] = - ftmp / omeg; 
		}
	}

	/* add fk and fdotk / i omega */
	for (ky = 0; ky < 2 * Ny; ky++) {
		for (kx = 0; kx < Nx; kx++) {
			fk[ky][kx] += fdotk[ky][kx];
		}
	}

	/* evolve */
	for (i = 0; i < nsteps; i++) {
		if (novideo) {
			t = dt;
		} else {
			t = i * dt;
		}
		/* apply evolution phase factor (store result in fdotk) */
		for (ky = 0; ky < Ny; ky++) {
			for (kx = 0; kx < Nx; kx++) {
				omeg = omega(kx - Nx / 2, ky - Ny / 2); 
				cosw = cos(omeg * t);
				sinw = sin(omeg * t);
				fdotkr[ky][kx] = cosw * fkr[ky][kx] + sinw * fki[ky][kx];
				fdotki[ky][kx] = cosw * fki[ky][kx] - sinw * fkr[ky][kx];
			}
		}		
		inverse_cfft(fdotk, Nx, Ny, f);
		writefitsplane((void *) f, fits);
		if (novideo) {
			/* apply evolution factor and multiply by i omega */
			for (ky = 0; ky < Ny; ky++) {
				for (kx = 0; kx < Nx; kx++) {
					omeg = omega(kx - Nx / 2, ky - Ny / 2);
					cosw = cos(omeg * t);
					sinw = sin(omeg * t);
					fdotkr[ky][kx] = omeg * (sinw * fkr[ky][kx] - cosw * fki[ky][kx]);
					fdotki[ky][kx] = omeg * (cosw * fkr[ky][kx] + sinw * fki[ky][kx]);
				}
			}
			inverse_cfft(fdotk, Nx, Ny, fdot);
			writefitsplane((void *) fdot, fits);
		}
	}	


	exit(0);
}


float	omega(int ikx, int iky)
{
	float	kx, ky, kk0;

	kx = ikx * dkx;
	ky = iky * dky;
	kk0 = 0.01 * dkx * dky;

	switch (dr) {
		case DR_DEBROGLIE :
			return (kk0 + kx * kx + ky * ky);
			break;
		case DR_OCEAN :
			return (pow(kk0 + kx * kx + ky * ky, 0.25));
			break;
		case DR_SCALAR :
			return (sqrt(kstar * kstar + kx * kx + ky * ky));
			break;
		default:
			return (sqrt(kk0 + kx * kx + ky * ky));
			break;
	}
}
