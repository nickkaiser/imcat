/*
 * phasetopsf.c --- generate realisations of psf
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils/error.h"
#include "utils/arrays.h"
#include "utils/args.h"
#include "imlib/fits.h"
#include "fftlib/myfft.h"


#define usage "\n\
NAME\n\
	phasetopsf - computes PSF from phase screen FITS image\n\
\n\
SYNOPSIS\n\
	phasetopsf [-N Nout] [-D D] [-p pupil]\n\
\n\
	Options are:\n\
		-N Nout		# output only central Nout x Nout subimage (N)\n\
		-D D		# mirror diameter in pixels (N)\n\
		-p pupilim	# supply pupil image.\n\
\n\
DESCRIPTION\n\
\n\
	phasetopsf reads a stream of N x N phase screen samples from stdin and sends to\n\
	stdout a stream of psfs.\n\
\n\
	With -p option we ignore any D value supplied with the -D option and we\n\
	read the pupil from 'pupilim'.  If this is 2N x 2N then it is treated as a\n\
	real amplitude transmission function $0 <= T <= 1$.  If it is 2 x 2N x 2N\n\
	the first plane is consisered to be the real amplitude transmission function\n\
	and the second plane is the phase in radians.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n\n"

main (int argc, char *argv[])
{
	int		N, M, Nout, xyoff, frame, nframes, ND;
	float		**phi, **T, **Tphi, **Tc, **Ts, **Tcr, **Tci, **Tsr, **Tsi, **g;
	int		ix, iy, ixx, iyy, usepupil, usepupilphase;
	double		a, b, r, x, y, Tsum;
	fft_type	Tck, Tsk;
	fitsheader	*fitsi, *fitso, *fitsp;
	char		*flag, *pupilfilename;
	FILE		*pupilf;

	/* defaults */
	Nout = 0;
	ND = 0;
	usepupil = 0;
	usepupilphase = 0;

	/* parse args */
	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'N':
				Nout = getargi();
				break;
			case 'D':
				ND = getargi();
				break;
			case 'p':
				usepupil = 1;
				pupilfilename = getargs();
				break;
			default:
				error_exit(usage);
		}
	}
	
	/* read the FITS image header */
	fitsi = readfitsheader(stdin);
	if (fitsi->ndim != 3) {
		error_exit("phasetopsf: input image must be 3D fits file\n");
	}
	if (fitsi->n[0] != fitsi->n[1]) {
		error_exit("phasetopsf: input image planes must be square\n");
	}
	N = fitsi->n[0];
	nframes = fitsi->n[2];

	/* set output frame size and pupil diameter unless supplied */
	Nout = (Nout ? Nout : N);
	ND = (ND ? ND : N);

	/* we compute in boxes this size */
	M = 2 * N;

	/* offset for output sub image */
	xyoff = (M - Nout) / 2;

	/* create the output header */
	fitso = copyfitsheader(fitsi);
	fitso->n[0] = fitso->n[1] = Nout;
	add_comment(argc, argv, fitso);
	writefitsheader(fitso);

	if (usepupil) {
		/* read the pupil file header */
		pupilf = fopen(pupilfilename, "r");
		if (!pupilf) {
			error_exit("phasetopsf : failed to open specified pupil file\n");
		}
		fitsp = readfitsheader(pupilf);
		if ((fitsp->n[0] != M) || (fitsp->n[1] != M)) {
			error_exit("phasetopsf : pupil image plane dimensions must be (2N) x (2N)\n");
		}
		switch (fitsp->ndim) {
			case 2:
				usepupilphase = 0;
				break;
			case 3:
				usepupilphase = 1;
				break;
			default:
				error_exit("phasetopsf : pupil image must be 2 or three dimensional\n");
		}
	}

	/* allocation of arrays */
	/* allocate the phase screen array */
	allocFloatArray(&phi, N, N);
	/* allocate the pupil function T, Tc, Ts, g etc */
	allocFloatArray(&T, M, M);
	if (usepupilphase) {
		allocFloatArray(&Tphi, M, M);
	}
	allocFloatArray(&Tc, M, M);
	allocFloatArray(&Ts, M, M);
	alloc_fft(&Tck, M, M);
	alloc_fft(&Tsk, M, M);
	allocFloatArray(&Tcr, M, 2 * M);
	Tci = Tcr + M;
	allocFloatArray(&Tsr, M, 2 * M);
	Tsi = Tsr + M;
	allocFloatArray(&g, Nout, Nout);

	/* make or read the pupil function T */
	if (usepupil) {
		readfitsplane((void *) T, fitsp);
		if (usepupilphase) {
			readfitsplane((void *) Tphi, fitsp);
		}
		fclose(pupilf);
	} else {
		for (iy = 0; iy < M; iy++) {
			y = iy - 0.5 * M;
			for (ix = 0; ix <  M; ix++) {
				x = ix - 0.5 * M;
				r = sqrt(x * x + y * y);
				T[iy][ix] = (r < 0.5 * ND ? 1.0 : 0.0);
				Tsum += T[iy][ix];
			}
		}
	}

	/* normalise it */
	Tsum = 0;
	for (iy = 0; iy < M; iy++) {
		for (ix = 0; ix <  M; ix++) {
			Tsum += T[iy][ix];
		}
	}
	Tsum = sqrt(Tsum);
	for (iy = 0; iy < M; iy++) {
		for (ix = 0; ix <  M; ix++) {
			T[iy][ix] /= Tsum;
		}
	}

	for (frame = 0; frame < nframes; frame++) {
		/* read the phase screen */
		readfitsplane((void *) phi, fitsi);
		/* compute T cos(phi) and T sin(phi) */
		for (iy = 0; iy < N; iy++) {
			iyy = iy + N / 2;
			for (ix = 0; ix <  N; ix++) {
				ixx = ix + N / 2;
				if (usepupilphase) {
					Tc[iyy][ixx] = T[iyy][ixx] * cos(phi[iy][ix] + Tphi[iyy][ixx]);
					Ts[iyy][ixx] = T[iyy][ixx] * sin(phi[iy][ix] + Tphi[iyy][ixx]);
				} else {
					Tc[iyy][ixx] = T[iyy][ixx] * cos(phi[iy][ix]);
					Ts[iyy][ixx] = T[iyy][ixx] * sin(phi[iy][ix]);
				}
			}
		}
		forward_fft(Tc, M, M, Tck);
		forward_fft(Ts, M, M, Tsk);
		get_fft(Tck, M, M, Tcr);
		get_fft(Tsk, M, M, Tsr);
		for (iy = 0; iy < Nout; iy++) {
			for (ix = 0; ix < Nout; ix++) {
				a = Tcr[xyoff + iy][xyoff + ix] - Tsi[xyoff + iy][xyoff + ix];
				b = Tsr[xyoff + iy][xyoff + ix] + Tci[xyoff + iy][xyoff + ix];
				g[Nout - 1 - iy][Nout - 1 - ix] = (a * a + b * b) / (M * M);
	    		}
		}
		writefitsplane((void *) g, fitso);
	}
	exit(0);
}

