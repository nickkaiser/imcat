/*
 * translate.c
 */

#define usage "\n\
NAME\n\
	translate_fft --- translate an image using 'sinc-interpolation'\n\
\n\
SYNOPSIS\n\
	translate_fft dy dx\n\
\n\
DESCRIPTION\n\
	'translate_fft' uses the FFT to translate an image.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n"

#include <stdio.h>
#include <math.h>

#include "../imlib/fits.h"
#include "../utils/arrays.h"
#include "../fftlib/myfft.h"

float	rfunc(float ki, float kj);
float	ifunc(float ki, float kj);
float	di, dj;

main (int argc, char *argv[])
{
	int		arg = 1, N1, N2, mode, i, j, pixtype;
	float		**f;
	fitsheader	*fits;
	fft_type 	fk;

	if (argc < 3) {
		fprintf(stderr, usage);
		exit(-1);
	}
	sscanf(argv[1], "%f", &dj);
	sscanf(argv[2], "%f", &di);

	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);
	alloc_fft(&fk, N1, N2);
	forward_fft(f, N1, N2, fk);
	cfilter(fk, N1, N2, rfunc, ifunc);
	inverse_fft(fk, N1, N2, f);
	
	add_comment(argc, argv, fits);
	write2Dfloatimage(f, fits);
	exit(0);
}


float	rfunc(float ki, float kj)
{
	return (cos(ki * di + kj * dj));
}


float	ifunc(float ki, float kj)
{
	return (sin(ki * di + kj * dj));
}




#undef EPS



