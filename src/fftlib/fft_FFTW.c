/*
 *
 * fft_FFTW.c	- FFTW fft routines
 *
 * note that in FFTW docs, the fast index is y, slow index is x
 * the FFT format is complex arrays C[Ny][Nx/2 + 1]
 *
 */


#include <stdio.h>
#include <stdlib.h>
#ifndef Darwin
#include <malloc.h>
#endif
#include <math.h>


#include <fftw.h>
#include <rfftw.h>

#include "myfft.h"

#ifndef PI
#define PI M_PI
#endif

void	myfftwcfft(float **fsrc, int nx, int ny, float **fdst, int direction);

void	forward_fft(float **f, int nx, int ny, fft_type fk)
{
	rfftwnd_plan p;	

	p = rfftw2d_create_plan(ny, nx, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
	rfftwnd_one_real_to_complex(p, &f[0][0], &fk[0][0]);
	rfftwnd_destroy_plan(p);
}

void	inverse_fft(fft_type fk, int nx, int ny, float **f) 
{
	rfftwnd_plan p;	
	float	scale;
	int	x, y;

	p = rfftw2d_create_plan(ny, nx, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
	rfftwnd_one_complex_to_real(p, &fk[0][0], &f[0][0]);
	rfftwnd_destroy_plan(p);
	scale = 1.0 / (nx * ny);
	for (y = 0; y < ny; y++) {
		for (x = 0; x < nx; x++) {
			f[y][x] *= scale;
		}
	}
}


void	filter(fft_type fk, int Nx, int Ny, float (*filterfunc)(float ki, float kj))
{
	int	ix, iy, n;
	float	T, kx, ky;
	fftw_complex **Fk;

	Fk = (fftw_complex **) fk;
	n = Nx / 2 + 1;
	for (iy = 0; iy < Ny; iy++) {
		ky = (iy < Ny / 2 ? iy : iy - Ny) * 2 * PI / (float) Ny;
		for (ix = 0; ix < n; ix++) {
			kx = ix * 2 * PI / (float) Nx;
			T = filterfunc(kx, ky);
			Fk[iy][ix].re *= T;
			Fk[iy][ix].im *= T;
		}
	}
}


void	cfilter(fft_type fk, int Nx, int Ny, float (*rfunc)(float ki, float kj), float (*ifunc)(float ki, float kj))
{
	int	ix, iy, n;
	float	Tr, Ti, fr, fi, kx, ky;
	fftw_complex **Fk;

	Fk = (fftw_complex **) fk;
	n = Nx / 2 + 1;
	for (iy = 0; iy < Ny; iy++) {
		ky = (iy < Ny / 2 ? iy : iy - Ny) * 2 * PI / (float) Ny;
		for (ix = 0; ix < n; ix++) {
			kx = ix * 2 * PI / (float) Nx;
			Tr = rfunc(kx, ky);
			Ti = ifunc(kx, ky);
			/* funny signs below are for compatibility */
			/* with NR convention. Our Fk is really conjugate */
			fr = Fk[iy][ix].re;
			fi = -Fk[iy][ix].im;
			Fk[iy][ix].re = fr * Tr - fi * Ti;
			Fk[iy][ix].im = -1.0 * (fr * Ti + fi * Tr);			
		}
	}
}



void	ccf(fft_type fk1, fft_type fk2, int Nx, int Ny, float **ff12, int x0, int y0)
{
	int	ix, iy, n;
	float	a1, b1, a2, b2;
	fftw_complex **Fk1, **Fk2;

	Fk1 = (fftw_complex **) fk1;
	Fk2 = (fftw_complex **) fk2;
	n = Nx / 2 + 1;
	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < n; ix ++) {
			a1 = Fk1[iy][ix].re;
			a2 = Fk2[iy][ix].re;
			b1 = Fk1[iy][ix].im;
			b2 = Fk2[iy][ix].im;
			Fk1[iy][ix].re = (a1 * a2 + b1 * b2) / (Nx * Ny);
			Fk1[iy][ix].im = (a1 * b2 - b1 * a2) / (Nx * Ny);
		}
	}
	inverse_fft(fk1, Nx, Ny, ff12);
	cycleimage(ff12, Nx, Ny, x0, y0);
}



void	power(fft_type fk, int Nx, int Ny, float **p, int x0, int y0)
{
	int	ix, iy, n, kx, ky;
	fftw_complex **Fk;
	float	scale;

	n = Nx / 2 + 1;
	scale = 1.0 / (Nx * Ny);

	/* fill the right hand half plus nyquist strip */
	for (iy = 0; iy < Ny; iy++) {
		ky = (iy + Ny / 2) % Ny;
		for (ix = 0; ix < n; ix++) {
			kx = (ix + Nx / 2) % Nx;
			p[ky][kx] = scale * (fk[iy][ix].re * fk[iy][ix].re + fk[iy][ix].im * fk[iy][ix].im);
		}
	}
	/* reflect to get the rest */
	for (ky = 0; ky < Ny; ky++) {
		iy = (Ny - ky) % Ny;
		for (kx = 1; kx < (Nx + 1) / 2; kx++) {
			p[ky][Nx / 2 - kx] = p[iy][Nx / 2 + kx];
		}
	}
}


void	alloc_fft(fft_type *fk, int Nx, int Ny)
{
	int 	n, y;

	n = Nx / 2 + 1;
	*fk = (fftw_complex **) calloc(Ny, sizeof(fftw_complex *));
	(*fk)[0] = (fftw_complex *) calloc(Ny * n, sizeof(fftw_complex));
	if (!(*fk) || !((*fk)[0]))
		error_exit("alloc_fft: failed to allocate memory\n");
	for (y = 1; y < Ny; y++) {
		(*fk)[y] = (*fk)[y - 1] + n;	
	}
}



void	free_fft(fft_type fk, int Nx, int Ny)
{
	free(fk[0]);
	free(fk);
}



void	copy_fft(fft_type fksrc, int Nx, int Ny, fft_type fkdest)
{
	int	x, y, n;

	n = Nx / 2 + 1;
	for (y = 0; y < Ny; y++) {
		for (x = 0; x < n; x++) {
			fkdest[y][x] = fksrc[y][x];
		}
	}
}


void	get_fft(fft_type fk, int Nx, int Ny, float **fdest)
{
	int	kx, ky, ix, iy, n;
	float	**fdestreal, **fdestimag;

	fdestreal = fdest;
	fdestimag = fdest + Ny;
	n = Nx / 2 + 1;

	/* fill the right hand half plus nyquist strip */
	for (iy = 0; iy < Ny; iy++) {
		ky = (iy + Ny / 2) % Ny;
		for (ix = 0; ix < n; ix++) {
			kx = (ix + Nx / 2) % Nx;
			fdestreal[ky][kx] = fk[iy][ix].re;
			fdestimag[ky][kx] = -fk[iy][ix].im;
		}
	}
	/* reflect to get the rest */
	for (ky = 0; ky < Ny; ky++) {
		iy = (Ny - ky) % Ny;
		for (kx = 1; kx < (Nx + 1) / 2; kx++) {
			fdestreal[ky][Nx / 2 - kx] = fdestreal[iy][Nx / 2 + kx];
			fdestimag[ky][Nx / 2 - kx] = -fdestimag[iy][Nx / 2 + kx];
		}
	}
}


void	set_fft(fft_type fk, int Nx, int Ny, float **fsrc)
{
	int	kx, ky, ix, iy, iyy, n;
	float	**freal, **fimag;

	freal = fsrc;
	fimag = fsrc + Ny;
	n = Nx / 2 + 1;

	/* get the right hand half plus nyquist strip */
	for (iy = 0; iy < Ny; iy++) {
		ky = (iy + Ny / 2) % Ny;
		for (ix = 0; ix < n; ix++) {
			kx = (ix + Nx / 2) % Nx;
			fk[iy][ix].re = freal[ky][kx];
			fk[iy][ix].im = -fimag[ky][kx];
		}
	}
}


void	forward_cfft(float **fsrc, int nx, int ny, float **fdst)
{
	myfftwcfft(fsrc, nx, ny, fdst, FFTW_FORWARD);
}

void	inverse_cfft(float **fsrc, int nx, int ny, float **fdst) 
{
	myfftwcfft(fsrc, nx, ny, fdst, FFTW_BACKWARD);
}

void	myfftwcfft(float **fsrc, int nx, int ny, float **fdst, int direction)
{
	fftw_complex 	*fk;
	fftwnd_plan 	p;
	float		**fr, **fi, scale;
	int		x, y;

	/* create the plan */
	p = fftw2d_create_plan(ny, nx, direction, FFTW_ESTIMATE | FFTW_IN_PLACE);

	/* allocate the complex array */
	fk = (fftw_complex *) calloc(nx * ny, sizeof(fftw_complex));

	/* pointers to sub-images */
	fr = fsrc;
	fi = fsrc + ny;

	if (direction == FFTW_BACKWARD) {
		/* shift the origin */
        	cycleimage(fr, nx, ny, -nx / 2, -ny / 2);
        	cycleimage(fi, nx, ny, -nx / 2, -ny / 2);
		scale = 1.0 / (nx * ny);
	} else {
		scale = 1.0;
	}

	/* fill the complex array */
	for (y = 0; y < ny; y++) {
		for (x = 0; x < nx; x++) {
			fk[nx * y + x].re = scale * fr[y][x];
			fk[nx * y + x].im = scale * fi[y][x];
		}
	}

	/* do the transform */
     	fftwnd_one(p, fk, NULL);

	/* extract the results */
	fr = fdst;
	fi = fdst + ny;
	for (y = 0; y < ny; y++) {
		for (x = 0; x < nx; x++) {
			fr[y][x] = fk[nx * y + x].re;
			fi[y][x] = fk[nx * y + x].im;
		}
	}

	if (direction == FFTW_FORWARD) {
		/* shift the origin */
        	cycleimage(fr, nx, ny, nx / 2, ny / 2);
        	cycleimage(fi, nx, ny, nx / 2, ny / 2);
	}

	/* clean up */
	fftwnd_destroy_plan(p);  
	free(fk);
}




