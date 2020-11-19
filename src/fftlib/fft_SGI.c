/*
 *
 * fft_sgi.c	- SGI fft routines
 *
 */

/*
 * SGI version uses unpacked notation so total size is ntotal = ldx * (Ny + 1)
 * with leading (i.e. fast) dimension ldx = 2 * ((Nx + 2) / 2)
 *
 */

#include <stdio.h>
#include <math.h>

#include "../../lib/sgi_fft/sgi_fft.h"
#include "myfft.h"

void	forward_fft(float **f, int Nx, int Ny, fft_type fk)
{
	int	x, y, ldx, ntotal;
	float	*work;

	ldx = 2 * ((Nx + 2) / 2);
	ntotal = ldx * (Ny + 1);
	work = 	(float *) calloc(Nx + 3 * Ny + 4 * FACTOR_SPACE, sizeof(float));

	work = sfft2dui(Nx, Ny, work);
	for (y = 0; y < Ny + 1; y++) {
		for (x = 0; x < ldx; x++) {
			if (0 <= y && y < Ny && 0 <= x && x <Nx)
 				((float *)fk)[x + y * ldx] = f[y][x];
			else
 				((float *)fk)[x + y * ldx] = 0;
		}
	}
	sfft2du(-1, Nx, Ny, (float *)fk, ldx, work);
	free(work);
}

void	inverse_fft(fft_type fk, int Nx, int Ny, float **f) 
{
	int	x, y, ldx;
	float	*work, *Fk;

	Fk = (float *) fk;
	ldx = 2 * ((Nx + 2) / 2);
	work = 	(float *) calloc(Nx + 3 * Ny + 4 * FACTOR_SPACE, sizeof(float));

	work = sfft2dui(Nx, Ny, work);
	sfft2du(1, Nx, Ny, Fk, ldx, work);
	for (y = 0; y < Ny; y++) {
		for (x = 0; x < Nx; x++) {
 			f[y][x] = Fk[x + y * ldx] / (float) (Nx * Ny);
		}
	}
	free(work);
}


void	filter(fft_type fk, int Nx, int Ny, float (*filterfunc)(float ki, float kj))
{
	int	ix, iy, ldx;
	float	T, kx, ky, *Fk;

	Fk = (float *) fk;
	ldx = 2 * ((Nx + 2) / 2);
	for (iy = 0; iy < Ny; iy++) {
		ky = (iy < Ny / 2 ? iy : iy - Ny) * 2 * PI / (float) Ny;
		for (ix = 0; ix < ldx / 2; ix ++) {
			kx = 2 * PI * ix / (float) Nx;
			T = filterfunc(kx, ky);
			Fk[iy * ldx + 2 * ix] *= T;
			Fk[iy * ldx + 2 * ix + 1] *= T;
		}
	}
}


void	cfilter(fft_type fk, int Nx, int Ny, float (*rfunc)(float ki, float kj), float (*ifunc)(float ki, float kj))
{
	int	ix, iy, ldx;
	float	Tr, Ti, fr, fi, kx, ky, *Fk;

	Fk = (float *) fk;
	ldx = 2 * ((Nx + 2) / 2);
	for (iy = 0; iy < Ny; iy++) {
		ky = (iy < Ny / 2 ? iy : iy - Ny) * 2 * PI / (float) Ny;
		for (ix = 0; ix < ldx / 2; ix ++) {
			kx = 2 * PI * ix / (float) Nx;
			Tr = rfunc(kx, ky);
			Ti = ifunc(kx, ky);
			fr = Fk[iy * ldx + 2 * ix];
			fi = Fk[iy * ldx + 2 * ix + 1];
			Fk[iy * ldx + 2 * ix] = fr * Tr - fi * Ti;
			Fk[iy * ldx + 2 * ix + 1] = fr * Ti + fi * Tr;			
		}
	}
}



void	ccf(fft_type fk1, fft_type fk2, int Nx, int Ny, float **ff12, int x0, int y0)
{
	int	ix, iy, ldx;
	float	a1, b1, a2, b2, *Fk1, *Fk2;

	Fk1 = (float *) fk1;
	Fk2 = (float *) fk2;
	ldx = 2 * ((Nx + 2) / 2);
	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < ldx / 2; ix ++) {
			a1 = Fk1[iy * ldx + 2 * ix];
			a2 = Fk2[iy * ldx + 2 * ix];
			b1 = Fk1[iy * ldx + 2 * ix + 1];
			b2 = Fk2[iy * ldx + 2 * ix + 1];
			Fk1[iy * ldx + 2 * ix] = 	(a1 * a2 + b1 * b2) / (Nx * Ny);
			Fk1[iy * ldx + 2 * ix + 1] = 	(a1 * b2 - b1 * a2) / (Nx * Ny);
		}
	}
	inverse_fft(fk1, Nx, Ny, ff12);
	cycleimage(ff12, Nx, Ny, x0, y0);
}



void	power(fft_type fk, int Nx, int Ny, float **p, int x0, int y0)
{
	int	ix, iy, ldx, x, y;
	float	a, b, *Fk;

	Fk = (float *) fk;
	ldx = 2 * ((Nx + 2) / 2);
	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < ldx / 2; ix ++) {
			a = Fk[iy * ldx + 2 * ix];
			b = Fk[iy * ldx + 2 * ix + 1];
			p[iy][ix] = (a * a + b * b) / (Nx * Ny);
		}
	}
	cycleimage(p, Nx, Ny, Nx / 2, Ny / 2);
/* 	we now have right half filled so now make other half by reflection */
	for (iy = 0; iy < Ny; iy++) {
		y = Ny - iy;
		for (ix = Nx / 2 + 1; ix < Nx; ix++) {
			x = Nx - ix;
			if (x >= 0 && x < Nx && y >= 0 && y < Ny)
				p[y][x] = p[iy][ix];
		}
	}
	cycleimage(p, Nx, Ny, Nx / 2, Ny / 2);
	cycleimage(p, Nx, Ny, x0, y0);
}




void	alloc_fft(fft_type *fk, int Nx, int Ny)
{
	int	ldx, ntotal;
	float	*Fk;

	ldx = 2 * ((Nx + 2) / 2);
	ntotal = ldx * (Ny + 1);
	Fk = (float *) calloc(ntotal, sizeof(float));
	*fk = (fft_type) Fk;
}



void	free_fft(fft_type fk, int Nx, int Ny)
{
	free(fk);
}




void	copy_fft(fft_type fksrc, int Nx, int Ny, fft_type fkdest)
{
	int	x, y, ldx;

	ldx = 2 * ((Nx + 2) / 2);
	for (y = 0; y < (Ny + 1); y++) {
		for (x = 0; x < ldx; x++) {
			fkdest[y * ldx + x] = fksrc[y * ldx + x];
		}
	}
}


/* minus signs on fimag below to make values agree with NR values */

void	get_fft(fft_type fk, int Nx, int Ny, float **fdest)
{
	int	kx, ky, ix, iy, iyy, ldx;
	float	**fdestreal, **fdestimag;

	fdestreal = fdest;
	fdestimag = fdest + Ny;
	ldx = 2 * ((Nx + 2) / 2);
	for (iy = 0; iy < Ny; iy++) {
		ky = (iy < Ny / 2 ? iy : iy - Ny);
		for (ix = 0; ix < Nx / 2; ix ++) {
			kx = ix;
			fdestreal[Ny / 2 + ky][Nx / 2 + kx] = fk[iy * ldx + 2 * ix];
			fdestimag[Ny / 2 + ky][Nx / 2 + kx] = -fk[iy * ldx + 2 * ix + 1];
			if (kx) {
				iyy = ((Ny / 2 - ky) < Ny ? Ny / 2 - ky : 0);
				fdestreal[iyy][Nx / 2 - kx] = fk[iy * ldx + 2 * ix];
				fdestimag[iyy][Nx / 2 - kx] = fk[iy * ldx + 2 * ix + 1];
			}
		}
		fdestreal[Ny / 2 + ky][0] = fk[iy * ldx + 2 * ix];
		fdestimag[Ny / 2 + ky][0] = -fk[iy * ldx + 2 * ix + 1];
	}
}


void	set_fft(fft_type fk, int Nx, int Ny, float **fsrc)
{
	int	kx, ky, ix, iy, iyy, ldx;
	float	**freal, **fimag;

	freal = fsrc;
	fimag = fsrc + Ny;
	ldx = 2 * ((Nx + 2) / 2);
	for (iy = 0; iy < Ny; iy++) {
		ky = (iy < Ny / 2 ? iy : iy - Ny);
		for (ix = 0; ix < Nx / 2; ix ++) {
			kx = ix;
			fk[iy * ldx + 2 * ix] = freal[Ny / 2 + ky][Nx / 2 + kx];
			fk[iy * ldx + 2 * ix + 1] = -fimag[Ny / 2 + ky][Nx / 2 + kx];
		}
		fk[iy * ldx + 2 * ix] = freal[Ny / 2 + ky][0];
		fk[iy * ldx + 2 * ix + 1] = -fimag[Ny / 2 + ky][0];
	}
}


