/*
 *
 * fft_FFTPACK.c	- FFTPACK fft routines
 *
 */


#include <stdio.h>
#include <stdlib.h>
#ifndef Darwin
#include <malloc.h>
#endif
#include <math.h>

#ifdef DUAL_PROC
#include "shmalloc.h"
#include "../imtools/fft_task.h"
static	int	fr_id, fk_id;	/* shared memory */
#endif

/*
#include "../fftpack/fftpack.h"
*/

#include "myfft.h"

#ifndef PI
#define PI M_PI
#endif


void	forward_fft(float **f, int nx, int ny, fft_type fk)
{
	int 	x, y, n;
	float	*r, *work1, *work2;
	char	sysstring[128];

	n = ny / 2;

	/* y-direction fft's */
#ifdef DUAL_PROC
	sprintf(sysstring, "fft_task %d %d %d %d %d", fr_id, fk_id, nx, ny, FORWARD_Y);
	if (system(sysstring)) {
		fprintf(stderr, "forward_fft: fft_task failed: freeing shared memory...\n");
                shmfree(fr_id);
                shmfree(fk_id);
		exit(-1);
	}
#else
	r = (float *) calloc(ny, sizeof(float));
	work1 = (float *) calloc(2 * ny + 15, sizeof(float));
	rffti_(&ny, work1);
	for (x = 0; x < nx; x++) {
		for (y = 0; y < ny; y++) {
			r[y] = f[y][x];
		}
		rfftf_(&ny, r, work1);
		fk[0][x].r = r[0];
		fk[0][x].i = 0.0;
		for (y = 1; y < n; y++) {
			fk[y][x].r = r[2 * y - 1];
			fk[y][x].i = r[2 * y];
		}
		fk[n][x].r = r[ny - 1];
		fk[n][x].i = 0.0;
	}
	free(work1);
	free(r);
#endif

	/* x-direction fft's */
#ifdef DUAL_PROC
	sprintf(sysstring, "fft_task %d %d %d %d %d", fr_id, fk_id, nx, ny, FORWARD_X);
	if (system(sysstring)) {
		fprintf(stderr, "forward_fft: fft_task failed: freeing shared memory...\n");
                shmfree(fr_id);
                shmfree(fk_id);
		exit(-1);
	}
#else
	work2 = (float *) calloc(4 * nx + 15, sizeof(float));
	cffti_(&nx, work2);
	for (y = 0; y <= n; y++) {
		cfftf_(&nx, fk[y], work2);
	}
	free(work2);
#endif
}

void	inverse_fft(fft_type fk, int nx, int ny, float **f) 
{
	int 	x, y, n;
	float	*r, *work1, *work2;
	char	sysstring[128];

	n = ny / 2;

	/* inverse x-transforms */
#ifdef DUAL_PROC
	sprintf(sysstring, "fft_task %d %d %d %d %d", fr_id, fk_id, nx, ny, INVERSE_X);
	if (system(sysstring)) {
		fprintf(stderr, "inverse_fft: fft_task failed: freeing shared memory...\n");
                shmfree(fr_id);
                shmfree(fk_id);
		exit(-1);
	}
#else
	work2 = (float *) calloc(4 * nx + 15, sizeof(float));
	cffti_(&nx, work2);
	for (y = 0; y <= n; y++) {
		cfftb_(&nx, fk[y], work2);
	}
	free(work2);
#endif

	/* inverse y-transforms  */
#ifdef DUAL_PROC
	sprintf(sysstring, "fft_task %d %d %d %d %d", fr_id, fk_id, nx, ny, INVERSE_Y);
	if (system(sysstring)) {
		fprintf(stderr, "inverse_fft: fft_task failed: freeing shared memory...\n");
                shmfree(fr_id);
                shmfree(fk_id);
		exit(-1);
	}
#else
	r = (float *) calloc(ny, sizeof(float));
	work1 = (float *) calloc(2 * ny + 15, sizeof(float));
	rffti_(&ny, work1);
	for (x = 0; x < nx; x++) {
		r[0] = fk[0][x].r;
		for (y = 1; y < n; y++) {
			r[2 * y - 1] = fk[y][x].r;
			r[2 * y] = fk[y][x].i;
		}
		r[ny - 1] = fk[n][x].r;
		rfftb_(&ny, r, work1);
		for (y = 0; y < ny; y++) {
			 f[y][x] = r[y] / (nx * ny);
		}
	}
	free(r);
	free(work1);
#endif
}


void	filter(fft_type fk, int Nx, int Ny, float (*filterfunc)(float ki, float kj))
{
	int	ix, iy, n;
	float	T, kx, ky;
	complex **Fk;

	Fk = (complex **) fk;
	n = Ny / 2;
	for (iy = 0; iy <= n; iy++) {
		ky = 2 * PI * iy / (float) Ny;
		for (ix = 0; ix < Nx; ix++) {
			kx = (ix < Nx / 2 ? ix : ix - Nx) * 2 * PI / (float) Nx;
			T = filterfunc(kx, ky);
			Fk[iy][ix].r *= T;
			Fk[iy][ix].i *= T;
		}
	}
}


void	cfilter(fft_type fk, int Nx, int Ny, float (*rfunc)(float ki, float kj), float (*ifunc)(float ki, float kj))
{
	int	ix, iy, n;
	float	Tr, Ti, fr, fi, kx, ky;
	complex **Fk;

	Fk = (complex **) fk;
	n = Ny / 2;
	for (iy = 0; iy <= n; iy++) {
		ky = 2 * PI * iy / (float) Ny;
		for (ix = 0; ix < Nx; ix++) {
			kx = (ix < Nx / 2 ? ix : ix - Nx) * 2 * PI / (float) Nx;
			Tr = rfunc(kx, ky);
			Ti = ifunc(kx, ky);
			/* funny signs below are for compatibility */
			/* with NR convention. Our Fk is really conjugate */
			fr = Fk[iy][ix].r;
			fi = -Fk[iy][ix].i;
			Fk[iy][ix].r = fr * Tr - fi * Ti;
			Fk[iy][ix].i = -1.0 * (fr * Ti + fi * Tr);			
		}
	}
}



void	ccf(fft_type fk1, fft_type fk2, int Nx, int Ny, float **ff12, int x0, int y0)
{
	int	ix, iy, n;
	float	a1, b1, a2, b2;
	complex **Fk1, **Fk2;

	Fk1 = (complex **) fk1;
	Fk2 = (complex **) fk2;
	n = Ny / 2;
	for (iy = 0; iy <= n; iy++) {
		for (ix = 0; ix < Nx; ix ++) {
			a1 = Fk1[iy][ix].r;
			a2 = Fk2[iy][ix].r;
			b1 = Fk1[iy][ix].i;
			b2 = Fk2[iy][ix].i;
			Fk1[iy][ix].r = (a1 * a2 + b1 * b2) / (Nx * Ny);
			Fk1[iy][ix].i = (a1 * b2 - b1 * a2) / (Nx * Ny);
		}
	}
	inverse_fft(fk1, Nx, Ny, ff12);
	cycleimage(ff12, Nx, Ny, x0, y0);
}



void	power(fft_type fk, int Nx, int Ny, float **p, int x0, int y0)
{
	int	ix, iy, n, x, y;
	float	a, b;
	complex **Fk;

	Fk = (complex **) fk;
	n = Ny / 2;
	for (iy = 0; iy <= n; iy++) {
		for (ix = 0; ix < Nx; ix ++) {
			a = Fk[iy][ix].r;
			b = Fk[iy][ix].i;
			p[iy][ix] = (a * a + b * b) / (Nx * Ny);
		}
	}
	cycleimage(p, Nx, Ny, Nx / 2, Ny / 2);
/* 	we now one half filled so now make other half by reflection */
	for (iy = n + 1; iy < Ny; iy++) {
		y = Ny - iy;
		for (ix = 0; ix < Nx; ix++) {
			x = Nx - ix;
			if (y >= 0 && y < Ny && x >= 0 && x < Nx)
				p[y][x] = p[iy][ix];
		}
	}
	cycleimage(p, Nx, Ny, Nx / 2, Ny / 2);
	cycleimage(p, Nx, Ny, x0, y0);
}




void	alloc_fft(fft_type *fk, int Nx, int Ny)
{
	int 	n, y;

	if (Ny % 2) {
		error_exit("alloc_fft: Ny must be even\n");
	}

	n = Ny / 2;
	*fk = (complex **) calloc(n + 1, sizeof(complex *));
	(*fk)[0] = (complex *) calloc(Nx * (n + 1), sizeof(complex));
	if (!(*fk) || !((*fk)[0]))
		error_exit("alloc_fft: failed to allocate memory\n");
	for (y = 1; y <= n; y++) {
		(*fk)[y] = (*fk)[y - 1] + Nx;	
	}
}



void	free_fft(fft_type fk, int Nx, int Ny)
{
	free(fk[0]);
	free(fk);
}


#ifdef DUAL_PROC

int	alloc_fft_shm(fft_type *fk, int Nx, int Ny)
{
	int 	n, y, shmid;

	if (Ny % 2) {
		error_exit("alloc_fft: Ny must be even\n");
	}

	n = Ny / 2;
	*fk = (complex **) calloc(n + 1, sizeof(complex *));
	shmid = shmalloc((void **) &((*fk)[0]), Nx * (n + 1) * sizeof(complex));
	if (!(*fk) || !((*fk)[0]))
		error_exit("alloc_fft: failed to allocate memory\n");
	for (y = 1; y <= n; y++) {
		(*fk)[y] = (*fk)[y - 1] + Nx;	
	}
	return(shmid);
}



void	free_fft_shm(fft_type fk, int shmid)
{
	shmfree(shmid);
	free(fk);
}


void	set_shmids(int fr_shmid, int fk_shmid)
{	
	fr_id = fr_shmid;
	fk_id = fk_shmid;
}

void	get_shmids(int *fr_shmid, int *fk_shmid)
{
	*fr_shmid  = fr_id;	
	*fk_shmid  = fk_id;
}

void	set_frshmid(int fr_shmid)
{
	fr_id = fr_shmid;
}


void	set_fkshmid(int fk_shmid)
{
	fk_id = fk_shmid;
}


#endif

void	copy_fft(fft_type fksrc, int Nx, int Ny, fft_type fkdest)
{
	int	x, y, n;

	n = Ny / 2;
	for (y = 0; y <= n; y++) {
		for (x = 0; x < Nx; x++) {
			fkdest[y][x] = fksrc[y][x];
		}
	}
}


void	get_fft(fft_type fk, int Nx, int Ny, float **fdest)
{
	int	kx, ky, ix, iy, ixx, n;
	float	**fdestreal, **fdestimag;

	fdestreal = fdest;
	fdestimag = fdest + Ny;
	n = Ny / 2;
	for (ix = 0; ix < Nx; ix ++) {
		kx = (ix < Nx / 2 ? ix : ix - Nx);
		for (iy = 0; iy < Ny / 2; iy++) {
			ky = iy;
			fdestreal[Ny / 2 + ky][Nx / 2 + kx] = fk[iy][ix].r;
			fdestimag[Ny / 2 + ky][Nx / 2 + kx] = -fk[iy][ix].i;
			if (ky) {
				ixx = ((Nx / 2 - kx) < Nx ? Nx / 2 - kx : 0);
				fdestreal[Ny / 2 - ky][ixx] = fk[iy][ix].r;
				fdestimag[Ny / 2 - ky][ixx] = fk[iy][ix].i;
			}
		}
		fdestreal[0][Nx / 2 + kx] = fk[iy][ix].r;
		fdestimag[0][Nx / 2 + kx] = -fk[iy][ix].i;
	}
}


void	set_fft(fft_type fk, int Nx, int Ny, float **fsrc)
{
	int	kx, ky, ix, iy, iyy, n;
	float	**freal, **fimag;

	freal = fsrc;
	fimag = fsrc + Ny;
	n = Ny / 2;
	for (ix = 0; ix < Nx; ix ++) {
		kx = (ix < Nx / 2 ? ix : ix - Nx);
		for (iy = 0; iy < Ny / 2; iy++) {
			ky = iy;
			fk[iy][ix].r = freal[Ny / 2 + ky][Nx / 2 + kx];
			fk[iy][ix].i = -fimag[Ny / 2 + ky][Nx / 2 + kx];
		}
		fk[iy][ix].r = freal[0][Nx / 2 + kx];
		fk[iy][ix].i = -fimag[0][Nx / 2 + kx];
	}
}


void	forward_cfft(float **fsrc, int nx, int ny, float **fdst)
{
	int 	x, y;
	float	*work, **fsr, **fsi, **fdr, **fdi;
	complex	*c;

	fsr = fsrc;
	fsi = fsrc + ny;
	fdr = fdst;
	fdi = fdst + ny;

	/* y-direction fft's */
	c = (complex *) calloc(ny, sizeof(complex));
	work = (float *) calloc(4 * ny + 15, sizeof(float));
	cffti_(&ny, work);
	for (x = 0; x < nx; x++) {
		for (y = 0; y < ny; y++) {
			c[y].r = fsr[y][x];
			c[y].i = fsi[y][x];
		}
		cfftf_(&ny, c, work);
		for (y = 0; y < ny; y++) {
			fdr[y][x] = c[y].r;
			fdi[y][x] = c[y].i;
		}
	}
	free(work);
	free(c);

	/* x-direction fft's */
	c = (complex *) calloc(nx, sizeof(complex));
	work = (float *) calloc(4 * nx + 15, sizeof(float));
	cffti_(&nx, work);
	for (y = 0; y < ny; y++) {
		for (x = 0; x < nx; x++) {
			c[x].r = fdr[y][x];
			c[x].i = fdi[y][x];
		}
		cfftf_(&nx, c, work);
		for (x = 0; x < nx; x++) {
			fdr[y][x] = c[x].r;
			fdi[y][x] = c[x].i;
		}
	}
	free(work);
	free(c);
	cycleimage(fdr, nx, ny, nx / 2, ny / 2);
	cycleimage(fdi, nx, ny, nx / 2, ny / 2);
}


void	inverse_cfft(float **fsrc, int nx, int ny, float **fdst) 
{
	int 	x, y;
	float	*work, **fsr, **fsi, **fdr, **fdi;
	complex	*c;

	fsr = fsrc;
	fsi = fsrc + ny;
	fdr = fdst;
	fdi = fdst + ny;

	cycleimage(fsr, nx, ny, -nx / 2, -ny / 2);
	cycleimage(fsi, nx, ny, -nx / 2, -ny / 2);

	/* inverse x-transforms */
	work = (float *) calloc(4 * nx + 15, sizeof(float));
	c = (complex *) calloc(nx, sizeof(complex));
	cffti_(&nx, work);
	for (y = 0; y < ny; y++) {
		for (x = 0; x < nx; x++) {
			c[x].r = fsr[y][x];
			c[x].i = fsi[y][x];
		}
		cfftb_(&nx, c, work);
		for (x = 0; x < nx; x++) {
			fdr[y][x] = c[x].r;
			fdi[y][x] = c[x].i;
		}
	}
	free(work);
	free(c);

	/* inverse y-transforms  */
	c = (complex *) calloc(ny, sizeof(complex));
	work = (float *) calloc(4 * ny + 15, sizeof(float));
	cffti_(&ny, work);
	for (x = 0; x < nx; x++) {
		for (y = 0; y < ny; y++) {
			c[y].r = fdr[y][x];
			c[y].i = fdi[y][x];
		}
		cfftb_(&ny, c, work);
		for (y = 0; y < ny; y++) {
			 fdr[y][x] = c[y].r / (nx * ny);
			 fdi[y][x] = c[y].i / (nx * ny);
		}
	}
	free(c);
	free(work);
}


