/*
 *
 * fft_nr.c	- numerical recipes fft routines
 *
 */


#include <stdio.h>
#include <math.h>

#include "../utils/nrutil.h"
#include "myfft.h"


#define PI M_PI

void 	rlft3(float ***data, float **speq, unsigned long nn1, unsigned long nn2,
        	unsigned long nn3, int isign);
void 	fft_size(int Nx, int Ny, int *NFFTx, int *NFFTy);
void	rotate(float *a, float *b, float phi);


void	forward_fft(float **f, int Nx, int Ny, fft_type fk)
{
	int	x, y, NFFTx, NFFTy;

	fft_size(Nx, Ny, &NFFTx, &NFFTy);
	for (y = 0; y < NFFTy; y++) {
		for (x = 0; x < NFFTx; x++) {
			if (0 <= y && y < Ny && 0 <= x && x <Nx)
 				fk->data[1][y+1][x+1] = f[y][x];
			else
 				fk->data[1][y+1][x+1] = 0;
		}
		fk->speq[1][2 * (y+1) - 1] = fk->speq[1][2 * (y+1)] = 0.0;
	}
	rlft3(fk->data, fk->speq, 1, NFFTy, NFFTx, 1);
}

void	inverse_fft(fft_type fk, int Nx, int Ny, float **f) 
{
	int	x, y, NFFTx, NFFTy;

	fft_size(Nx, Ny, &NFFTx, &NFFTy);
	rlft3(fk->data, fk->speq, 1, NFFTy, NFFTx, -1);
	for (y = 0; y < Ny; y++) {
		for (x = 0; x < Nx; x++) {
			f[y][x] = 2 * fk->data[1][y+1][x+1] / (NFFTx * NFFTy);
		}
	}
}


void	filter(fft_type fk, int Nx, int Ny, float (*filterfunc)(float ki, float kj))
{
	int	i2, i3, NFFTx, NFFTy;
	float	T, kx, ky;

	fft_size(Nx, Ny, &NFFTx, &NFFTy);
	for (i2 = 1; i2 <= NFFTy; i2++) {
		ky = ((i2-1) < NFFTy / 2 ? (i2-1) : (i2-1) - NFFTy) * 2 * PI / (float) NFFTy;
		for (i3 = 1; i3 <= NFFTx / 2; i3++) {
			kx = 2 * PI * (i3-1) / (float) NFFTx;
			T = filterfunc(kx, ky);
			fk->data[1][i2][2 * i3 - 1] *= T;
			fk->data[1][i2][2 * i3] *= T;
		}
		kx = PI;
		T = filterfunc(kx, ky);
		fk->speq[1][2 * i2 - 1] *= T;
		fk->speq[1][2 * i2] *= T;
	}
}


void	cfilter(fft_type fk, int Nx, int Ny, float (*rfunc)(float ki, float kj),  float (*ifunc)(float ki, float kj))
{
	int	i2, i3, NFFTx, NFFTy;
	float	T, Tr, Ti, kx, ky, fr, fi;

	fft_size(Nx, Ny, &NFFTx, &NFFTy);
	for (i2 = 1; i2 <= NFFTy; i2++) {
		ky = ((i2-1) < NFFTy / 2 ? (i2-1) : (i2-1) - NFFTy) * 2 * PI / (float) NFFTy;
		for (i3 = 1; i3 <= NFFTx / 2; i3++) {
			kx = 2 * PI * (i3-1) / (float) NFFTx;
			Tr = rfunc(kx, ky);
			Ti = ifunc(kx, ky);
			fr = fk->data[1][i2][2 * i3 - 1];
			fi = fk->data[1][i2][2 * i3];
			fk->data[1][i2][2 * i3 - 1] = (fr * Tr - fi * Ti);
			fk->data[1][i2][2 * i3] = (fi * Tr + fr * Ti);
		}
		kx = PI;
		Tr = rfunc(kx, ky);
		Ti = ifunc(kx, ky);
		fr = fk->speq[1][2 * i2 - 1];
		fi = fk->speq[1][2 * i2];
		fk->speq[1][2 * i2 - 1] = (fr * Tr - fi * Ti);
		fk->speq[1][2 * i2] = (fi * Tr + fr * Ti);
	}
}



void	ccf(fft_type fk1, fft_type fk2, int Nx, int Ny, float **ff12, int x0, int y0)
{
	int	i2, i3, NFFTx, NFFTy;
	float	a1, b1, a2, b2, kx, ky, phi;

	fft_size(Nx, Ny, &NFFTx, &NFFTy);
	for (i2 = 1; i2 <= NFFTy; i2++) {
		ky = ((i2-1) < NFFTy / 2 ? (i2-1) : (i2-1) - NFFTy) * 2 * PI / (float) NFFTy;
		for (i3 = 1; i3 <= NFFTx / 2; i3++) {
			kx = 2 * PI * (i3-1) / (float) NFFTx;
			phi = x0 * kx + y0 * ky;
			a1 = fk1->data[1][i2][2 * i3 - 1];
			b1 = fk1->data[1][i2][2 * i3];
			a2 = fk2->data[1][i2][2 * i3 - 1];
			b2 = fk2->data[1][i2][2 * i3];
			fk1->data[1][i2][2 * i3 - 1] = 	(a1 * a2 + b1 * b2) / (NFFTx * NFFTy);
			fk1->data[1][i2][2 * i3] = 	(a1 * b2 - b1 * a2) / (NFFTx * NFFTy);
			rotate(&(fk1->data[1][i2][2 * i3 - 1]), &(fk1->data[1][i2][2 * i3]), phi);  
		}
		i3 = NFFTx / 2 + 1;
		kx = 2 * PI * (i3-1) / (float) NFFTx;
		phi = x0 * kx + y0 * ky;
		a1 = fk1->speq[1][2 * i2 - 1];
		b1 = fk1->speq[1][2 * i2];
		a2 = fk2->speq[1][2 * i2 - 1];
		b2 = fk2->speq[1][2 * i2];
		fk1->speq[1][2 * i2 - 1] = (a1 * a2 + b1 * b2) / (Nx * Ny);
		fk1->speq[1][2 * i2] = (a1 * b2 - b1 * a2) / (Nx * Ny);
		rotate(&(fk1->speq[1][2 * i2 - 1]), &(fk1->speq[1][2 * i2]), phi);  
	}
	inverse_fft(fk1, Nx, Ny, ff12);
}



void	power(fft_type fk, int Nx, int Ny, float **p, int x0, int y0)
{
	int	i2, i3, NFFTx, NFFTy, kx, ky, x, y, yy;
	float	a, b;

	fft_size(Nx, Ny, &NFFTx, &NFFTy);
	for (i2 = 1; i2 <= NFFTy; i2++) {
		ky = ((i2-1) < NFFTy / 2 ? (i2-1) : (i2-1) - NFFTy);
		y = Ny / 2 + ky;
		if (y < 0 || y >= Ny)
			continue;
		for (i3 = 1; i3 <= NFFTx / 2; i3++) {
			kx = i3 - 1;
			a = fk->data[1][i2][2 * i3 - 1];
			b = fk->data[1][i2][2 * i3];
			x = Nx / 2 + kx;
			if (x >= 0 && x < Nx)
				p[y][x] = (a * a + b * b) / (Nx * Ny);
			x = Nx / 2 - kx;
			yy = Ny / 2 - ky;
			if (x >= 0 && x < Nx / 2)
				if (yy >= 0 && yy < Ny)
					p[yy][x] = (a * a + b * b) / (Nx * Ny);
		}
	}
	/* now we have P(k) centred on the grid */
	cycleimage(p, Nx, Ny, Nx / 2, Ny / 2);
	cycleimage(p, Nx, Ny, x0, y0);
}


void	free_fft(fft_type fk, int Nx, int Ny)
{
	int	NFFTx, NFFTy;

	fft_size(Nx, Ny, &NFFTx, &NFFTy);
	free_matrix(fk->speq, 1, 1, 1, 2 * NFFTy);
	free_f3tensor(fk->data, 1, 1, 1, NFFTy, 1, NFFTx);
	free(fk);
}



void 	fft_size(int Nx, int Ny, int *NFFTx, int *NFFTy)
{
	*NFFTx = 1;
	while (*NFFTx < Nx)
		(*NFFTx) *= 2;
	*NFFTy = 1;
	while (*NFFTy < Ny)
		(*NFFTy) *= 2;
}



void	rotate(float *a, float *b, float phi)
{
	float	A, B, c, s;

	c = cos(phi);
	s = sin(-phi);
	A = (c * *a + s * *b);
	B = (c * *b - s * *a);
	*a = A;
	*b = B;
}


void	alloc_fft(fft_type *fk, int Nx, int Ny)
{
	int	NFFTx, NFFTy;

	fft_size(Nx, Ny, &NFFTx, &NFFTy);
	(*fk) = (transform *) calloc(1, sizeof(transform));
	(*fk)->data = f3tensor(1, 1, 1, NFFTy, 1, NFFTx);
	(*fk)->speq = matrix(1, 1, 1, 2*NFFTy);
}


void	copy_fft(fft_type fksrc, int Nx, int Ny, fft_type fkdest)
{
	int	i2, i3, NFFTx, NFFTy;

	fft_size(Nx, Ny, &NFFTx, &NFFTy);
	for (i2 = 1; i2 <= NFFTy; i2++) {
		for (i3 = 1; i3 <= NFFTx / 2; i3++) {
			fkdest->data[1][i2][2 * i3 - 1] = fksrc->data[1][i2][2 * i3 - 1];
			fkdest->data[1][i2][2 * i3] = fksrc->data[1][i2][2 * i3];
		}
		fkdest->speq[1][2 * i2 - 1] = fksrc->speq[1][2 * i2 - 1];
		fkdest->speq[1][2 * i2] = fksrc->speq[1][2 * i2];
	}
}



void	get_fft(fft_type fk, int Nx, int Ny, float **fdest)
{
	int	i2, i3, NFFTx, NFFTy, kx, ky, ix, iy;
	float	fr, fi, **fdestreal, **fdestimag;

	fft_size(Nx, Ny, &NFFTx, &NFFTy);
	if ((Nx != NFFTx) || (Ny != NFFTy)) {
		fprintf(stderr, "get_fft: image dimension != 2^N: bailing out\n");
		exit(-1);
	}
	fdestreal = fdest;
	fdestimag = fdest + Ny;
	for (i2 = 1; i2 <= NFFTy; i2++) {
		ky = ((i2-1) < NFFTy / 2 ? (i2-1) : (i2-1) - NFFTy);
		for (i3 = 1; i3 <= NFFTx / 2; i3++) {
			kx = (i3-1);
			fdestreal[Ny / 2 + ky][Nx / 2 + kx] = fk->data[1][i2][2 * i3 - 1];
			fdestimag[Ny / 2 + ky][Nx / 2 + kx] = fk->data[1][i2][2 * i3];
			if (kx) {
				iy = ((Ny / 2 - ky) < Ny ? Ny / 2 - ky : 0);
				fdestreal[iy][Nx / 2 - kx] = fk->data[1][i2][2 * i3 - 1];
				fdestimag[iy][Nx / 2 - kx] = - fk->data[1][i2][2 * i3];
			}
		}
		fdestreal[Ny / 2 + ky][0] = fk->speq[1][2 * i2 - 1];
		fdestimag[Ny / 2 + ky][0] = fk->speq[1][2 * i2];
	}
}





void	set_fft(fft_type fk, int Nx, int Ny, float **fsrc)
{
	int	i2, i3, NFFTx, NFFTy, kx, ky, ix, iy;
	float	fr, fi, **fsrcreal, **fsrcimag;

	fft_size(Nx, Ny, &NFFTx, &NFFTy);
	if ((Nx != NFFTx) || (Ny != NFFTy)) {
		fprintf(stderr, "get_fft: image dimension != 2^N: bailing out\n");
		exit(-1);
	}
	fsrcreal = fsrc;
	fsrcimag = fsrc + Ny;
	for (i2 = 1; i2 <= NFFTy; i2++) {
		ky = ((i2-1) < NFFTy / 2 ? (i2-1) : (i2-1) - NFFTy);
		for (i3 = 1; i3 <= NFFTx / 2; i3++) {
			kx = (i3-1);
			fk->data[1][i2][2 * i3 - 1] = fsrcreal[Ny / 2 + ky][Nx / 2 + kx];
			fk->data[1][i2][2 * i3] = fsrcimag[Ny / 2 + ky][Nx / 2 + kx];
		}
		fk->speq[1][2 * i2 - 1] = fsrcreal[Ny / 2 + ky][0];
		fk->speq[1][2 * i2] = fsrcimag[Ny / 2 + ky][0];
	}
}


