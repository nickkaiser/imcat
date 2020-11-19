/*
 * myfft.c
 */

#include <stdio.h>
#include <limits.h>
#include <math.h>
#include "myfft.h"
#include "../utils/error.h"
#include "../imlib/fits.h"

#define PI M_PI

/* 
 * cycleimage() --- shifts the origin of a periodic image to x0, y0	
 */

void	cycleimage(float **f, int Nx, int Ny, int x0, int y0)
{
	int	x, y, xp, yp;
	float	*F;

	/* shift the origin into the image */
	while (x0 < 0) {
		x0 += Nx;
	}
	while (y0 < 0) {
		y0 += Ny;
	}
	
	F = (float *) calloc(Nx, sizeof(float));
	for (y = 0; y < Ny; y++) {
		for (x = 0; x < Nx; x++) {
			F[x] = f[y][x];
		}
		for (x = 0; x < Nx; x++) {
			f[y][(x + x0) % Nx] = F[x];
		}
	}
	free(F);
	F = (float *) calloc(Ny, sizeof(float));
	for (x = 0; x < Nx; x++) {
		for (y = 0; y < Ny; y++) {
			F[y] = f[y][x];
		}
		for (y = 0; y < Ny; y++) {
			f[(y + y0) % Ny][x] = F[y];
		}
	}
	free(F);
}


/* old style non-contiguous memory version
void	cycleimage(float **f, int Nx, int Ny, int x0, int y0)
{
	int	x, y;
	float	*F, **FF;

	if (x0 < 0 || x0 >= Nx || y0 < 0 || y0 >= Ny)
		error_exit("cycleimage: shift must be to point in image\n");

	F = (float *) calloc(Nx, sizeof(float));
	FF = (float **) calloc(Ny, sizeof(float *));
	for (y = 0; y < Ny; y++) {
		for (x = 0; x < Nx; x++) {
			F[x] = f[y][x];
		}
		for (x = 0; x < Nx - x0; x++) {
			f[y][x + x0] = F[x];
		}
		for (x = Nx - x0; x < Nx; x++) {
			f[y][x + x0 - Nx] = F[x];
		}
	}
	for (y = 0; y < Ny; y++) {
		FF[y] = f[y];
	}
	for (y = 0; y < Ny - y0; y++) {
		f[y + y0] = FF[y];
	}
	for (y = Ny - y0; y < Ny; y++) {
		f[y + y0 - Ny] = FF[y];
	}
	free(F);
	free(FF);
}
*/

void	substitute(float **f, int Nx, int Ny, float magicsubstitute)
{
	int	x, y;

	for (y = 0; y < Ny; y++) {
		for (x = 0; x < Nx; x++) {
			if (f[y][x] == (float) FLOAT_MAGIC)
				f[y][x] = magicsubstitute;
		}
	}
}
