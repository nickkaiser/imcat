/*
 * fixedges.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "imlib/fits.h"
#include "fixedges.h"

typedef struct pixel {
	int	x, y, fixed;
	float	f;
} pixel;


int	getpixcoords(pixel *pix, int npix, float **f, int Nx, int Ny);
int	fixpix(pixel *pix, int npix, float **f, int Nx, int Ny);


void	fixedges(float **f, int Nx, int Ny)
{
	int	x, y;
	float	f0;

	for (y = 0; y < Ny; y++) {
		for (x = 0; x < Nx; x++) {
			if (f[y][x] != FLOAT_MAGIC) {
				f0 = f[y][x];
				while(x--) {
					f[y][x] = f0;
				}
				break;
			}
		}
		for (x = Nx - 1; x >= 0; x--) {
			if (f[y][x] != FLOAT_MAGIC) {
				f0 = f[y][x];
				while(x++ < Nx - 1) {
					f[y][x] = f0;
				}
				break;
			}
		}
	}
	for (x = 0; x < Nx; x++) {
		for (y = 0; y < Ny; y++) {
			if (f[y][x] != FLOAT_MAGIC) {
				f0 = f[y][x];
				while(y--) {
					f[y][x] = f0;
				}
				break;
			}
		}
		for (y = Ny - 1; y >= 0; y--) {
			if (f[y][x] != FLOAT_MAGIC) {
				f0 = f[y][x];
				while(y++ < Ny - 1) {
					f[y][x] = f0;
				}
				break;
			}
		}
	}
}


void	fixholes(float **f, int Nx, int Ny)
{
	int	x, y, nmagic, ipix;
	pixel	*pix;

	/* count the magic pixels in the original image */
	nmagic = 0;
	for (x = 0; x < Nx; x++) {
		for (y = 0; y < Ny; y++) {
			if (f[y][x] == FLOAT_MAGIC) {
				nmagic++;
			}
		}
	} 

	/* allocate magic pixels */
	pix = (pixel *) calloc(nmagic, sizeof(pixel));
	getpixcoords(pix, nmagic, f, Nx, Ny);
	while (nmagic = fixpix(pix, nmagic, f, Nx, Ny)) {
	}
	free(pix);
}

int	getpixcoords(pixel *pix, int npix, float **f, int Nx, int Ny)
{
	int	x, y, ipix;

	ipix = 0;
	for (x = 0; x < Nx; x++) {
		for (y = 0; y < Ny; y++) {
			if (f[y][x] == FLOAT_MAGIC) {
				pix[ipix].x = x;
				pix[ipix].y = y;
				ipix++;
				if (ipix == npix) {
					return(1);
				}
			}
		}
	}
	return(0);
}


int	fixpix(pixel *pix, int npix, float **f, int Nx, int Ny)
{
	int	x, y, p, pp;

	/* flag the pixels to be fixed */
	for (p = 0; p < npix; p++) {
		for (x = pix[p].x - 1; x <= pix[p].x + 1; x += 2) {
			for (y = pix[p].y - 1; y <= pix[p].y + 1; y += 2) {
				if (x >= 0 && x < Nx && y >= 0 && y < Ny) {
					if (f[y][x] != FLOAT_MAGIC) {
						pix[p].f = f[y][x];
						pix[p].fixed = 1;
						break;
					}
				}
			}
			if (pix[p].fixed) {
				break;
			}
		}
	}
	/* change the fixed pixels in the image */
	for (p = 0; p < npix; p++) {
		if (pix[p].fixed) {
			f[pix[p].y][pix[p].x] = pix[p].f;
		}
	}
	/* rearrange the pixel list */
	pp = 0;
	for (p = 0; p < npix; p++) {
		if (!pix[p].fixed) {
			pix[pp++] = pix[p];
		}
	}
	/* return the number of pix to be fixed */
	return(pp);
}


