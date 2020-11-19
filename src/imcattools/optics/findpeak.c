#include <math.h>
#include "findpeak.h"

#define TINY	1.e-20


int	findpeak(float **f, int N1, int N2, int *ixpk, int *iypk, float *xpk, float *ypk, float *fpk)
{
	int	ix, iy, x, y;
	float	fmax, Fp, Fpp;
	
	/* find the hottest pixel */
	fmax = 0.0;
	for (y = 0; y < N2; y++) {
		for (x = 0; x < N1; x++) {
			if (f[y][x] > fmax) {
				fmax = f[y][x];
				ix = x;
				iy = y;
			}
		}
	}
	*ixpk = ix;
	*iypk = iy;
	*xpk = 0.5 + (float) ix;
	*ypk = 0.5 + (float) iy;
	*fpk = fmax;
	/* refine the peak */
	if (ix > 0 && ix < N1 - 1 && iy > 0 && iy < N2 - 1) {
		Fp  = 0.5 * (f[iy][ix+1] - f[iy][ix-1]);
        	Fpp = f[iy][ix+1] - 2 * f[iy][ix] + f[iy][ix-1];
                if (fabs(Fpp) > TINY) {
                	*xpk -= Fp / Fpp;
		}
                Fp  = 0.5 * (f[iy+1][ix] - f[iy-1][ix]);
                Fpp = f[iy+1][ix] - 2 * f[iy][ix] + f[iy-1][ix];
                if (fabs(Fpp) > TINY) {
                	*ypk -= Fp / Fpp;
		}
	}
	return(1);
}

