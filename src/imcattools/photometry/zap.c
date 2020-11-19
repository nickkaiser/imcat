/*
 * zap.c ---- functions to zap/restore image
 */

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include "zap.h"
#include "../../imlib/fits.h"

#define MAGIC FLOAT_MAGIC

void	zap(int mode, double rmax, double *x, float **f, float **fzap, short **nzap, int N1, int N2)
{
	/* set fzap to MAGIC in disk and increment nzap */
	int	ix, iy, ix0, iy0, d;
	double	rr, rrmax;
	
	d = (int) ceil(rmax);
	rrmax = rmax * rmax;
	ix0 = (int) floor(0.5 + x[0]);
	iy0 = (int) floor(0.5 + x[1]);
	for (iy = iy0 - d; iy <= iy0 + d; iy++) {
		if (iy < 0 || iy >= N2)
			continue;
		for (ix = ix0 - d; ix <= ix0 + d; ix++) {
			if (ix < 0 || ix >= N1)
				continue;
			rr = (ix - ix0) * (ix - ix0) + (iy - iy0) * (iy - iy0);
			if (rr <= rrmax) {
				switch (mode) {
					case ZAP_MODE:
						fzap[iy][ix] = MAGIC;
						nzap[iy][ix]++;
						break;
					case UNZAP_MODE:
						if (nzap[iy][ix] == 1)
							fzap[iy][ix] = f[iy][ix];
						nzap[iy][ix]--;
						break;
					default:
						error_exit("zap: bad zapmode\n");
						break;
				}
			}
		}
	}
}



