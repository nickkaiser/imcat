/**
 ** scrunch_stream
 **	- works in stream mode 
 **	- tries to patch bad pixels
 **/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
 
#include	"fits.h"
#include	"scrunch_stuff.h"
#include	"../utils/error.h"

int		xfloatcmp();

void		scrunch_stream(fitsheader *fitsin, fitsheader *fitsout, int mode)
{
	int		i, j, k, N1in, N2in, N1out, N2out, imin, imax, jmin, jmax;
	int		ngood;
	float 		*fin[2], *fout, f[4];
	double		fsum;

	N1in = fitsin->n[0];
	N2in = fitsin->n[1];
	N1out = N1in / 2;
	N2out = N2in / 2;

	fin[0] = (float *) calloc(N1in, sizeof(float));
	fin[1] = (float *) calloc(N1in, sizeof(float));
	fout = (float *) calloc(N1out, sizeof(float));
	for (i = 0; i < N2in; i += 2) {
		readfitsline(fin[0], fitsin);
		readfitsline(fin[1], fitsin);
		for (j = 0; j < N1out; j++) {
			ngood = 0;
			if (fin[0][2 * j] != FLOAT_MAGIC)
				f[ngood++] = fin[0][2 * j];
			if (fin[0][2 * j + 1] != FLOAT_MAGIC)
				f[ngood++] = fin[0][2 * j + 1];
			if (fin[1][2 * j] != FLOAT_MAGIC)
				f[ngood++] = fin[1][2 * j];
			if (fin[1][2 * j + 1] != FLOAT_MAGIC)
				f[ngood++] = fin[1][2 * j + 1];
			switch (mode) {
				case	MEAN:
					if (ngood > 0) {
						fsum = 0;
						for (k = 0; k < ngood; k++)
							fsum += f[k];
						fout[j] = fsum / ngood;
					} else	
						fout[j] = FLOAT_MAGIC;
					break;
				case	CMEAN:
					if (ngood == 4) {
						fout[j] = (f[0] + f[1] + f[2] + f[3]) / ngood;
					} else {
						fout[j] = FLOAT_MAGIC;
					}
					break;
				case	MEDIAN:
					switch (ngood) {
						case 0:
							fout[j] = FLOAT_MAGIC;
							break;
						case 1:
							fout[j] = f[0];
							break;
						case 2: 
							fout[j] = 0.5 * (f[0] + f[1]);
							break;
						case 3:
							qsort(f, 3, sizeof(float), xfloatcmp);
							fout[j] = f[1];
							break;	
						case 4:
							qsort(f, 4, sizeof(float), xfloatcmp);
							fout[j] = 0.5 * (f[1] + f[2]);
							break;
						default:
							error_exit("scrunch_stream: this cannot happen");
						break;
					}
					break;
				default:
					error_exit("scrunch_stream: this cannot happen");
					break;
			}
		}
		writefitsline(fout, fitsout);
	}
	free(fin[0]); free(fin[1]), free(fout);
}



int		xfloatcmp(float *s1, float *s2) 
{
	if (*s1 < *s2)
		return (-1);
	else if (*s1 > *s2)
		return (1);
	else
		return (0);
}

