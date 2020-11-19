/*
 * avgpixstack.h
 *
 * definition of avgpixstack() function
 *
 * Nick Kaiser
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "avgpixstack.h"
#include "utils/error.h"
#include "imlib/fits.h"

int		thefloatcmp();
float   	median(float *f, int nel);
float   	rankedimage(float *f, int nel, int rank);

int	avgpixstack(int opmode, int rank, int nplanes, float *f, float *weight, float *ftemp, float *favg, float *weightsum)
{
	int	i, j;
	float	fsum, wsum;

	if (opmode < 2) {
		fsum = wsum = 0;
		for (i = 0; i < nplanes; i++) {
			if (weight[i] > 0.0) {
				if (opmode == 1) {
					fsum += f[i] * weight[i];
					wsum += weight[i];
				} else {
					fsum += f[i];
					wsum += 1.0;
				}
			}
		}
		if (wsum > 0.0) {
			*favg = fsum / wsum;
		} else {
			*favg = FLOAT_MAGIC;
		}
		*weightsum = wsum;
		return (wsum > 0.0);
	} else {
		/* copy f values with non-zero weight to ftemp */
		j = 0;
		for (i = 0; i < nplanes; i++) {
			if (weight[i] > 0.0) {
				ftemp[j++] = f[i];
			}
		}
		if (j) {
			switch(opmode) {
				case APS_MEDIAN:
					*favg = median(ftemp, j);
					*weightsum = 1.0;
					return(1);
				case APS_RANK:
					*favg = rankedimage(ftemp, j, rank);
					*weightsum = 1.0;
					return(1);
			}
		} else {
			*favg = FLOAT_MAGIC;
			*weightsum = 0.0;
			return(0);
		}
	}
}


float   median(float *f, int nel)
{
        float   result;
        int     thefloatcmp();

        if (!nel)
                error_exit("median: zero length array\n");
        qsort(f, nel, sizeof(float), thefloatcmp);
        if (2 * (nel / 2) == nel)       /* nel is even */
                result = 0.5 * (f[nel/2] + f[nel/2 - 1]);
        else                            /* nel is odd */
                result = f[nel/2];
        return (result);
}


float   	rankedimage(float *f, int nel, int rank)
{
        float   result;
        int     thefloatcmp();

        if (!nel)
                error_exit("rankedimage: zero length array\n");
        if (!rank)
                error_exit("rankedimage: rank must be non-zero\n");
        qsort(f, nel, sizeof(float), thefloatcmp);
	if (rank > 0) {
		rank--;
		return (rank < nel ? f[rank] : f[nel - 1]);
	} else {
		return(nel + rank >= 0 ? f[nel + rank] : f[0]);
	}
}


int		thefloatcmp(float *s1, float *s2) 
{
	if (*s1 < *s2)
		return (-1);
	else if (*s1 > *s2)
		return (1);
	else
		return (0);
}
