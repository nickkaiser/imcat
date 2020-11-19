/*
 * combine_stuff.c - image combining routines 
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


#include "../utils/error.h"
#include "../utils/stats_stuff.h"
#include "combine_stuff.h"
#include "../imlib/fits.h"


int		thefloatcmp();

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

float   mean(float *f, int nel)
{
        float   result = 0.0;
	int	i;

        if (!nel)
                error_exit("mean: zero length array\n");
	
	for (i = 0; i < nel; i++)
		result += (float) f[i];
	result /= nel;
	
        return (result);
}


float   avsigclip(float *f, int nel, float clip)
{
        float   result = 0.0;
	float	uquart, lquart, median;
        int     i, lowi, medi, uppi, nused = 0;
	float	sigma;

        if (!nel)
                error_exit("avsigclip: zero length array\n");
	lowi = floor(0.5 + 0.25 * nel);
	medi = floor(0.5 + 0.50 * nel);
	uppi = floor(0.5 + 0.75 * nel);
        qsort(f, nel, sizeof(float), thefloatcmp);
	lquart = f[lowi];
	median = f[medi];
	uquart = f[uppi];

	/* now estimate sigma */
	sigma = (float) (uquart - lquart) / 1.36;

	/* now calculate mean of points within clip * sigma of median */
	for (i = 0; i < nel; i++) {
		if (f[i] > median - clip * sigma && f[i] < median + clip * sigma) {
			nused++;
			result += f[i];
		}
	}
        return ((nused ? result / nused : FLOAT_MAGIC));
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



float   	avsigclip2(float *f, float *sigma, int nel, float clip, float frac1, float frac2, float *sigmaout)
{
	float	fmedian, fsum, wtsum, f1, f2, df;
	int	i;
	static	float	*fcopy = NULL;
	static	lastnel = 0;

	/* allocate fcopy if necessary */
	if (nel > lastnel) {
		lastnel = nel;
		if (fcopy) {
			free(fcopy);
		}
		fcopy = (float *) calloc(nel, sizeof(float));
	}

	/* make copy to sort */
	for (i = 0; i < nel; i++) {
		fcopy[i] = f[i];
	}
	if (clip > 0.0) {
		if (nel < 3) {
			*sigmaout = FLOAT_MAGIC;
			return (FLOAT_MAGIC);
		}
		fmedian = median(fcopy, nel);
		f1 = frac1 * fabs(fmedian);
		f2 = frac2 * fabs(fmedian);
	}
	fsum = wtsum = 0.0;
	for (i = 0; i < nel; i++) {
		df = f[i] - fmedian;
		if (clip > 0.0) {
			if ((fabs(df) > clip * sigma[i]) && (df < f1 || df > f2)) {
				continue;
			}
		}
		fsum += f[i] / (sigma[i] * sigma[i]);
		wtsum += 1.0 / (sigma[i] * sigma[i]);
	}
	if (wtsum > 0.0) {
		*sigmaout = 1.0 / sqrt(wtsum);
		return ((float) (fsum / wtsum));
	} else {
		*sigmaout = FLOAT_MAGIC;
		return (FLOAT_MAGIC);
	}
}


float   	getsigma(float *sigmavec, int nel)
{
	float	wtsum;
	int	i;

	wtsum = 0.0;
	for (i = 0; i < nel; i++) {
		wtsum += 1.0 / (sigmavec[i] * sigmavec[i]);
	}
	if (wtsum > 0.0) {
		return(1.0 / sqrt(wtsum));
	} else {
		return (FLOAT_MAGIC);
	}
}
