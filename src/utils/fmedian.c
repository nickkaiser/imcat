/* function to return the median of an array */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "fmedian.h"
#include "error.h"


float	fmedian(float *f, int nel)
{
	float	result;

	if (!nel)
		error_exit("fmedian: zero length array\n");
	qsort(f, nel, sizeof(float), descfloatcmp);
	if (2 * (nel / 2) == nel) 	/* nel is even */
		result = 0.5 * (f[nel/2] + f[nel/2 - 1]);
	else 				/* nstars is odd */
		result = f[nel/2];
	return (result);
}


int		descfloatcmp(void *p1, void *p2)
{
	float	*f1, *f2;

	f1 = (float *) p1;
	f2 = (float *) p2;
	return (*f1 < *f2 ? 1 : (*f1 == *f2 ? 0 : -1));	/* descending order */
}

