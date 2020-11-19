#include "stdio.h"
#include "tonry3d.h"
#include "../../imlib/fits.h"
#include "../../utils/error.h"
#include "psutils.h"

int	surfaceplot(fitsheader *fits, float width, float height, float alt, float az, float zfac, float zoff)
{
	float	**f;

	allocFloatArray(&f, fits->n[0], fits->n[1]);
	readfitsplane((void **) f, fits);
	mgoplt3d_(f[0], &(fits->n[0]), &(fits->n[1]), &alt, &az, &zfac, &zoff, width, height);
}


int	mgoline_(float *x1, float *y1, float *x2, float *y2)
{
	psline((double) *x1, (double) *y1, (double) *x2, (double) *y2);
}