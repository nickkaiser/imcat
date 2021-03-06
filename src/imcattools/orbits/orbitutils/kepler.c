/*
 * wrapper for Jim Heasley's routines
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kepler.h"
#include "../orbitutils/vectors.h"

#define EMU     0.0002959122370203
#define K	0.01720209895

void 	assignkeplerorbitfrombuffer(keplerorbit *orbit, double *buffer)
{
	orbit->a = buffer[0];
	orbit->e = buffer[1];
	orbit->i = M_PI * buffer[2] / 180.0;
	orbit->omega = M_PI * buffer[3] / 180.0;
	orbit->Omega = M_PI * buffer[4] / 180.0;
	orbit->M = M_PI * buffer[5] / 180.0;
}

void 	fillbufferfromkeplerorbit(keplerorbit *orbit, double *buffer)
{
	buffer[0] = orbit->a;
	buffer[1] = orbit->e;
	buffer[2] = 180.0 * orbit->i / M_PI;
	buffer[3] = 180.0 * orbit->omega / M_PI;
	buffer[4] = 180.0 * orbit->Omega / M_PI;
	buffer[5] = 180.0 * orbit->M / M_PI;
}

int	cartesiantokepler(double *r, double *v, keplerorbit *theorbit)
{
	double	t = 0.0, enne, gm = K * K;
	static double	*elemin = NULL, *elemout = NULL;
	static	char	intype[4] = "CAR", outtype[4] = "KEP";

	if (!elemin) {
		elemin = (double *) calloc(6, sizeof(double));
		elemout = (double *) calloc(6, sizeof(double));
	}

	/* copy r, v so we don't modify them */
	copy(elemin, r);
	copy(elemin + 3, v);

	/* scale the velocity from AU per dynamical time to AU per day */
	scale(elemin + 3, K);

	/* call orbfit routine */
	coocha_(elemin, intype, &gm, elemout, outtype, &enne);

	/* pack elem[] into theorbit */
	theorbit->a = elemout[0];
	theorbit->e = elemout[1];
	theorbit->i = elemout[2];
	theorbit->Omega = elemout[3];
	theorbit->omega = elemout[4];
	theorbit->M = elemout[5];

	return(0);
}



int	keplertocartesian(keplerorbit *theorbit, double *r, double *v)
{
	static double	*elemin = NULL, *elemout = NULL, gm = K * K, enne;
	static	char	intype[4] = "KEP", outtype[4] = "CAR";

	if (!elemin) {
		elemin = (double *) calloc(6, sizeof(double));
		elemout = (double *) calloc(6, sizeof(double));
	}

	/* pack theorbit into elem[] */
	elemin[0] = theorbit->a;
	elemin[1] = theorbit->e;
	elemin[2] = theorbit->i;
	elemin[3] = theorbit->Omega;
	elemin[4] = theorbit->omega;
	elemin[5] = theorbit->M;

	/* call orbfit routine */
	coocha_(elemin, intype, &gm, elemout, outtype, &enne);

	copy(r, elemout);
	copy(v, elemout + 3);

	/* scale the velocity from AU per day to AU per dynamical time */
	scale(v, 1.0 / K);

	return(0);
}

