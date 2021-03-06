/*
 * wrapper for Jim Heasley's routines
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kepler.h"
#include "../orbitutils/vectors.h"

#define EMU 	0.0002959122370203
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
	double	t = 0.0, emu = EMU;
	static double	*elem = NULL, *rp, *vp;

	if (!elem) {
		elem = (double *) calloc(6, sizeof(double));
		rp = (double *) calloc(3, sizeof(double));
		vp = (double *) calloc(3, sizeof(double));
	}

	/* make copies of r, v so we don't modify them */
	copy(rp, r);
	copy(vp, v);

	/* scale the velocity from AU per dynamical time to AU per day */
	scale(vp, K);

	/* call Jim's routine */
	if (ilmnts_(elem, &emu, &t, rp, vp)) {
		fprintf(stderr, "cartesiantokepler : warning : Heasley routine elmnts_() failed\n");
		return(1);
	}

	/* pack elem[] into theorbit */
	theorbit->a = elem[0];
	theorbit->e = elem[1];
	theorbit->i = elem[2];
	theorbit->Omega = elem[3];
	theorbit->omega = elem[4];
	theorbit->M = - sqrt(EMU / theorbit->a) * elem[5] / theorbit->a;

	return(0);
}



int	keplertocartesian(keplerorbit *theorbit, double *r, double *v)
{
	double	t = 0.0, emu = EMU;
	static double *elem = NULL;

	if (!elem) {
		elem = (double *) calloc(6, sizeof(double));
	}

	/* pack theorbit into elem[] */
	elem[0] = theorbit->a;
	elem[1] = theorbit->e;
	elem[2] = theorbit->i;
	elem[3] = theorbit->Omega;
	elem[4] = theorbit->omega;
	elem[5] = - theorbit->M * theorbit->a / sqrt(EMU / theorbit->a);

	/* call Jim's routine */
	if (irrdot_(elem, &emu, &t, r, v)) {
		fprintf(stderr, "keplertocartesian : warning : Heasley routine rrdot_() failed\n");
		return(1);
	}

	/* scale the velocity from AU per day to AU per dynamical time */
	scale(v, 1.0 / K);

	return(0);
}

