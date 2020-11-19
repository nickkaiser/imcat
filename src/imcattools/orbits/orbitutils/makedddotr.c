/*
 * makedddotr.c - definition of makedddotr()
 */

#include <math.h>
#include <stdio.h>
#include "vectors.h"

void	makedddotr(double *dddotr, double *r, double *v)
{
	int 	i;
	double	rr, rrr, rdotv;
	
	rr = dot(r, r);
	rrr = pow(rr, 1.5);
	rdotv = dot(r, v);
	for (i = 0; i < 3; i++) {
		dddotr[i] = (3 * rdotv * r[i] / rr - v[i]) / rrr;
	}
}
