/*
 * Ffunc.c - definitions of F(d, mu), and its derivative
 */

#include <math.h>
#include "Ffunc.h"


double	F(double d, double mu)
{
	double	tmp;

	tmp = 1 + 2 * d * mu + d * d;
	return ((1 - pow(tmp, -1.5)) / d);
}

double	dFdd(double d, double mu)
{
	double	tmp;

	tmp = 1 + 2 * d * mu + d * d;
	return (3 * pow(tmp, -2.5) * (mu + d) / d - (1 - pow(tmp, -1.5)) / (d * d));
}
