/*
 * op_math.c
 */

#include <stdio.h>
#include <math.h>


double	times(double x, double y)
{
	return (x * y);
}


double	plus(double x, double y)
{
	return (x + y);
}


double	divide(double x, double y)
{
	if (y == 0.0) {
		return (0.0);
	} else {
		return (x / y);
	}
}


double	minus(double x, double y)
{
	return (x - y);
}

double	max(double x, double y)
{
	return (x > y ? x : y);
}

double	min(double x, double y)
{
	return (x < y ? x : y);
}